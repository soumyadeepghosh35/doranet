#!/usr/bin/env python3
"""
Generate a DORAnet enzymatic reaction network and find pathways from the
starter(s) to every generated compound.

Usage:
    python runDORAnet.py config.yaml
"""

import csv
import sys
from pathlib import Path

import yaml
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import time

# Lazy, per-process singletons: ComponentContribution() takes 10-20s to load
# and must not be rebuilt on every reaction. Each parallel worker process
# (e.g. under pathway_ranking's multiprocessing.Pool) builds its own once.
_localCompoundCache = None
_componentContribution = None


def _getEquilibrator(compoundsDbPath):
    global _localCompoundCache, _componentContribution
    if _componentContribution is None:
        from equilibrator_api import ComponentContribution
        from equilibrator_assets.local_compound_cache import LocalCompoundCache

        _localCompoundCache = LocalCompoundCache()
        dbFile = Path(compoundsDbPath)
        if not dbFile.exists():
            # Confirmed from official docs: one-time, ~1.3 GB download from
            # Zenodo. Run this yourself ahead of time (e.g. on a login node
            # with internet access) rather than let the first real job
            # trigger it.
            _localCompoundCache.generate_local_cache_from_default_zenodo(str(dbFile))
        else:
            # NOT verified against the real API (couldn't test -- Zenodo is
            # blocked in my environment). The docs show
            # generate_local_cache_from_default_zenodo() interactively
            # prompting to overwrite an existing file, which won't work in a
            # non-interactive script -- so this line is my best guess at how
            # to point LocalCompoundCache at an already-built file instead.
            # Please confirm the correct call yourself (check
            # `help(LocalCompoundCache)` or the equilibrator_cache source)
            # before relying on this branch.
            _localCompoundCache.ccache = _localCompoundCache.ccache.__class__(str(dbFile))
        _componentContribution = ComponentContribution(ccache=_localCompoundCache.ccache)
    return _componentContribution, _localCompoundCache


def buildRxnDg(compoundsDbPath):
    """Return a rxn_thermo_calculator: dict -> standard dG' (kcal/mol), via
    eQuilibrator. Returns None (reaction rejected) if any compound can't be
    resolved/registered, or the reaction isn't atomically balanced.

    This returns a closure (compoundsDbPath is baked in from config), which
    would normally be a problem: earlier in this project, a closure-based
    calculator broke pathway_ranking's multiprocessing.Pool, since standard
    pickle can't serialize closures. It's safe here specifically because (a)
    generate_network() -- the only place this is used -- never uses
    multiprocessing at all, and (b) molecule_thermo_calculator is passed as
    None to one_step()/pathway_ranking() below, so there's nothing for that
    function's Pool to pickle either. If you ever pass this calculator into
    something that runs it via multiprocessing, it will need to become a
    top-level function again.
    """
    def rxnDg(rxn):
        return rxnDgImpl(rxn, compoundsDbPath)
    return rxnDg


def rxnDgImpl(rxn, compoundsDbPath):
    """The actual eQuilibrator calculation, called from rxnDg's closure above."""
    from equilibrator_api import Reaction

    try:
        cc, lc = _getEquilibrator(compoundsDbPath)

        allSmiles = list(rxn["reactants"]) + list(rxn["products"])
        compounds = lc.get_compounds(allSmiles)
        if any(c is None for c in compounds):
            return None
        smilesToCompound = dict(zip(allSmiles, compounds))

        stoich = {}
        for smi in rxn["reactants"]:
            cpd = smilesToCompound[smi]
            stoich[cpd] = stoich.get(cpd, 0) - 1
        for smi in rxn["products"]:
            cpd = smilesToCompound[smi]
            stoich[cpd] = stoich.get(cpd, 0) + 1

        reaction = Reaction(stoich)
        if not reaction.is_balanced():
            return None

        dgPrime = cc.standard_dg_prime(reaction)
        return dgPrime.value.m_as("kJ/mol") / 4.184  # kJ/mol -> kcal/mol, consistent with the rest of this project
    except Exception:
        return None


def moleculeProps(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "N/A", 0.0, 0
    return (
        rdMolDescriptors.CalcMolFormula(mol),
        round(rdMolDescriptors.CalcExactMolWt(mol), 4),
        mol.GetNumHeavyAtoms(),
    )


def writeMoleculesCsv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["SMILES", "Is_Starter", "MolFormula", "MolWeight", "NumHeavyAtoms"])
        for smi, isStarter in rows:
            formula, weight, heavy = moleculeProps(smi)
            writer.writerow([smi, isStarter, formula, weight, heavy])


def extractPathwayTargets(pathwaysTxtPath):
    """Return the set of molecules that had at least one successful pathway
    in *_pathways.txt: the final reaction's product side in each block."""
    pathwaysTxtPath = Path(pathwaysTxtPath)
    if not pathwaysTxtPath.exists():
        return set()

    blocks, current = [], []
    for line in pathwaysTxtPath.read_text(encoding="utf-8").splitlines():
        if line.startswith("pathway number ") and current:
            blocks.append(current)
            current = [line]
        else:
            current.append(line)
    if current:
        blocks.append(current)

    reached = set()
    for block in blocks:
        rxnLines = [line.strip() for line in block if ">>" in line]
        if rxnLines:
            reached.update(rxnLines[-1].split(">>")[-1].split("."))
    return reached


def main():
    start_time = time.time()
    if len(sys.argv) != 2:
        print("Usage: python runDORAnet.py config.yaml")
        sys.exit(1)

    config = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8")) or {}

    sys.path.insert(0, config["doranetPath"])
    import doranet.modules.enzymatic as enzymatic
    import doranet.modules.post_processing as post_processing

    jobName = config.get("jobName", "doranet_job")
    starters = set(config["starters"])
    helpers = set(config["helpers"])
    generations = int(config.get("generations", 2))
    searchDepth = int(config.get("searchDepth", generations))
    maxNumRxns = int(config.get("maxNumRxns", generations))
    minRxnAtomEconomy = float(config.get("minRxnAtomEconomy", 0.0))

    computeThermodynamics = bool(config.get("computeThermodynamics", False))
    # 0, not 15: DORAnet's own example uses a dG cutoff of 0 kcal/mol for
    # enzymatic reactions -- keep the reaction only if it's exergonic as
    # written. This is a different quantity and convention than the enthalpy
    # cutoff (15 kcal/mol) used for chemical/synthetic-module reactions.
    maxRxnThermoChange = float(config.get("maxRxnThermoChange", 0))
    compoundsDbPath = config.get("compoundsDbPath", "compounds.sqlite")
    activeRxnDg = buildRxnDg(compoundsDbPath) if computeThermodynamics else None

    print(f"Starters: {starters}")
    print(f"Generated network length: {generations}")
    if computeThermodynamics:
        print(f"Thermodynamic filtering: ON (eQuilibrator dG', cutoff {maxRxnThermoChange} kcal/mol)")
    else:
        print("Thermodynamic filtering: OFF (No_Thermo, unfiltered)")

    network = enzymatic.generate_network(
        job_name=jobName,
        starters=starters,
        gen=generations,
        max_atoms=config["maxAtoms"],
        direction="forward",
        ruleset=config.get("ruleset", "JN3604IMT"),
        rxn_thermo_calculator=activeRxnDg,
        max_rxn_thermo_change=maxRxnThermoChange,
    )
    print(f"Network: {len(network.mols)} molecules, {len(network.rxns)} reactions")

    smilesList = list(starters) + [m.uid for m in network.mols if m.uid not in starters]
    allTargets = set(smilesList) - starters - helpers
    writeMoleculesCsv(f"{jobName}_molecules_all.csv", [(s, s in starters) for s in smilesList])

    if not allTargets:
        print("No candidate molecules generated -- nothing to search pathways for.")
        return

    # DORAnet's own default weights, from post_processing.one_step. Merged
    # with any user override so a partial dict (e.g. just {by_product_number: 0}
    # to work around the Byproduct_index crash) doesn't KeyError on the keys
    # it didn't specify -- one_step() itself only fills in this whole dict
    # when weights is None, it doesn't merge partial overrides.
    defaultWeights = {
        "reaction_thermo": 2,
        "number_of_steps": 4,
        "by_product_number": 2,
        "atom_economy": 1,
        "salt_score": 0,
        "in_reaxys": 0,
        "coolness": 0,
    }
    weights = {**defaultWeights, **(config.get("weights") or {})}

    post_processing.one_step(
        networks={network},
        total_generations=generations,
        starters=starters,
        helpers=helpers,
        target=allTargets,
        job_name=jobName,
        search_depth=searchDepth,
        max_num_rxns=maxNumRxns,
        min_rxn_atom_economy=minRxnAtomEconomy,
        # No per-molecule calculator: eQuilibrator's dG is reaction-level only.
        # molecule_thermo_calculator here only ever fed the (off-by-default)
        # enol-transform correction and part of by-product scoring -- leaving
        # it None is a no-op unless you've set transformEnolsFlag elsewhere.
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=maxRxnThermoChange,
        weights=weights,
    )

    # When computeThermodynamics is True, the network was already dG-filtered
    # at generation time, so any target reached here is reachable via a
    # pathway made entirely of exergonic (or below-cutoff) steps.
    reachedTargets = extractPathwayTargets(f"{jobName}_pathways.txt") & allTargets
    writeMoleculesCsv(
        f"{jobName}_molecules_with_pathways.csv",
        [(s, True) for s in starters] + [(s, False) for s in sorted(reachedTargets)],
    )
    label = "thermodynamically-supported (eQuilibrator dG')" if computeThermodynamics else "generated"
    print(f"Compounds with a {label} pathway: {len(reachedTargets)} / {len(allTargets)}")
    print(f"Total time: {time.time() - start_time:.2f}s")


if __name__ == "__main__":
    main()