#!/usr/bin/env python

"""
Enzymatic network expansion using a custom cofactor table.

DORAnet ships its cofactors in all_cofactors.tsv and builds its cofactor
lookup tables from that file at import. This script instead reads a custom
cofactor TSV (same #ID / Name / SMILES schema) whose path is given in
config.yaml, then rebuilds those tables so expansion uses only the cofactors
defined in that file.

Any shipped cofactor absent from the custom file is excluded automatically,
which drops every rule that would reference it and keeps the filters from
looking up a missing token.

Usage:
    python runDORAnet.py config.yaml
"""

import csv
import sys
import time
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def loadConfig(configPath):
    with open(configPath) as handle:
        return yaml.safe_load(handle)


def resolvePath(rawPath, baseDir):
    path = Path(rawPath)
    return path if path.is_absolute() else (baseDir / path)


def canonicalSmiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def loadCofactorsIntoModule(enzModule, cofactorFile, extraExcluded):
    """Rebuild the enzymatic module's cofactor tables from a custom TSV.

    enzModule : the doranet.modules.enzymatic.generate_network module.
    cofactorFile : path to a TSV with columns '#ID', 'Name', 'SMILES'.
    extraExcluded : cofactor ids present in the file that should still be
        excluded (e.g. CARBONYL_CoF / AMINO_CoF, which stock DORAnet excludes).

    The '#ID' must match a cofactor token used by the rule set for that
    cofactor to occupy a reaction slot.
    """
    shippedIds = set(enzModule.cofactors_dict.keys())

    table = pd.read_csv(cofactorFile, sep="\t")
    if not {"#ID", "SMILES"}.issubset(table.columns):
        raise ValueError(
            f"{cofactorFile} must have '#ID' and 'SMILES' columns"
        )

    customDict = {}
    customSet = set()
    for cofId, smiles in zip(table["#ID"], table["SMILES"]):
        canon = canonicalSmiles(smiles)
        if canon is None:
            print(f"WARNING invalid SMILES for cofactor '{cofId}', skipping")
            continue
        if cofId not in shippedIds:
            print(
                f"WARNING cofactor id '{cofId}' is not a token used by any "
                f"shipped rule, so it will never occupy a cofactor slot"
            )
        customDict[cofId] = canon
        customSet.add(canon)

    customIds = set(customDict.keys())
    excluded = (shippedIds - customIds) | set(extraExcluded)

    cleanDict = {}
    cleanSet = set()
    for cofId, smiles in zip(table["#ID"], table["SMILES"]):
        if cofId in customDict and cofId not in excluded:
            cleaned = enzModule.clean_SMILES(smiles)
            cleanDict[cofId] = cleaned
            cleanSet.add(cleaned)

    enzModule.cofactors_path = Path(cofactorFile)
    enzModule.cofactors_dict = customDict
    enzModule.cofactors_set = customSet
    enzModule.excluded_cofactors = tuple(excluded)
    enzModule.cofactors_clean_dict = cleanDict
    enzModule.cofactors_clean = cleanSet

    activeTokens = customIds - set(excluded)
    print(f"Cofactor file: {cofactorFile}")
    print(f"Cofactor pool in use ({len(activeTokens)} tokens): "
          f"{sorted(activeTokens)}")


def writeMoleculeCsv(smilesList, userStarters, outputPath):
    with open(outputPath, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["SMILES", "isStarter", "molFormula", "molWeight", "numHeavyAtoms"]
        )
        for smi in smilesList:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                molWeight = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
                numHeavy = mol.GetNumHeavyAtoms()
            else:
                formula, molWeight, numHeavy = "N/A", 0, 0
            writer.writerow(
                [smi, smi in userStarters, formula, molWeight, numHeavy]
            )


def main():
    configPath = Path(sys.argv[1] if len(sys.argv) > 1 else "config.yaml")
    config = loadConfig(configPath)
    baseDir = configPath.resolve().parent

    sys.path.insert(0, str(resolvePath(config["doranetPath"], baseDir)))

    import doranet.modules.enzymatic as enzymatic
    import doranet.modules.enzymatic.generate_network as enzGenMod
    import doranet.modules.post_processing as post_processing

    startTime = time.time()

    jobConfig = config["job"]
    jobName = jobConfig["name"]
    userStarters = set(config["starters"])
    userHelpers = set(config["helpers"])

    cofConfig = config["cofactors"]
    cofactorFile = resolvePath(cofConfig["file"], baseDir)
    loadCofactorsIntoModule(
        enzGenMod,
        cofactorFile,
        cofConfig.get("excludedCofactors", []),
    )

    print(f"Starter: {list(userStarters)[0]}")

    forwardNetwork = enzymatic.generate_network(
        job_name=jobName,
        starters=userStarters,
        gen=jobConfig["generations"],
        max_atoms=jobConfig.get("maxAtoms"),
        direction=jobConfig.get("direction", "forward"),
        ruleset=jobConfig.get("ruleset", "JN3604IMT"),
    )

    smilesList = list(userStarters) + [
        mol.uid for mol in forwardNetwork.mols if mol.uid not in userStarters
    ]
    print(f"\nGenerated {len(smilesList) - len(userStarters)} new molecules "
          f"+ {len(userStarters)} starters")

    allTargets = set(smilesList) - userStarters - userHelpers
    print(f"Using {len(allTargets)} generated molecules as targets")

    outputPath = Path(f"{jobName}_molecules.csv")
    writeMoleculeCsv(smilesList, userStarters, outputPath)
    print(f"Saved molecules to {outputPath}")

    if allTargets:
        post_processing.one_step(
            networks={forwardNetwork},
            total_generations=jobConfig["generations"],
            starters=userStarters,
            helpers=userHelpers,
            target=allTargets,
            job_name=jobName,
        )
        print(f"Pathway files generated with prefix: {jobName}")
    else:
        print("No targets found for pathway generation")

    print(f"Time: {time.time() - startTime:.2f} s")


if __name__ == "__main__":
    main()