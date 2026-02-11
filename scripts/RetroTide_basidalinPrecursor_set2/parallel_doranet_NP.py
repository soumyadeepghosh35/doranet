#!/usr/bin/env python

import sys
import os
import argparse
import time
import csv
from pathlib import Path
from multiprocessing import Process, Queue, cpu_count
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
from rdkit import RDLogger
import yaml

RDLogger.DisableLog('rdApp.*')

DORANET_PATH = Path("/users/sghosh6/DTRA_project/MACAW/doranet")
sys.path.insert(0, str(DORANET_PATH))

import doranet.modules.enzymatic as enzymatic
import doranet.modules.post_processing as post_processing

HELPERS = {
    'O', 'O=O', '[H][H]', 'O=C=O', 'C=O', '[C-]#[O+]', 'Br', '[Br][Br]',
    'CO', 'C=C', 'O=S(O)O', 'N', 'O=S(=O)(O)O', 'O=NO', 'N#N',
    'O=[N+]([O-])O', 'NO', 'C#N', 'S', 'O=S=O', 'N#CO', '[H+]', 'OO',
    'Cl', 'I', 'O=C(O)O', 'O=P(O)(O)O', 'O=P(O)(O)OP(=O)(O)O', 'C',
    'CC', 'CC=O', 'CC(=O)O', 'CCC(=O)O'
}

# =========================================================================
# Natural product databases (hardcoded, not in config)
# =========================================================================
NP_DATABASES = [
    ('../natural_products_data/NPAtlas_np_data.sdf', 'npatlas'),
    ('../natural_products_data/coconut_np_data.sdf', 'coconut')
]


def loadNaturalProductInchiKeys():
    """Load all natural product InChIKeys from SDF databases."""
    npInchiKeys = set()
    npSmilesMap = {}  # InChIKey -> canonical SMILES for reverse lookup

    print("\n" + "=" * 70)
    print("LOADING NATURAL PRODUCT DATABASES")
    print("=" * 70)

    for sdfPath, dbName in NP_DATABASES:
        print(f"Loading {dbName} from: {sdfPath}")

        if not os.path.exists(sdfPath):
            print(f"  WARNING: File not found: {sdfPath}, skipping...")
            continue

        supplier = Chem.SDMolSupplier(sdfPath)
        validCount = 0
        invalidCount = 0

        for mol in supplier:
            if mol is not None:
                try:
                    canonicalSmiles = Chem.MolToSmiles(mol)
                    inchi = MolToInchi(mol)
                    if inchi is not None:
                        inchiKey = InchiToInchiKey(inchi)
                        npInchiKeys.add(inchiKey)
                        npSmilesMap[inchiKey] = canonicalSmiles
                        validCount += 1
                    else:
                        invalidCount += 1
                except:
                    invalidCount += 1
            else:
                invalidCount += 1

        print(f"  {dbName}: {validCount} valid molecules, {invalidCount} invalid/skipped")

    print(f"\nTotal unique NP InChIKeys loaded: {len(npInchiKeys)}")
    print("=" * 70)

    return npInchiKeys, npSmilesMap


def smilesToInchiKey(smiles):
    """Convert a SMILES string to an InChIKey."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        try:
            inchi = MolToInchi(mol)
            if inchi is not None:
                return InchiToInchiKey(inchi)
        except:
            pass
    return None


def findNPTargets(smilesList, starterSmiles, npInchiKeys):
    """
    From generated molecules, find those that match natural products
    (excluding starters and helpers) using InChIKey matching.
    Returns a set of SMILES that are natural product matches.
    """
    npTargets = set()
    userStarters = {starterSmiles}

    for smi in smilesList:
        if smi in userStarters or smi in HELPERS:
            continue

        inchiKey = smilesToInchiKey(smi)
        if inchiKey is not None and inchiKey in npInchiKeys:
            npTargets.add(smi)

    return npTargets


def loadConfig(configPath):
    with open(configPath, 'r') as f:
        config = yaml.safe_load(f)

    defaults = {
        'input': {
            'SMILESfile': None,
            'smiles_column': 'Canonical_SMILES',
            'start_index': 0,
            'num_smiles': None
        },
        'output': {
            'directory': 'doranet_output'
        },
        'parallel': {
            'num_workers': 4
        },
        'network': {
            'generations': 2,
            'ruleset': 'JN3604IMT',
            'direction': 'forward'
        },
        'max_atoms': {
            'C': 12,
            'N': 3,
            'O': 5,
            'S': 3
        }
    }

    for section, values in defaults.items():
        if section not in config:
            config[section] = values
        elif isinstance(values, dict):
            for key, defaultVal in values.items():
                if key not in config[section]:
                    config[section][key] = defaultVal

    return config


def validateConfig(config):
    errors = []

    if not config['input']['SMILESfile']:
        errors.append("input.SMILESfile is required")
    elif not os.path.exists(config['input']['SMILESfile']):
        errors.append(f"Input file not found: {config['input']['SMILESfile']}")

    if config['parallel']['num_workers'] < 1:
        errors.append("parallel.num_workers must be at least 1")

    if config['network']['generations'] < 1:
        errors.append("network.generations must be at least 1")

    if errors:
        print("Configuration errors:")
        for err in errors:
            print(f"  - {err}")
        sys.exit(1)

    config['parallel']['num_workers'] = min(
        config['parallel']['num_workers'],
        cpu_count()
    )

    return config


def printConfig(config):
    print("=" * 70)
    print("CONFIGURATION (Natural Product Target Mode)")
    print("=" * 70)
    print(f"Input file: {config['input']['SMILESfile']}")
    print(f"SMILES column: {config['input']['smiles_column']}")
    print(f"Start index: {config['input']['start_index']}")
    print(f"Number to process: {config['input']['num_smiles'] or 'all'}")
    print(f"Output directory: {config['output']['directory']}_NP")
    print(f"Workers: {config['parallel']['num_workers']}")
    print(f"Generations: {config['network']['generations']}")
    print(f"Ruleset: {config['network']['ruleset']}")
    print(f"Direction: {config['network']['direction']}")
    print(f"Max atoms: {config['max_atoms']}")
    print(f"NP Databases: {[db[1] for db in NP_DATABASES]}")
    print("=" * 70)


def readSmilesFromCsv(csvPath, smilesColumn, startIdx, numSmiles):
    smilesList = []

    with open(csvPath, 'r') as f:
        reader = csv.DictReader(f)

        if smilesColumn not in reader.fieldnames:
            available = ', '.join(reader.fieldnames)
            raise ValueError(f"Column '{smilesColumn}' not found. Available: {available}")

        for idx, row in enumerate(reader):
            if idx < startIdx:
                continue
            if numSmiles is not None and len(smilesList) >= numSmiles:
                break

            smi = row[smilesColumn]
            if smi and smi.strip():
                smilesList.append(smi.strip())

    return smilesList


def processSingleStarter(starterIdx, starterSmiles, config, npInchiKeys):
    """
    Process a single starter molecule.
    Uses natural product database matches as targets instead of all generated molecules.
    """
    jobName = f"starter_{starterIdx:05d}_NP"
    jobOutputDir = Path(config['output']['directory'] + "_NP") / jobName
    jobOutputDir.mkdir(parents=True, exist_ok=True)

    originalDir = os.getcwd()
    os.chdir(jobOutputDir)

    result = {
        'starter_idx': starterIdx,
        'starter_smiles': starterSmiles,
        'status': 'failed',
        'num_molecules': 0,
        'num_np_targets': 0,
        'num_total_generated': 0,
        'time_seconds': 0,
        'error': None
    }

    startTime = time.time()

    try:
        userStarters = {starterSmiles}

        forwardNetwork = enzymatic.generate_network(
            job_name=jobName,
            starters=userStarters,
            gen=config['network']['generations'],
            max_atoms=config['max_atoms'],
            direction=config['network']['direction'],
            ruleset=config['network']['ruleset']
        )

        smilesList = [mol.uid for mol in forwardNetwork.mols]
        result['num_molecules'] = len(smilesList)

        # Find all generated molecules (excluding starters and helpers)
        allGenerated = set(smilesList) - userStarters - HELPERS
        result['num_total_generated'] = len(allGenerated)

        # Find NP targets using InChIKey matching
        npTargets = findNPTargets(smilesList, starterSmiles, npInchiKeys)
        result['num_np_targets'] = len(npTargets)

        # Save molecules CSV with NP flag
        outputPath = Path(f"{jobName}_molecules.csv")
        with open(outputPath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'SMILES', 'Is_Starter', 'Is_Helper', 'IsNaturalProduct',
                'MolFormula', 'MolWeight', 'NumHeavyAtoms', 'InChIKey'
            ])
            for smi in smilesList:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    formula = rdMolDescriptors.CalcMolFormula(mol)
                    molWeight = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
                    numHeavy = mol.GetNumHeavyAtoms()
                    inchiKey = smilesToInchiKey(smi) or "N/A"
                else:
                    formula, molWeight, numHeavy, inchiKey = "N/A", 0, 0, "N/A"

                isNP = smi in npTargets
                isHelper = smi in HELPERS
                isStarter = smi in userStarters

                writer.writerow([
                    smi, isStarter, isHelper, isNP,
                    formula, molWeight, numHeavy, inchiKey
                ])

        # Run post-processing only on NP targets
        if npTargets:
            post_processing.one_step(
                networks={forwardNetwork},
                total_generations=config['network']['generations'],
                starters=userStarters,
                helpers=HELPERS,
                target=npTargets,
                job_name=jobName,
            )

        result['status'] = 'success'

    except Exception as e:
        result['status'] = 'failed'
        result['error'] = str(e)

    finally:
        os.chdir(originalDir)

    result['time_seconds'] = round(time.time() - startTime, 2)

    return result


def workerProcess(workerId, taskQueue, resultQueue, config, totalTasks, npInchiKeys):
    """
    Worker process that continuously picks up tasks from the queue.
    Now receives npInchiKeys for NP matching.
    """
    while True:
        task = taskQueue.get()

        if task is None:
            break

        starterIdx, starterSmiles = task

        result = processSingleStarter(starterIdx, starterSmiles, config, npInchiKeys)

        displayIdx = starterIdx + 1
        displayWorkerId = workerId + 1

        print(f"[Worker {displayWorkerId}] [{displayIdx:05d}/{totalTasks:05d}] "
              f"{result['status'].upper()} | "
              f"Molecules: {result['num_molecules']} | "
              f"NP Targets: {result['num_np_targets']}/{result['num_total_generated']} | "
              f"Time: {result['time_seconds']}s | "
              f"SMILES: {starterSmiles[:40]}...")

        resultQueue.put(result)


def runParallelProcessing(smilesList, config):
    numWorkers = config['parallel']['num_workers']
    totalTasks = len(smilesList)
    
    print(f"\nStarting parallel processing")
    print(f"  Total SMILES: {totalTasks}")
    print(f"  Workers: {numWorkers}")
    print("-" * 70)
    
    # Create task queue and populate with all tasks
    taskQueue = Queue()
    for idx, smi in enumerate(smilesList):
        taskQueue.put((idx, smi))
    
    # Add poison pills (one per worker) to signal completion
    for _ in range(numWorkers):
        taskQueue.put(None)
    
    # Create result queue
    resultQueue = Queue()
    
    # Create and start worker processes
    workers = []
    for workerId in range(numWorkers):
        p = Process(
            target=workerProcess,
            args=(workerId, taskQueue, resultQueue, config, totalTasks)
        )
        p.daemon = False
        workers.append(p)
        p.start()
    
    # Collect results as they come in (don't wait for join first)
    results = []
    for _ in range(totalTasks):
        try:
            result = resultQueue.get(timeout=600)  # 10 min timeout per result
            results.append(result)
        except Exception as e:
            print(f"Warning: Timeout or error collecting result: {e}")
            break
    
    # Now wait for all workers to finish
    for p in workers:
        p.join(timeout=60)
        if p.is_alive():
            print(f"Warning: Worker {p.pid} did not exit cleanly, terminating...")
            p.terminate()
            p.join(timeout=10)
    
    # Sort results by starter index
    results.sort(key=lambda x: x['starter_idx'])
    
    print(f"\nCollected {len(results)} / {totalTasks} results")
    
    return results


def saveSummary(results, outputDir):
    summaryPath = Path(outputDir) / "processing_summary_NP.csv"

    with open(summaryPath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'starter_idx', 'starter_smiles', 'status', 'num_molecules',
            'num_total_generated', 'num_np_targets', 'time_seconds', 'error'
        ])

        for r in results:
            writer.writerow([
                r['starter_idx'],
                r['starter_smiles'],
                r['status'],
                r['num_molecules'],
                r['num_total_generated'],
                r['num_np_targets'],
                r['time_seconds'],
                r['error'] or ''
            ])

    return summaryPath


def printSummary(results, totalTime):
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']

    totalMolecules = sum(r['num_molecules'] for r in successful)
    totalGenerated = sum(r['num_total_generated'] for r in successful)
    totalNPTargets = sum(r['num_np_targets'] for r in successful)
    avgTime = sum(r['time_seconds'] for r in results) / len(results) if results else 0
    startersWithNP = sum(1 for r in successful if r['num_np_targets'] > 0)

    npPercent = (totalNPTargets / totalGenerated * 100) if totalGenerated > 0 else 0

    print("\n" + "=" * 70)
    print("PROCESSING SUMMARY (Natural Product Target Mode)")
    print("=" * 70)
    print(f"Total starters processed: {len(results)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")
    print(f"Total molecules in networks: {totalMolecules}")
    print(f"Total generated (excl. starters/helpers): {totalGenerated}")
    print(f"Total NP targets found: {totalNPTargets} ({npPercent:.2f}%)")
    print(f"Starters with ≥1 NP target: {startersWithNP}/{len(successful)}")
    print(f"Average time per starter: {avgTime:.2f}s")
    print(f"Total wall time: {totalTime:.2f}s")

    if failed:
        print(f"\nFailed starters:")
        for r in failed[:10]:
            print(f"  [{r['starter_idx']}] {r['starter_smiles'][:50]}... | Error: {r['error']}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more")

    print("=" * 70)


def saveConfigCopy(config, outputDir):
    configCopyPath = Path(outputDir) / "config_used_NP.yaml"
    with open(configCopyPath, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    return configCopyPath


def main():
    parser = argparse.ArgumentParser(
        description='Parallel DORAnet pathway generation with Natural Product targets'
    )
    parser.add_argument(
        '-c', '--config',
        default='config.yaml',
        help='Path to configuration YAML file (default: config.yaml)'
    )
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: Configuration file '{args.config}' not found")
        sys.exit(1)

    print(f"Loading configuration from: {args.config}")
    config = loadConfig(args.config)
    config = validateConfig(config)
    printConfig(config)

    # Append _NP to output directory
    outputDir = Path(config['output']['directory'] + "_NP")
    outputDir.mkdir(parents=True, exist_ok=True)

    configCopy = saveConfigCopy(config, outputDir)
    print(f"Configuration saved to: {configCopy}")

    # Load natural product databases
    npInchiKeys, npSmilesMap = loadNaturalProductInchiKeys()

    if not npInchiKeys:
        print("Error: No natural product InChIKeys loaded. Check database paths.")
        sys.exit(1)

    print("\nReading SMILES from CSV...")
    smilesList = readSmilesFromCsv(
        config['input']['SMILESfile'],
        config['input']['smiles_column'],
        config['input']['start_index'],
        config['input']['num_smiles']
    )

    if not smilesList:
        print("Error: No valid SMILES found in input file")
        sys.exit(1)

    print(f"Loaded {len(smilesList)} SMILES")

    totalStartTime = time.time()

    results = runParallelProcessing(smilesList, config, npInchiKeys)

    totalTime = time.time() - totalStartTime

    summaryPath = saveSummary(results, outputDir)
    print(f"\nSummary saved to: {summaryPath}")

    printSummary(results, totalTime)


if __name__ == "__main__":
    main()