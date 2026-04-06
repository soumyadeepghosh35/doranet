import os

os.environ["OMP_NUM_THREADS"]      = "1"
os.environ["MKL_NUM_THREADS"]      = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"]  = "1"

import math
import time
import argparse
import yaml
import pandas as pd
from tqdm import tqdm
from DORA_XGB import DORA_XGB


def loadConfig(configPath):
    if not os.path.exists(configPath):
        raise FileNotFoundError(f"Config file not found: {configPath}")
    with open(configPath, "r") as f:
        config = yaml.safe_load(f)
    print(f"Config loaded from: {configPath}")
    return config


def resolveWorkingPaths(inputConfig, outputConfig):
    cwd            = os.getcwd()
    fileNamePrefix = outputConfig["fileNamePrefix"]
    paths = {
        "inputCSV"      : os.path.join(cwd, inputConfig["reactionCSV"]),
        "outputFullCSV" : os.path.join(cwd, f"{fileNamePrefix}_reactionDF_DORAXGB.csv"),
        "outputHighCSV" : os.path.join(cwd, f"{fileNamePrefix}_reactionDF_DORAXGB_highFeasibility.csv"),
    }
    print(f"Working directory : {cwd}")
    print(f"Input CSV         : {paths['inputCSV']}")
    print(f"Output full CSV   : {paths['outputFullCSV']}")
    print(f"Output high CSV   : {paths['outputHighCSV']}")
    if not os.path.exists(paths["inputCSV"]):
        raise FileNotFoundError(f"Input CSV not found: {paths['inputCSV']}")
    return paths


def loadReactionDF(inputCsvPath, testMode, testRows):
    reactionDF = pd.read_csv(inputCsvPath)
    if "reactionString" not in reactionDF.columns:
        raise ValueError(
            f"'reactionString' column not found in {inputCsvPath}. "
            f"Available columns: {list(reactionDF.columns)}"
        )
    if testMode:
        reactionDF = reactionDF.iloc[:testRows].copy()
        print(f"Test mode enabled: using first {testRows} rows only")
    print(f"Reactions loaded : {len(reactionDF)}")
    print(f"Columns          : {list(reactionDF.columns)}")
    return reactionDF


def loadModels():
    models = {
        "byDescMW"   : DORA_XGB.feasibility_classifier(cofactor_positioning="by_descending_MW"),
        "byAscMW"    : DORA_XGB.feasibility_classifier(cofactor_positioning="by_ascending_MW"),
        "addConcat"  : DORA_XGB.feasibility_classifier(cofactor_positioning="add_concat"),
        "addSubtract": DORA_XGB.feasibility_classifier(cofactor_positioning="add_subtract"),
    }
    print("Models loaded successfully")
    return models


def getFeasibilityScoresAndLabels(rxnStr, models):
    return pd.Series({
        "feasibilityScore_rule1" : models["byDescMW"].predict_proba(rxnStr),
        "feasibilityLabel_rule1" : models["byDescMW"].predict_label(rxnStr),
        "feasibilityScore_rule2" : models["byAscMW"].predict_proba(rxnStr),
        "feasibilityLabel_rule2" : models["byAscMW"].predict_label(rxnStr),
        "feasibilityScore_rule3" : models["addConcat"].predict_proba(rxnStr),
        "feasibilityLabel_rule3" : models["addConcat"].predict_label(rxnStr),
        "feasibilityScore_rule4" : models["addSubtract"].predict_proba(rxnStr),
        "feasibilityLabel_rule4" : models["addSubtract"].predict_label(rxnStr),
    })


def runFeasibilityScoring(reactionDF, models, chunkSize):
    reactionDF = reactionDF.copy()
    reactionDF["rxn_str"] = (
        reactionDF["reactionString"]
        .astype(str)
        .str.replace(" ", "", regex=False)
    )

    totalRows    = len(reactionDF)
    numChunks    = math.ceil(totalRows / chunkSize)
    chunkResults = []

    print(f"Total reactions : {totalRows}")
    print(f"Chunk size      : {chunkSize}")
    print(f"Total chunks    : {numChunks}")

    overallStartTime = time.time()

    for chunkIdx in tqdm(range(numChunks), desc="Scoring chunks", unit="chunk"):
        chunkStart = chunkIdx * chunkSize
        chunkEnd   = min(chunkStart + chunkSize, totalRows)
        chunkDF    = reactionDF.iloc[chunkStart:chunkEnd].copy()

        chunkStartTime   = time.time()
        feasibilityChunk = chunkDF["rxn_str"].apply(
            lambda rxnStr: getFeasibilityScoresAndLabels(rxnStr, models)
        )
        chunkDF      = pd.concat([chunkDF, feasibilityChunk], axis=1)
        chunkResults.append(chunkDF)

        chunkElapsed      = time.time() - chunkStartTime
        overallElapsed    = time.time() - overallStartTime
        rowsRemaining     = totalRows - chunkEnd
        estimatedTimeLeft = (overallElapsed / chunkEnd) * rowsRemaining if chunkEnd > 0 else 0

        tqdm.write(
            f"Chunk {chunkIdx + 1:>4}/{numChunks} | "
            f"Rows {chunkStart:>7} - {chunkEnd:<7} | "
            f"Chunk time: {chunkElapsed:>6.1f}s | "
            f"Elapsed: {overallElapsed:>7.1f}s | "
            f"ETA: {estimatedTimeLeft:>7.1f}s"
        )

    reactionDF_DORAXGB = pd.concat(chunkResults, ignore_index=True)

    totalElapsed = time.time() - overallStartTime
    print(f"Done in {totalElapsed:.1f}s | "
          f"{totalRows / totalElapsed:.0f} reactions/sec | "
          f"Shape: {reactionDF_DORAXGB.shape}")

    return reactionDF_DORAXGB


def saveOutputs(reactionDF, paths, filteringConfig):
    rule     = filteringConfig["highFeasibilityRule"]
    label    = filteringConfig["highFeasibilityLabel"]
    scoreCol = f"feasibilityScore_{rule}"
    labelCol = f"feasibilityLabel_{rule}"

    reactionDF.to_csv(paths["outputFullCSV"], index=False, encoding="utf-8")
    print(f"Saved full scored reactions   : {paths['outputFullCSV']}")

    highFeasibilityDF = (
        reactionDF[
            (reactionDF[labelCol] == label) &
            (reactionDF[scoreCol].notna())
        ]
        .sort_values(by=scoreCol, ascending=False)
        .reset_index(drop=True)
    )
    highFeasibilityDF.to_csv(paths["outputHighCSV"], index=False, encoding="utf-8")
    print(f"Saved high-feasibility subset : {paths['outputHighCSV']}")
    print(f"High-feasibility reactions    : {len(highFeasibilityDF)} / {len(reactionDF)}")


def main():
    parser = argparse.ArgumentParser(
        description="Run DORA XGB feasibility scoring on a reaction CSV."
    )
    parser.add_argument(
        "config",
        type=str,
        help="Path to YAML config file (e.g. doraxgb_config.yaml)"
    )
    args   = parser.parse_args()
    config = loadConfig(args.config)

    scoringConfig   = config["scoring"]
    testMode        = scoringConfig.get("testMode", False)
    testRows        = scoringConfig.get("testRows", 100)
    chunkSize       = scoringConfig.get("chunkSize", 5000)

    print("---------------------------------")
    print("Step 1: File paths")
    print("---------------------------------")
    paths = resolveWorkingPaths(config["input"], config["output"])

    print("---------------------------------")
    print("Step 2: Load reactions")
    print("---------------------------------")
    reactionDF = loadReactionDF(paths["inputCSV"], testMode, testRows)

    print("---------------------------------")
    print("Step 3: Load models")
    print("---------------------------------")
    models = loadModels()

    print("---------------------------------")
    print("Step 4: Run feasibility scoring")
    print("---------------------------------")
    reactionDF_DORAXGB = runFeasibilityScoring(
        reactionDF = reactionDF,
        models     = models,
        chunkSize  = chunkSize,
    )

    print("---------------------------------")
    print("Step 5: Save outputs")
    print("---------------------------------")
    saveOutputs(
        reactionDF      = reactionDF_DORAXGB,
        paths           = paths,
        filteringConfig = config["filtering"],
    )

    print("---------------------------------")
    print(" !! Job completed !!")
    print("---------------------------------")


if __name__ == "__main__":
    main()