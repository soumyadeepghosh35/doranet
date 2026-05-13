#!/usr/bin/env python3
import os
import argparse
import time
import yaml
import pandas as pd
from tqdm import tqdm
from functools import lru_cache
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

# -------------------------------------------------------------------
# MUST be set before model/numpy-heavy imports
# -------------------------------------------------------------------
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"
# Optional hard-disable dynamic threading in some libs
os.environ["OMP_DYNAMIC"] = "FALSE"
os.environ["MKL_DYNAMIC"] = "FALSE"

from DORA_XGB import DORA_XGB  # noqa: E402


RULE_TO_POSITIONING = {
    "rule1": "by_descending_MW",
    "rule2": "by_ascending_MW",
    "rule3": "add_concat",
    "rule4": "add_subtract",
}

# Worker globals
_W_MODELS = None
_W_RULES = None
_W_THRESHOLD = 0.5


# -------------------- Affinity helpers (Linux) --------------------
def set_affinity_or_warn(core_ids):
    """
    Restrict current process to given CPU core IDs (Linux only).
    """
    if core_ids is None:
        return
    try:
        os.sched_setaffinity(0, set(core_ids))
    except AttributeError:
        print("[WARN] CPU affinity not supported on this OS. Cannot strictly pin cores.")
    except Exception as e:
        print(f"[WARN] Failed to set CPU affinity: {e}")


def get_default_core_ids(max_cores):
    avail = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else list(range(os.cpu_count() or 1))
    if max_cores > len(avail):
        raise ValueError(f"Requested max_cores={max_cores}, but only {len(avail)} CPUs available in affinity mask.")
    return avail[:max_cores]


# -------------------- Worker init / scoring --------------------
def init_worker(rules_to_run, threshold, core_ids):
    global _W_MODELS, _W_RULES, _W_THRESHOLD
    # Pin each worker to allowed cores set
    set_affinity_or_warn(core_ids)

    _W_RULES = rules_to_run
    _W_THRESHOLD = threshold
    _W_MODELS = {
        r: DORA_XGB.feasibility_classifier(cofactor_positioning=RULE_TO_POSITIONING[r])
        for r in rules_to_run
    }


def score_one_worker(rxn_str):
    out = {"rxn_str": rxn_str}
    for r in _W_RULES:
        p = _W_MODELS[r].predict_proba(rxn_str)
        out[f"feasibilityScore_{r}"] = p
        out[f"feasibilityLabel_{r}"] = int(p >= _W_THRESHOLD)
    return out


# -------------------- Core pipeline --------------------
def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def validate_rules(rules):
    if not rules:
        raise ValueError("scoring.rulesToRun is empty.")
    bad = [r for r in rules if r not in RULE_TO_POSITIONING]
    if bad:
        raise ValueError(f"Invalid rules: {bad}")


def load_reactions(path, test_mode=False, test_rows=100):
    df = pd.read_csv(path)
    if "reactionString" not in df.columns:
        raise ValueError("'reactionString' column not found.")
    if test_mode:
        df = df.iloc[:test_rows].copy()
    df["rxn_str"] = df["reactionString"].astype(str).str.replace(" ", "", regex=False)
    return df


def score_unique_single(unique_rxns, rules, threshold, cache_size):
    models = {r: DORA_XGB.feasibility_classifier(cofactor_positioning=RULE_TO_POSITIONING[r]) for r in rules}

    @lru_cache(maxsize=cache_size)
    def pred(rxn):
        d = {}
        for r in rules:
            p = models[r].predict_proba(rxn)
            d[f"feasibilityScore_{r}"] = p
            d[f"feasibilityLabel_{r}"] = int(p >= threshold)
        return d

    rows = []
    for rxn in tqdm(unique_rxns, desc="Scoring (single)", unit="rxn"):
        d = pred(rxn)
        d["rxn_str"] = rxn
        rows.append(d)
    return pd.DataFrame(rows)


def score_unique_mp(unique_rxns, rules, threshold, max_cores, mp_chunk_size, core_ids):
    # Use spawn for stricter process behavior
    ctx = mp.get_context("spawn")

    rows = []
    with ProcessPoolExecutor(
        max_workers=max_cores,
        mp_context=ctx,
        initializer=init_worker,
        initargs=(rules, threshold, core_ids),
    ) as ex:
        it = ex.map(score_one_worker, unique_rxns, chunksize=mp_chunk_size)
        for d in tqdm(it, total=len(unique_rxns), desc=f"Scoring ({max_cores} cores)", unit="rxn"):
            rows.append(d)

    return pd.DataFrame(rows)


def run_scoring(df, rules, threshold, use_mp, max_cores, mp_chunk_size, cache_size, core_ids):
    t0 = time.time()
    unique_rxns = pd.Series(df["rxn_str"].dropna().unique()).tolist()

    if use_mp:
        scored_u = score_unique_mp(
            unique_rxns=unique_rxns,
            rules=rules,
            threshold=threshold,
            max_cores=max_cores,
            mp_chunk_size=mp_chunk_size,
            core_ids=core_ids,
        )
    else:
        scored_u = score_unique_single(unique_rxns, rules, threshold, cache_size)

    out = df.merge(scored_u, on="rxn_str", how="left")
    dt = time.time() - t0
    print(f"Done in {dt:.1f}s | {len(df)/max(dt,1e-9):.1f} rows/s")
    return out


def save_outputs(df, in_cfg, out_cfg, filt_cfg, rules):
    cwd = os.getcwd()
    prefix = out_cfg["fileNamePrefix"]
    full_csv = os.path.join(cwd, f"{prefix}_reactionDF_DORAXGB.csv")
    high_csv = os.path.join(cwd, f"{prefix}_reactionDF_DORAXGB_highFeasibility.csv")

    rule = filt_cfg["highFeasibilityRule"]
    label = filt_cfg["highFeasibilityLabel"]
    if rule not in rules:
        raise ValueError(f"filtering rule {rule} not in rulesToRun {rules}")

    s_col = f"feasibilityScore_{rule}"
    l_col = f"feasibilityLabel_{rule}"

    df.to_csv(full_csv, index=False)

    high = (
        df[(df[l_col] == label) & (df[s_col].notna())]
        .sort_values(s_col, ascending=False)
        .reset_index(drop=True)
    )
    high.to_csv(high_csv, index=False)
    print(f"Saved: {full_csv}")
    print(f"Saved: {high_csv}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=str)
    args = ap.parse_args()

    cfg = load_config(args.config)
    sc = cfg.get("scoring", {})

    rules = sc.get("rulesToRun", ["rule1"])
    validate_rules(rules)

    use_mp = bool(sc.get("useMultiprocessing", True))
    max_cores = int(sc.get("max_cores", 8))
    if max_cores < 1:
        raise ValueError("scoring.max_cores must be >= 1")

    # Explicit core pinning config (optional). If omitted, first max_cores available CPUs are used.
    core_ids = sc.get("core_ids", None)
    if core_ids is None:
        core_ids = get_default_core_ids(max_cores)
    else:
        if len(core_ids) != max_cores:
            raise ValueError("len(scoring.core_ids) must equal scoring.max_cores")
        core_ids = [int(x) for x in core_ids]

    # Pin main process too
    set_affinity_or_warn(core_ids)
    print(f"CPU affinity pinned to cores: {core_ids}")
    print(f"max_cores={max_cores}, useMultiprocessing={use_mp}")

    input_csv = os.path.join(os.getcwd(), cfg["input"]["reactionCSV"])
    if not os.path.exists(input_csv):
        raise FileNotFoundError(input_csv)

    df = load_reactions(
        input_csv,
        test_mode=bool(sc.get("testMode", False)),
        test_rows=int(sc.get("testRows", 100)),
    )

    out = run_scoring(
        df=df,
        rules=rules,
        threshold=float(sc.get("labelThreshold", 0.5)),
        use_mp=use_mp,
        max_cores=max_cores,
        mp_chunk_size=int(sc.get("mpChunkSize", 64)),
        cache_size=int(sc.get("cacheSize", 200000)),
        core_ids=core_ids,
    )

    save_outputs(out, cfg["input"], cfg["output"], cfg["filtering"], rules)
    print("Job completed.")


if __name__ == "__main__":
    main()