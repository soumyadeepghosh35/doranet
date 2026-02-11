import csv
import subprocess
import os
import shutil
import sys
import time

start_time = time.time()

# === CONFIGURE HERE ===
CLI_PATH = "./retropath_1.0.0-alpha_linux-x64/RetroPath.Cli"  # path to the CLI binary
RULES_FILE = "retrorules_rr01_rp2_flat_all.csv"
SINK_FILE = "sink.csv"  # empty sink file
MASTER_SOURCE = "source_rhea1000.csv"
PATHWAY_LENGTH = 1
WORKDIR = "per_source_runs"  # temporary per-source working dirs
COMBINED_RESULTS = "combined_results.csv"
SKIP_LOG = "skipped_sources.log"  # molecules skipped due to errors
OTHER_ERRORS_LOG = "other_failures.log"  # non-specific failures
TIMEOUT_SECONDS = 600  # timeout for each subprocess.run
# ========================

# ensure workdir exists
os.makedirs(WORKDIR, exist_ok=True)

# read master source list
sources = []
with open(MASTER_SOURCE, newline="") as f:
    reader = csv.reader(f)
    for row in reader:
        if not row or all(not cell.strip() for cell in row):
            continue
        # skip header heuristically
        first = row[0].strip().lower()
        if (
            "name" in first
            and len(row) > 1
            and row[1].strip().lower().startswith("inchi")
        ) or first.startswith("inchi"):
            continue
        if len(row) == 1:
            name = f"src_{len(sources)}"
            inchi = row[0].strip()
        else:
            name = row[0].strip() or f"src_{len(sources)}"
            inchi = row[1].strip()
        sources.append((name, inchi))

print("number of sources", len(sources))
if not sources:
    print(f"[ERROR] No source entries parsed from {MASTER_SOURCE}", file=sys.stderr)
    sys.exit(1)


# helpers
def safe_name(s):
    # sanitize name for filesystem usage
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)[:80]


header_written = False

# open skip/other error logs
skip_f = open(SKIP_LOG, "w")
other_f = open(OTHER_ERRORS_LOG, "w")

for idx, (name, inchi) in enumerate(sources, 1):
    short = safe_name(name)
    run_tag = f"{idx:04d}_{short}"
    src_file = os.path.join(WORKDIR, f"source_{run_tag}.csv")
    outdir = os.path.join(WORKDIR, f"out_{run_tag}")
    # write single-molecule source file (name,InChI)
    with open(src_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Name", "InChI"])
        writer.writerow([name, inchi])

    print(f"[{idx}/{len(sources)}] Running for source '{name}'")

    # build command
    cmd = [
        CLI_PATH,
        RULES_FILE,
        src_file,
        SINK_FILE,
        str(PATHWAY_LENGTH),
        "--output-dir",
        outdir,
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=TIMEOUT_SECONDS
        )
    except subprocess.TimeoutExpired:
        msg = f"{name},{inchi},TIMEOUT_{TIMEOUT_SECONDS}s,{run_tag}\n"
        skip_f.write(msg)
        skip_f.flush()
        print(f"  -> Skipped (timeout after {TIMEOUT_SECONDS}s).")
        # cleanup
        if os.path.isdir(outdir):
            shutil.rmtree(outdir, ignore_errors=True)
        os.remove(src_file)
        continue

    combined_output = proc.stdout + "\n" + proc.stderr

    # check for the specific error to skip
    if "Sequence contains more than one element" in combined_output:
        msg = f"{name},{inchi},SKIPPED_SEQUENCE_MULTIPLE,{run_tag}\n"
        skip_f.write(msg)
        skip_f.flush()
        print(f"  -> Skipped (sequence contains more than one element).")
        # clean up this run's folder if exists
        if os.path.isdir(outdir):
            shutil.rmtree(outdir, ignore_errors=True)
        # remove source file
        try:
            os.remove(src_file)
        except OSError:
            pass
        continue

    # other non-zero return code
    if proc.returncode != 0:
        msg = f"{name},{inchi},ERROR_RC{proc.returncode},{run_tag},stderr={proc.stderr.strip().replace(os.linesep,' ')}\n"
        other_f.write(msg)
        other_f.flush()
        print(f"  -> Failed with return code {proc.returncode}, logged.")
        # keep the folder for inspection, but continue
        continue

    # expect results.csv under outdir
    result_path = os.path.join(outdir, "results.csv")
    if not os.path.exists(result_path):
        # fallback: maybe inside a 'results' subfolder
        alt = os.path.join(outdir, "results", "results.csv")
        if os.path.exists(alt):
            result_path = alt
        else:
            msg = f"{name},{inchi},NO_RESULTS,{run_tag}\n"
            other_f.write(msg)
            other_f.flush()
            print(f"  -> No results.csv found, logged.")
            continue

    # merge into combined_results.csv
    try:
        with open(result_path, newline="") as rf:
            reader = csv.reader(rf)
            rows = list(reader)
            if not rows:
                print(f"  -> results.csv empty for {name}, skipping.")
                continue
            # header
            header = rows[0]
            data = rows[1:]
            if not header_written:
                with open(COMBINED_RESULTS, "w", newline="") as wf:
                    writer = csv.writer(wf)
                    writer.writerow(header)
                    writer.writerows(data)
                header_written = True
            else:
                with open(COMBINED_RESULTS, "a", newline="") as wf:
                    writer = csv.writer(wf)
                    writer.writerows(data)
        print(f"  -> Appended {len(data)} rows from result.")
    except Exception as e:
        other_f.write(f"{name},{inchi},EXCEPTION_MERGE,{run_tag},{e}\n")
        other_f.flush()
        print(f"  -> Exception merging result: {e}", file=sys.stderr)

    # cleanup per-molecule artifacts
    try:
        shutil.rmtree(outdir, ignore_errors=True)
    except Exception:
        pass
    try:
        os.remove(src_file)
    except Exception:
        pass

# close logs
skip_f.close()
other_f.close()

print("\n=== Summary ===")
print(f"Total sources processed: {len(sources)}")
if os.path.exists(COMBINED_RESULTS):
    print(f"Combined results written to: {COMBINED_RESULTS}")
else:
    print("No successful results were merged.")
print(f"Skipped (Sequence contains more than one element): see {SKIP_LOG}")
print(f"Other failures: see {OTHER_ERRORS_LOG}")

end_time = time.time()
elapsed_time = (end_time - start_time) / 60
print("time used:", "{:.2f}".format(elapsed_time), " minutes")
