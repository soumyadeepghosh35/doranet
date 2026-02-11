import zipfile
import pandas as pd

# Paths
ZIP_PATH = "ci6b00564_si_002.zip"
INNER_CSV = "data/dataSetB.csv"
OUTPUT_CSV = "dataset_subset.csv"


def extract_and_match_manual(
    zip_path: str, inner_csv: str, cols: list[int], output_path: str
):
    with zipfile.ZipFile(zip_path) as z:
        with z.open(inner_csv) as f:
            df = pd.read_csv(f, dtype=str)

    subset = df.iloc[:, cols]

    subset.iloc[:, 1] = subset.iloc[:, 1].str.replace(", ", "$", regex=False)

    subset.to_csv(output_path, index=False)


if __name__ == "__main__":
    extract_and_match_manual(ZIP_PATH, INNER_CSV, [2, 3], OUTPUT_CSV)
