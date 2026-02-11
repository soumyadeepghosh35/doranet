from rdkit import Chem
import json
from multiprocessing import Pool


def clean_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=False)
    else:
        return "O"


def sanitize_smiles(smi):
    return clean_smiles(smi)


def main():
    n_cpus = 8

    # Load all files 1.smi through 12.smi
    smiles_list = []
    file_count = 0
    for i in range(1, 13):
        file_count += 1
        filename = f"{i}.smi"
        with open(filename, "r") as f:
            smiles_list.extend([line.strip() for line in f])
    print(file_count, " GDB files loaded")
    print("Total number of GDB molecules:", len(smiles_list))

    # Raw C1–6 filtering & distribution
    GDB_C1_6 = set()
    GDB_C1_6_dis = {i: 0 for i in range(1, 7)}
    for smi in smiles_list:
        C_count = smi.count("C") + smi.count("c")
        if 0 < C_count < 7:
            GDB_C1_6.add(smi)
            GDB_C1_6_dis[C_count] += 1

    with open("GDB13_C1-6_raw_smiles.json", "w", encoding="utf-8") as f:
        json.dump(list(GDB_C1_6), f, indent=4)
    print("GDB13 C1-6 carbon num distribution:")
    print(GDB_C1_6_dis)

    # Sanitization using multiprocessing
    with Pool(processes=n_cpus) as pool:
        cleaned_list = pool.map(sanitize_smiles, smiles_list)
    GDB_smiles = set(cleaned_list)
    print(f"Sanitized {len(GDB_smiles)} GDB SMILES")

    # Cleaned C1–6 filtering
    GDB_cleaned_C1_6_set = set()
    for smi in GDB_smiles:
        C_count = smi.count("C") + smi.count("c")
        if 0 < C_count < 7:
            GDB_cleaned_C1_6_set.add(smi)

    with open("GDB13_C1-6_clean_smiles.json", "w", encoding="utf-8") as f:
        json.dump(list(GDB_cleaned_C1_6_set), f, indent=4)
    return GDB_cleaned_C1_6_set


if __name__ == "__main__":
    # Load DORAnet result
    with open("3gen_SMILES_C1-6.json", "r") as file:
        data = json.load(file)
    dora_set = set()
    for i in data:
        dora_set.add(Chem.MolToSmiles(Chem.MolFromSmiles(i)))

    print("total number of dora molecules", len(dora_set))

    # Load GDB
    GDB_cleaned_C1_6_set = main()

    # Compare
    know_mol = 0
    for i in dora_set:
        if i in GDB_cleaned_C1_6_set:
            know_mol += 1

    print("Known molecules: ", know_mol)
    print("Novel percentage: ", 1 - know_mol / len(dora_set))

    total_dict = dict()  # carbon num distribution in DORAnet mols
    unknown_mols = {
        1: list(),
        2: list(),
        3: list(),
        4: list(),
        5: list(),
        6: list(),
    }

    for i in dora_set:
        carbon_num = len(
            Chem.MolFromSmiles(i).GetSubstructMatches(Chem.MolFromSmarts("[#6]"))
        )
        if carbon_num in total_dict:
            total_dict[carbon_num] = total_dict[carbon_num] + 1
        else:
            total_dict[carbon_num] = 1

        if i not in GDB_cleaned_C1_6_set:
            unknown_mols[carbon_num].append(i)

    print("DORAnet total:", total_dict)

    print("Novel percentage:")
    for key in total_dict:
        print(key, len(unknown_mols[key]) / total_dict[key])
