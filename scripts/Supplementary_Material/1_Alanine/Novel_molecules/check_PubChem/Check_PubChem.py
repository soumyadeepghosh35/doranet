from rdkit import Chem
import json
import pandas as pd


with open("3gen_SMILES_C1-6.json", "r") as file:
    data = json.load(file)

mol_set = set()

for i in data:
    mol_set.add(Chem.MolToSmiles(Chem.MolFromSmiles(i)))

print("Total number of DORAnet molecules", len(mol_set))


dataframe1 = pd.read_excel("PubChem_compound_formulaquery_H-C1-6O0-5N0-4Cl0-4.xlsx")

pub_smiles = dataframe1["canonicalsmiles"]

pub_set = set()


def clean_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        Chem.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=False)

    else:
        return "O"


for i in pub_smiles:
    pub_set.add(clean_smiles(i))
    if "." not in i:
        if "[O-]" in i:
            i = i.replace("[O-]", "O")
        pub_set.add(clean_smiles(i))
    else:
        mols = i.split(".")
        for j in mols:
            if "HH" not in j:
                if "[O-]" in j:
                    j = j.replace("[O-]", "O")
                pub_set.add(clean_smiles(j))


know_mol = 0
for i in mol_set:
    if i in pub_set:
        know_mol += 1


print("Novel percentage: ", 1 - know_mol / len(mol_set))


result_dict = dict()
total_dict = dict()
unknown_mols = {
    0: list(),
    1: list(),
    2: list(),
    3: list(),
    4: list(),
    5: list(),
    6: list(),
}

for i in mol_set:
    carbon_num = len(
        Chem.MolFromSmiles(i).GetSubstructMatches(Chem.MolFromSmarts("[#6]"))
    )
    if carbon_num in total_dict:
        total_dict[carbon_num] = total_dict[carbon_num] + 1
    else:
        total_dict[carbon_num] = 1

    if i in pub_set:
        if carbon_num in result_dict:
            result_dict[carbon_num] = result_dict[carbon_num] + 1
        else:
            result_dict[carbon_num] = 1
    else:
        unknown_mols[carbon_num].append(i)


print("Novel percentage:")

for key in result_dict:
    print(key, len(unknown_mols[key]) / total_dict[key])
