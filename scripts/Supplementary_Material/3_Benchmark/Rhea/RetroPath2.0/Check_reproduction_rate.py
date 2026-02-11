import pandas as pd
from rdkit import Chem
import json


def load_reaction_smiles(folder_path):
    df = pd.read_csv(folder_path, usecols=["Reaction SMILES"])
    reaction_smiles = df["Reaction SMILES"].dropna().astype(str).str.strip().tolist()
    print("Number of reaction SMILES loaded", len(reaction_smiles))
    return reaction_smiles


def remove_stereo_and_explicit_hydrogens(smiles: str, canonical: bool = True) -> str:
    mol = Chem.MolFromSmiles(smiles)
    # strip stereochemistry
    Chem.RemoveStereochemistry(mol)
    # remove any explicit hydrogens
    mol = Chem.RemoveHs(mol)
    # return plain SMILES
    return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=canonical)


def standard_inchi_transform(smiles):  # return RetroPath normalized tautomers
    mol = Chem.MolFromSmiles(smiles)
    inchi = Chem.MolToInchi(mol)
    mol = Chem.MolFromInchi(inchi, sanitize=True)
    if mol is not None:
        return Chem.MolToSmiles(mol)
    else:
        return "None"


if __name__ == "__main__":
    rxn_smiles = load_reaction_smiles("combined_results.csv")
    unique_smiles = list()
    for i in rxn_smiles:
        if i not in unique_smiles:
            unique_smiles.append(i)
    print("Number of unique SMILES", len(unique_smiles))

    retropath_predictions = dict()
    for rxn in unique_smiles:
        clean_rea = remove_stereo_and_explicit_hydrogens(rxn.split(">>")[0])
        if clean_rea not in retropath_predictions:
            retropath_predictions[clean_rea] = list()

        pros = rxn.split(">>")[1].split(".")
        clean_pros = [remove_stereo_and_explicit_hydrogens(pro) for pro in pros]
        retropath_predictions[clean_rea].append(clean_pros)

    # check against rhea subset

    with open("Rhea_random_1000.json", "r") as f:
        Rhea_reactions = json.load(f)

    success_count = 0
    skipped_count = 0
    for rxn in Rhea_reactions:
        rhea_rea_st = standard_inchi_transform(rxn.split(">>")[0])
        rhea_pros = rxn.split(">>")[1].split(".")
        rhea_pros_st = [standard_inchi_transform(p) for p in rhea_pros]
        if rhea_rea_st not in retropath_predictions:
            skipped_count += 1
            continue

        pridicted_sets = retropath_predictions[rhea_rea_st]
        for pridicted_set in pridicted_sets:
            if set(rhea_pros_st).issubset(set(pridicted_set)):
                success_count += 1
                break

    print("Success count", success_count)
    print("Success rate", success_count / len(Rhea_reactions))
    print("Reactant without result", skipped_count)
