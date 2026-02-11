from syntheseus import Molecule
from syntheseus.reaction_prediction.inference import LocalRetroModel, ChemformerModel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
import numpy as np
import time


def clean_molecule(smiles):
    if smiles == "OO":
        return "O=O"
    if Chem.MolFromSmiles(smiles).HasSubstructMatch(
        Chem.MolFromSmarts("[C](=[O])[O][OH]")
    ):
        return "O=O"

    neutral = Uncharger().uncharge(Chem.MolFromSmiles(smiles))
    return Chem.MolToSmiles(neutral)


def get_predicted_precursors_smiles(product_smiles: str, model, num_results):
    mol = Molecule(product_smiles)
    [predictions] = model(
        [mol], num_results=num_results
    )  # unpack the single-batch result

    return [
        sorted([clean_molecule(r.smiles) for r in pred.reactants])
        for pred in predictions
    ]


data = np.genfromtxt(
    "Sanitized_reactions_onlyCHO.csv",
    comments="?",
    dtype=str,
    delimiter=",",
    skip_header=0,
)

if __name__ == "__main__":
    start_time = time.time()

    model = LocalRetroModel()

    print("Loaded reaction number:", len(data))

    correct_num = 0
    for rxn in data:
        reactants, products = rxn.split(">>")
        reactants = reactants.split(".")
        reactants.sort()

        target = products
        precursor_lists = get_predicted_precursors_smiles(target, model, num_results=5)

        for prec in precursor_lists:
            if prec == reactants:
                correct_num += 1
                break

    correction_rate = correct_num / len(data)
    end_time = time.time()
    elapsed_time = end_time - start_time

    print("Time used: " + "{:.2f}".format(elapsed_time) + " seconds")
    print("Correction number:", correct_num)
    print("Correction rate:", correction_rate)
