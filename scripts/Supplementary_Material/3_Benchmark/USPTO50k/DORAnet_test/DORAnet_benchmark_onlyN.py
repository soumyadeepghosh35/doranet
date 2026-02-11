import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import time
from datetime import datetime
from Reaction_Smarts import op_smarts
import doranet as dn


def predictions(
    smarts_list,
    starters_list,
):
    engine = dn.create_engine()
    network = engine.new_network()

    for smiles in starters_list:
        network.add_mol(engine.mol.rdkit(smiles))

    for smarts in smarts_list:
        network.add_op(
            engine.op.rdkit(
                smarts.smarts,
                drop_errors=True,
            ),
        )

    strat = engine.strat.cartesian(network)

    strat.expand(
        num_iter=1,
    )

    # get reactions from network
    total_sets_of_products_smiles = []
    for rxn in network.rxns:
        products = rxn.products
        name_value = "name"
        products_smiles = [network.mols[i].uid for i in products]

        total_sets_of_products_smiles.append((name_value, products_smiles))

    return_list = []
    for name, sets in total_sets_of_products_smiles:
        sets = list(set(sets))
        sets.sort()
        if (name, sets) not in return_list:
            return_list.append((name, sets))
    return return_list


print("Testing N reactions")

failed_list = []

start_time = time.time()

helper_smiles = [
    "[H][H]",
    "O",
    "Cl",
    "Br",
    "O=O",
    "I",
    "F",
    "[Cl][Cl]",
]


for idx, helper in enumerate(helper_smiles):
    mol = Chem.MolFromSmiles(helper)
    helper_smiles[idx] = Chem.MolToSmiles(mol)


data = np.genfromtxt(
    "Sanitized_reactions_onlyN.csv",
    comments="?",
    dtype=str,
    delimiter=",",
    skip_header=0,
)

Reaction_smiles_list = data
correct_num = 0
checked_reactions = 0
nitrogen_num = 0

unchecked_rxns = list()
rxn_name_dict = dict()
total_calls = 0

for sm_idx, smarts in enumerate(Reaction_smiles_list):
    mol_error = 0
    rxn_smiles = str()

    reactants, products = smarts.split(">>")
    reactants = reactants.split(".")
    products = products.split(".")

    checked_reactions += 1

    products.sort()
    correct_flag = False
    correct_rxn_name = str()
    generated_product_sets_list = predictions(op_smarts, reactants)  # --
    total_calls += 1

    for each_name, each_generated_set in generated_product_sets_list:
        if each_generated_set == products:
            correct_flag = True
            correct_rxn_name = each_name
            break

    Use_helper = True
    if correct_flag is False and Use_helper is True:
        for helper in helper_smiles:
            if correct_flag is False:
                if helper not in products:
                    prodcut_with_helper = products + [helper]
                    prodcut_with_helper.sort()
                    for each_name, each_generated_set in generated_product_sets_list:
                        if each_generated_set == prodcut_with_helper:
                            correct_flag = True
                            correct_rxn_name = each_name
                            break

    if correct_flag is False and Use_helper is True:
        for helper in helper_smiles:
            if helper not in reactants:
                generated_product_sets_list = predictions(
                    op_smarts, reactants + [helper]
                )
                total_calls += 1
                for each_name, each_generated_set in generated_product_sets_list:
                    if each_generated_set == products:
                        correct_flag = True
                        correct_rxn_name = each_name
                        break

                if correct_flag is False and (helper == "O=O" or helper == "[H][H]"):
                    for each_name, each_generated_set in generated_product_sets_list:
                        prodcut_with_helper = products + ["O"]

                        prodcut_with_helper.sort()

                        if each_generated_set == prodcut_with_helper:
                            correct_flag = True
                            correct_rxn_name = each_name
                            break

                if correct_flag is False and helper == "[H][H]":
                    for each_name, each_generated_set in generated_product_sets_list:
                        for h in ["F", "Cl", "Br", "I"]:
                            prodcut_with_helper = products + [h]
                            prodcut_with_helper.sort()

                            if each_generated_set == prodcut_with_helper:
                                correct_flag = True
                                correct_rxn_name = each_name
                                break

                if correct_flag is False and (
                    helper == "F" or helper == "Cl" or helper == "Br" or helper == "I"
                ):
                    for each_name, each_generated_set in generated_product_sets_list:
                        prodcut_with_helper = products + ["O"]

                        prodcut_with_helper.sort()

                        if each_generated_set == prodcut_with_helper:
                            correct_flag = True
                            correct_rxn_name = each_name
                            break

    if correct_flag is True:
        correct_num += 1
        rxn_name_dict[correct_rxn_name] = rxn_name_dict.get(correct_rxn_name, 0) + 1

    if correct_flag is False:
        failed_list.append(rxn_smiles)


correction_rate = correct_num / checked_reactions
end_time = time.time()
elapsed_time = (end_time - start_time) / 60


print("N reactions reproduction rate:", correction_rate)
print("time used:", "{:.2f}".format(elapsed_time), " minutes")
print("Total doranet calls:", total_calls)
print()

f = open("Result_onlyN_rxns.txt", "a")

f.write("\n")
f.write("Nitrogen reaction test result: " + str(datetime.now()))
f.write("  Time used: " + "{:.2f}".format(elapsed_time) + " minutes")
f.write("  Checked reactions: " + str(checked_reactions))
f.write("  Matched reactions: " + str(correct_num))
f.write("  Correction rate: " + str(correction_rate))
f.write("\n")
f.write("END")
f.close()
