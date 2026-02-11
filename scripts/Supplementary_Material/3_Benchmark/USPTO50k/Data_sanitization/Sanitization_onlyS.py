import numpy as np
from rdkit import Chem
import time
import re

print("Sanitizing dataset for S reactions")
start_time = time.time()

filtered_list = []

data = np.genfromtxt(
    "dataset_subset.csv",
    comments="?",
    dtype=str,
    delimiter=",",
    skip_header=1,
)

Reaction_smiles_list = data[:, 0]
ReactantSet = data[:, 1]
ReactantSet_idx = []


for i in ReactantSet:
    all_numbers = re.findall(r"\[(.*?)\]", i)
    Reactant_idx = all_numbers[0].split("$")
    for idx, j in enumerate(Reactant_idx):
        Reactant_idx[idx] = int(j)
    ReactantSet_idx.append(Reactant_idx)

correct_num = 0
if len(ReactantSet) != len(Reaction_smiles_list):
    print("Number of reactions error")
checked_reactions = 0
nitrogen_num = 0

for sm_idx, smarts in enumerate(Reaction_smiles_list):
    if len(smarts.split(">")) == 3:
        mol_error = 0
        rxn_smiles = str()

        reactants, agents, products = smarts.split(">")

        reactants = reactants.split(".")
        temp_reactants = []
        for rea_idx in ReactantSet_idx[sm_idx]:
            temp_reactants.append(reactants[rea_idx])
        reactants = temp_reactants

        for idx, reactant in enumerate(reactants):
            if Chem.MolFromSmiles(reactant).HasSubstructMatch(
                Chem.MolFromSmarts("[C](=[O])[O][OH]")
            ):
                reactant = "O=O"
            if reactant == "OO":
                reactant = "O=O"

            mol = Chem.MolFromSmiles(reactant)
            if mol is not None:
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)
                reactants[idx] = Chem.MolToSmiles(mol, isomericSmiles=False)
            else:
                mol_error += 1
                reactants.remove(reactant)

        reactants = list(set(reactants))
        ########################

        n_flag = False
        remove_ioin = False  # <--
        O_minus_add_H = True  # <--
        if O_minus_add_H is True:
            for idx, rea in enumerate(reactants):
                if "[O-]" in rea:
                    reactants[idx] = rea.replace("[O-]", "O")

        for rea in reactants:
            # if "N" in rea or "n" in rea:
            if (
                "N" in rea
                or "n" in rea
                or "S" in rea
                or "s" in rea
                or "B" in rea
                or "P" in rea
            ):  # --N and S
                n_flag = True

        if remove_ioin is True:
            temp = []
            for rea in reactants:
                if "-" not in rea and "+" not in rea:
                    temp.append(rea)
            reactants = temp
        for rea in reactants:
            rxn_smiles = rxn_smiles + "." + rea
        rxn_smiles = rxn_smiles[1:]
        rxn_smiles = rxn_smiles + ">>"

        Keep_N = True  # <--
        if n_flag is False or Keep_N is True:
            products = products.split(".")
            for idx, product in enumerate(products):
                mol = Chem.MolFromSmiles(product)
                if mol is not None:
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
                    products[idx] = Chem.MolToSmiles(mol, isomericSmiles=False)
                else:
                    mol_error += 1
                    products.remove(product)

            if O_minus_add_H is True:
                for idx, pro in enumerate(products):
                    if "[O-]" in pro:
                        products[idx] = pro.replace("[O-]", "O")

            if remove_ioin is True:
                temp = []
                for pro in products:
                    if "-" not in pro and "+" not in pro:
                        temp.append(pro)
                products = temp

            for pro in products:
                rxn_smiles = rxn_smiles + "." + pro

            check_balance = True  # <--
            left_carbon = 0
            right_carbon = 0
            if check_balance is True:
                for i in reactants:
                    left_carbon += len(
                        Chem.MolFromSmiles(i).GetSubstructMatches(
                            Chem.MolFromSmarts("[#6]")
                        )
                    )
                for i in products:
                    right_carbon += len(
                        Chem.MolFromSmiles(i).GetSubstructMatches(
                            Chem.MolFromSmarts("[#6]")
                        )
                    )

            if (
                reactants
                and products
                and left_carbon == right_carbon
                and reactants != products
            ):
                rea_N_H_num = 0
                pro_N_H_num = 0
                rea_N_num = 0
                pro_N_num = 0
                rea_N_d_num = 0
                pro_N_d_num = 0
                rea_N_elem = []
                pro_N_elem = []

                rea_P_H_num = 0
                pro_P_H_num = 0
                rea_P_num = 0
                pro_P_num = 0
                rea_P_d_num = 0
                pro_P_d_num = 0
                rea_P_elem = []
                pro_P_elem = []

                rea_S_H_num = 0
                pro_S_H_num = 0
                rea_S_num = 0
                pro_S_num = 0
                rea_S_d_num = 0
                pro_S_d_num = 0
                rea_S_elem = []
                pro_S_elem = []

                rea_B_H_num = 0
                pro_B_H_num = 0
                rea_B_num = 0
                pro_B_num = 0
                rea_B_d_num = 0
                pro_B_d_num = 0
                rea_B_elem = []
                pro_B_elem = []

                for rea_smiles in reactants:
                    mol = Chem.MolFromSmiles(rea_smiles)
                    for atom in mol.GetAtoms():
                        if Chem.Atom.GetAtomicNum(atom) == 7:  # N
                            rea_N_H_num += atom.GetTotalNumHs()
                            rea_N_num += 1
                            rea_N_d_num += Chem.Atom.GetTotalDegree(atom)
                            for nei in Chem.Atom.GetNeighbors(atom):
                                rea_N_elem.append(Chem.Atom.GetAtomicNum(nei))

                        if Chem.Atom.GetAtomicNum(atom) == 15:  # P
                            rea_P_H_num += atom.GetTotalNumHs()
                            rea_P_num += 1
                            rea_P_d_num += Chem.Atom.GetTotalDegree(atom)
                            for nei in Chem.Atom.GetNeighbors(atom):
                                rea_P_elem.append(Chem.Atom.GetAtomicNum(nei))

                        if Chem.Atom.GetAtomicNum(atom) == 16:  # S
                            rea_S_H_num += atom.GetTotalNumHs()
                            rea_S_num += 1
                            rea_S_d_num += Chem.Atom.GetTotalDegree(atom)
                            for nei in Chem.Atom.GetNeighbors(atom):
                                rea_S_elem.append(Chem.Atom.GetAtomicNum(nei))

                        if Chem.Atom.GetAtomicNum(atom) == 5:  # B
                            rea_B_H_num += atom.GetTotalNumHs()
                            rea_B_num += 1
                            rea_B_d_num += Chem.Atom.GetTotalDegree(atom)
                            for nei in Chem.Atom.GetNeighbors(atom):
                                rea_B_elem.append(Chem.Atom.GetAtomicNum(nei))

                for pro_smiles in products:
                    mol = Chem.MolFromSmiles(pro_smiles)
                    for atom in mol.GetAtoms():
                        if Chem.Atom.GetAtomicNum(atom) == 7:
                            pro_N_H_num += atom.GetTotalNumHs()
                            pro_N_num += 1
                            pro_N_d_num += Chem.Atom.GetTotalDegree(atom)
                            for nei in Chem.Atom.GetNeighbors(atom):
                                pro_N_elem.append(Chem.Atom.GetAtomicNum(nei))

                        if Chem.Atom.GetAtomicNum(atom) == 15:
                            pro_P_H_num += atom.GetTotalNumHs()
                            pro_P_num += 1
                            pro_P_d_num += Chem.Atom.GetTotalDegree(atom)
                            for nei in Chem.Atom.GetNeighbors(atom):
                                pro_P_elem.append(Chem.Atom.GetAtomicNum(nei))

                        if Chem.Atom.GetAtomicNum(atom) == 16:
                            pro_S_H_num += atom.GetTotalNumHs()
                            pro_S_num += 1
                            pro_S_d_num += Chem.Atom.GetTotalDegree(atom)
                            for nei in Chem.Atom.GetNeighbors(atom):
                                pro_S_elem.append(Chem.Atom.GetAtomicNum(nei))

                        if Chem.Atom.GetAtomicNum(atom) == 5:
                            pro_B_H_num += atom.GetTotalNumHs()
                            pro_B_num += 1
                            pro_B_d_num += Chem.Atom.GetTotalDegree(atom)
                            for nei in Chem.Atom.GetNeighbors(atom):
                                pro_B_elem.append(Chem.Atom.GetAtomicNum(nei))

                rea_N_elem.sort()
                pro_N_elem.sort()
                rea_P_elem.sort()
                pro_P_elem.sort()
                rea_S_elem.sort()
                pro_S_elem.sort()
                rea_B_elem.sort()
                pro_B_elem.sort()

                if all(
                    [
                        rea_N_H_num == pro_N_H_num,
                        rea_P_H_num == pro_P_H_num,
                        #      rea_S_H_num == pro_S_H_num,
                        rea_B_H_num == pro_B_H_num,
                        rea_N_num == pro_N_num,
                        rea_P_num == pro_P_num,
                        rea_S_num == pro_S_num,
                        rea_B_num == pro_B_num,
                        rea_N_d_num == pro_N_d_num,
                        rea_P_d_num == pro_P_d_num,
                        #       rea_S_d_num == pro_S_d_num,
                        rea_B_d_num == pro_B_d_num,
                        rea_N_elem == pro_N_elem,
                        rea_P_elem == pro_P_elem,
                        #       rea_S_elem == pro_S_elem,
                        rea_B_elem == pro_B_elem,
                    ]
                ) and any(
                    [
                        rea_S_H_num != pro_S_H_num,
                        rea_S_d_num != pro_S_d_num,
                        rea_S_elem != pro_S_elem,
                    ]
                ):
                    checked_reactions += 1

                    if n_flag is True:
                        nitrogen_num += 1

                    correct_num += 1
                    filtered_list.append(rxn_smiles)


filtered_list2 = list()

for sm_idx, smarts in enumerate(filtered_list):
    if len(smarts.split(">>.")) == 2:
        mol_error = 0
        rxn_smiles = str()
        reactants, products = smarts.split(">>.")
        reactants = reactants.split(".")

        for idx, reactant in enumerate(reactants):
            if Chem.MolFromSmiles(reactant).HasSubstructMatch(
                Chem.MolFromSmarts("[C](=[O])[O][OH]")
            ):
                reactant = "O=O"
            if reactant == "OO":
                reactant = "O=O"

            mol = Chem.MolFromSmiles(reactant)
            if mol is not None:
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)
                reactants[idx] = Chem.MolToSmiles(mol, isomericSmiles=False)
            else:
                mol_error += 1
                reactants.remove(reactant)

        reactants = list(set(reactants))
        products = products.split(".")
        products = list(set(products))

        reactants, products = (
            [i for i in reactants if i not in products],
            [j for j in products if j not in reactants],
        )

        n_flag = False
        Mg_flag = False
        remove_ioin = False  # <--
        O_minus_add_H = True  # <--
        if O_minus_add_H is True:
            for idx, rea in enumerate(reactants):
                if "[S-]" in rea:
                    reactants[idx] = rea.replace("[S-]", "S")

        for rea in reactants:
            # if "N" in rea or "n" in rea:
            if (
                "N" in rea
                or "n" in rea
                or "S" in rea
                or "s" in rea
                or "B" in rea
                or "P" in rea
            ):  # --
                n_flag = True
            if "Mg" in rea or "Zn" in rea:  # --
                Mg_flag = True
            if "[C-]#N" in rea or "Cu" in rea or "Sn" in rea:  # --
                Mg_flag = True

        if remove_ioin is True:
            temp = []
            for rea in reactants:
                if "-" not in rea and "+" not in rea:
                    temp.append(rea)
            reactants = temp
        for rea in reactants:
            rxn_smiles = rxn_smiles + "." + rea
        rxn_smiles = rxn_smiles[1:]
        rxn_smiles = rxn_smiles + ">>"

        Keep_Mg = False  # <--
        if Mg_flag is False or Keep_Mg is True:
            for idx, product in enumerate(products):
                mol = Chem.MolFromSmiles(product)
                if mol is not None:
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
                    products[idx] = Chem.MolToSmiles(mol, isomericSmiles=False)
                else:
                    mol_error += 1
                    products.remove(product)

            if O_minus_add_H is True:
                for idx, pro in enumerate(products):
                    if "[S-]" in rea:
                        reactants[idx] = rea.replace("[S-]", "S")

            if remove_ioin is True:
                temp = []
                for pro in products:
                    if "-" not in pro and "+" not in pro:
                        temp.append(pro)
                products = temp

            for pro in products:
                rxn_smiles = rxn_smiles + "." + pro

            # remove a smiles appeared on both sides
            for j in reactants:
                if j in products:
                    reactants.remove(j)
                    products.remove(j)

            check_balance = True  # <--
            left_carbon = 0
            right_carbon = 0
            if check_balance is True:
                for i in reactants:
                    left_carbon += len(
                        Chem.MolFromSmiles(i).GetSubstructMatches(
                            Chem.MolFromSmarts("[#6]")
                        )
                    )
                for i in products:
                    right_carbon += len(
                        Chem.MolFromSmiles(i).GetSubstructMatches(
                            Chem.MolFromSmarts("[#6]")
                        )
                    )

            if (
                reactants
                and products
                and left_carbon == right_carbon
                and reactants != products
            ):
                filtered_list2.append(".".join(reactants) + ">>" + ".".join(products))

end_time = time.time()
elapsed_time = (end_time - start_time) / 60

print("S reactions subset ready, number of reactinos:", len(filtered_list2))
print()

np.savetxt("Sanitized_reactions_onlyS.csv", filtered_list2, delimiter=", ", fmt="% s")
