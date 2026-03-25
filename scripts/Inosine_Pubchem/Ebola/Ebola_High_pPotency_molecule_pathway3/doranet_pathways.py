#!/usr/bin/env python
import os
import sys
import time
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import csv

# Restrict threading per process
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

DORANET_PATH = Path(r"/users/sghosh6/DTRA_project/MACAW/doranet")
sys.path.insert(0, str(DORANET_PATH))

import doranet.modules.enzymatic as enzymatic
import doranet.modules.post_processing as post_processing

start_time = time.time()

job_name = "Ebola_High_pPotency_molecule_pathway3"

user_starters = {'C1=NC2=C(C(=O)N1)N=CN2[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O'}

user_helpers = {'N', 'O=[N+]([O-])O', 'N#CO', '[C-]#[O+]', 'O=C=O', 'Br', 'C=C', 'CO', 'O', 'N#N', 'S', '[Br][Br]', '[H][H]', 'O=NO', 'C#N', 'C=O', 'O=S(O)O', 'O=O', 'O=S(=O)(O)O', 'O=S=O', 'NO'}

user_target = {"CC(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O)[C@@H](O)C(=O)NCCC(=O)NCCS"}

max_atoms = {'C': 15, 'N': 6, 'O': 8, 'S': 3}
generations = 3
ruleset = "JN3604IMT"

print(f"Starter: {list(user_starters)[0]}")
print(f"Target:  {list(user_target)[0]}")
print(f"Job:     {job_name}")

forward_network = enzymatic.generate_network(
    job_name=job_name,
    starters=user_starters,
    gen=generations,
    max_atoms=max_atoms,
    direction="forward",
    targets=user_target,
    ruleset=ruleset
)

smiles_list = list(user_starters) + [mol.uid for mol in forward_network.mols if mol.uid not in user_starters]
print(f"Generated {len(smiles_list) - len(user_starters)} new molecules + {len(user_starters)} starters")

all_targets = set(smiles_list) - user_starters - user_helpers
print(f"Using {len(all_targets)} generated molecules as targets for pathway finding")

output_path = Path(f"{job_name}_molecules.csv")
with open(output_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['SMILES', 'Is_Starter', 'MolFormula', 'MolWeight', 'NumHeavyAtoms'])
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mol_weight = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
            num_heavy = mol.GetNumHeavyAtoms()
        else:
            formula, mol_weight, num_heavy = "N/A", 0, 0
        writer.writerow([smi, smi in user_starters, formula, mol_weight, num_heavy])

print(f"Saved molecules to {output_path}")

if all_targets:
    post_processing.one_step(
        networks={forward_network},
        total_generations=generations,
        starters=user_starters,
        helpers=user_helpers,
        target=all_targets,
        job_name=job_name,
    )
    print(f"Pathway files generated with prefix: {job_name}")
else:
    print("No targets found for pathway generation")

print(f"Time: {time.time() - start_time:.2f} s")
