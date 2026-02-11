import sys
from pathlib import Path
from rdkit import Chem
import csv

# Add the local doranet path
DORANET_PATH = Path("/users/sghosh6/DTRA_project/MACAW/doranet")
sys.path.insert(0, str(DORANET_PATH))

import doranet.modules.enzymatic as enzymatic
import doranet.modules.post_processing as post_processing

# Phenylalanine → Tyrosine (single step: aromatic hydroxylation)
user_starters = {'NC(Cc1ccccc1)C(=O)O'}  # L-Phenylalanine
user_target = {'NC(Cc1ccc(O)cc1)C(=O)O'}  # L-Tyrosine
user_helpers = {'O', 'O=O', '[H][H]'}

job_name = "phe_to_tyr_test"

print(f"Starter: {list(user_starters)[0]}")
print(f"Target:  {list(user_target)[0]}")

forward_network = enzymatic.generate_network(
    job_name=job_name,
    starters=user_starters,
    gen=2,
    max_atoms={'C': 12, 'N': 3, 'O': 5},
    direction="forward",
    targets=user_target,
    ruleset="JN3604IMT"
)

smiles_list = [mol.uid for mol in forward_network.mols]
target_smiles = list(user_target)[0]
target_found = target_smiles in smiles_list

print(f"\nGenerated {len(smiles_list)} molecules")
print(f"Target found: {target_found}")

for i, smi in enumerate(smiles_list[:15]):
    marker = " <-- TARGET" if smi == target_smiles else ""
    marker = " <-- STARTER" if smi in user_starters else marker
    print(f"  {i+1}. {smi}{marker}")

output_path = Path(f"{job_name}_molecules.csv")
with open(output_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['SMILES', 'Is_Target', 'Is_Starter'])
    for smi in smiles_list:
        writer.writerow([smi, smi == target_smiles, smi in user_starters])

print(f"\nSaved to {output_path}")

post_processing.one_step(
    networks={forward_network},
    total_generations=2,
    starters=user_starters,
    helpers=user_helpers,
    target=user_target,
    job_name=job_name,
)