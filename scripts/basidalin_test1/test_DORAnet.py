import sys
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.QED import qed
import csv

# Add the local doranet_ycMOD to the path (before pip-installed doranet)
DORANET_PATH = Path("/users/sghosh6/DTRA_project/MACAW/doranet")
sys.path.insert(0, str(DORANET_PATH))

import doranet.modules.enzymatic as enzymatic
import doranet.modules.synthetic as synthetic
import doranet.modules.post_processing as post_processing

# copy in SMILES string of the molecule you want to modify
user_starters = {'C1=C(C(=CC=O)OC1=O)N'}

# I typically use these for synthetic chemistry reactions
# for enzymatic reactions, the cofactor/ helper molecules needed are already in DORAnet by default
user_helpers = {'O','O=O','[H][H]','O=C=O','C=O','[C-]#[O+]','Br','[Br][Br]','CO',
                'C=C','O=S(O)O','N','O=S(=O)(O)O','O=NO','N#N','O=[N+]([O-])O','NO',
                'C#N','S','O=S=O','N#CO'}

# copy in SMILES string of the molecule you want to reach
# if you are only modifying the starting molecule, you don't need to specify a target
user_target = {'OC1=CC=CC=C1'}       

job_name = "basidalin_syhthesis"

forward_network = enzymatic.generate_network(
    job_name = job_name,
    starters = user_starters,
    gen = 2,
    direction = "forward",
    #targets = user_target,
    ruleset = "JN3604IMT")

'''
forward_network = synthetic.generate_network(
     job_name = job_name,
     starters = user_starters,
     helpers = user_helpers,
     gen = 1,
     direction = "forward",
)
'''

for mol in forward_network.mols:
    print(mol.uid)

# Extract SMILES from network object
smiles_list = [mol.uid for mol in forward_network.mols]

# Save to CSV
output_path = Path("doranet_generated_molecules.csv")
with open(output_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['SMILES', 'Source'])
    for smi in smiles_list:
        writer.writerow([smi, 'enzymatic_network'])

print(f"Saved {len(smiles_list)} molecules to {output_path}")

'''
#------------------
# Function to calculate drug likeliness properties
def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    
    return {
        'MW': round(mw, 2),
        'LogP': round(logp, 2),
        'HBD': hbd,
        'HBA': hba,
        'Lipinski_Violations': violations,
        'QED': round(qed(mol), 3)
    }

# Extract and process molecules
smiles_list = [mol.uid for mol in forward_network.mols]

results = []
for smi in smiles_list:
    props = calculate_properties(smi)
    if props:
        props['SMILES'] = smi
        props['Is_Starter'] = smi in user_starters
        props['Is_Helper'] = smi in user_helpers
        props['Is_Target'] = smi in user_target
        results.append(props)

df = pd.DataFrame(results)

# Create output directory
output_dir = Path("output")
output_dir.mkdir(exist_ok=True)

# Save all molecules
df.to_csv(output_dir / f"{job_name}_all_molecules.csv", index=False)

# Filter for drug-like molecules (exclude helpers/cofactors)
df_products = df[~df['Is_Helper']].copy()
df_druglike = df_products[df_products['Lipinski_Violations'] <= 1].sort_values('QED', ascending=False)
df_druglike.to_csv(output_dir / f"{job_name}_druglike_candidates.csv", index=False)

# Print summary
print(f"\n{'='*50}")
print(f"RESULTS SUMMARY")
print(f"{'='*50}")
print(f"Total molecules: {len(df)}")
print(f"Products (excluding helpers): {len(df_products)}")
print(f"Drug-like candidates: {len(df_druglike)}")
print(f"Target found: {df['Is_Target'].any()}")
print(f"\nFiles saved to: {output_dir}")
print(f"\nTop 5 drug-like candidates:")
print(df_druglike[['SMILES', 'MW', 'LogP', 'QED']].head().to_string())
#------------------
'''

post_processing.one_step(
    networks = {
        forward_network,
        },
    total_generations = 1,
    starters = user_starters,
    helpers = user_helpers,
    target = user_target,
    job_name = job_name,
    )
