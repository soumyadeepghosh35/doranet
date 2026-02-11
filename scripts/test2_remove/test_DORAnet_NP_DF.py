#!/usr/bin/env python
"""
Extended debugging script to find why pathway counts differ.
"""

import time
import sys
import warnings
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
from rdkit import RDLogger
import pandas as pd

#!/usr/bin/env python

import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.inchi import MolFromInchi, MolToInchi, InchiToInchiKey
import csv
import pandas as pd

# Add the local doranet path
DORANET_PATH = Path("/users/sghosh6/DTRA_project/MACAW/doranet")
sys.path.insert(0, str(DORANET_PATH))

import doranet.modules.enzymatic as enzymatic
import doranet.modules.post_processing as post_processing


def canonicalize_smiles(smi: str) -> str:
    """Convert SMILES to canonical form for comparison."""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
    except:
        pass
    return None


def smiles_to_inchikey(smi: str) -> str:
    """Convert SMILES to InChIKey for robust comparison."""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            inchi = MolToInchi(mol)
            if inchi:
                return InchiToInchiKey(inchi)
    except:
        pass
    return None


def load_natural_products_database(db_path: str, smiles_column: str = None) -> dict:
    """
    Load natural products database.
    
    Returns dict mapping canonical SMILES -> original SMILES and metadata
    """
    print(f"Loading natural products database from: {db_path}")
    
    # Determine file type
    if db_path.endswith('.csv'):
        df = pd.read_csv(db_path, low_memory=False)
    elif db_path.endswith('.tsv'):
        df = pd.read_csv(db_path, sep='\t', low_memory=False)
    elif db_path.endswith('.parquet'):
        df = pd.read_parquet(db_path)
    else:
        df = pd.read_csv(db_path, low_memory=False)
    
    print(f"  Database has {len(df)} entries")
    print(f"  Columns: {df.columns.tolist()}")
    
    # Find SMILES column
    if smiles_column and smiles_column in df.columns:
        smi_col = smiles_column
    else:
        possible_cols = ['SMILES', 'smiles', 'Canonical_SMILES', 'canonical_smiles', 
                         'sugar_free_smiles', 'SUGAR_FREE_SMILES', 'molecule_smiles']
        smi_col = None
        for col in possible_cols:
            if col in df.columns:
                smi_col = col
                break
        if smi_col is None:
            raise ValueError(f"No SMILES column found. Available: {df.columns.tolist()}")
    
    print(f"  Using SMILES column: {smi_col}")
    
    # Find name/ID column if available
    name_col = None
    for col in ['name', 'Name', 'compound_name', 'ID', 'id', 'coconut_id', 'npatlas_id']:
        if col in df.columns:
            name_col = col
            break
    
    # Build lookup dictionary with InChIKey for robust matching
    np_database = {}
    valid_count = 0
    
    for idx, row in df.iterrows():
        original_smi = row[smi_col]
        if pd.isna(original_smi):
            continue
        
        canonical = canonicalize_smiles(str(original_smi))
        inchikey = smiles_to_inchikey(str(original_smi))
        
        if canonical:
            valid_count += 1
            entry = {
                'original_smiles': str(original_smi),
                'canonical_smiles': canonical,
                'inchikey': inchikey,
                'name': row[name_col] if name_col and not pd.isna(row.get(name_col)) else None
            }
            
            # Index by both canonical SMILES and InChIKey
            np_database[canonical] = entry
            if inchikey:
                np_database[inchikey] = entry
    
    print(f"  Loaded {valid_count} valid natural products")
    
    return np_database


def find_natural_product_matches(generated_smiles: list, np_database: dict) -> dict:
    """
    Find which generated molecules match natural products.
    
    Returns dict: generated_smiles -> np_info
    """
    matches = {}
    
    for smi in generated_smiles:
        canonical = canonicalize_smiles(smi)
        inchikey = smiles_to_inchikey(smi)
        
        # Try canonical SMILES match first
        if canonical and canonical in np_database:
            matches[smi] = np_database[canonical]
            continue
        
        # Try InChIKey match (more robust)
        if inchikey and inchikey in np_database:
            matches[smi] = np_database[inchikey]
    
    return matches


def run_doranet_with_np_matching(
    starter_smiles: str,
    np_database_path: str,
    job_name: str,
    generation: int = 2,
    max_atoms: dict = None,
    helpers: set = None,
    np_smiles_column: str = None
):
    """
    Run DORAnet network generation, match against natural products,
    and generate pathways to matches.
    """
    
    if max_atoms is None:
        max_atoms = {'C': 12, 'N': 5, 'O': 8, 'S': 3}
    
    if helpers is None:
        helpers = {
            'O', 'O=O', '[H][H]', 'O=C=O', 'C=O', '[C-]#[O+]', 'Br', '[Br][Br]', 'CO',
            'C=C', 'O=S(O)O', 'N', 'O=S(=O)(O)O', 'O=NO', 'N#N', 'O=[N+]([O-])O', 'NO',
            'C#N', 'S', 'O=S=O', 'N#CO'
        }
    
    user_starters = {starter_smiles}
    
    print("=" * 70)
    print("DORAnet Natural Products Pathway Finder")
    print("=" * 70)
    print(f"Starter: {starter_smiles}")
    print(f"Generation: {generation}")
    print(f"Max atoms: {max_atoms}")
    
    # Step 1: Load natural products database
    print("\n" + "-" * 70)
    print("Step 1: Loading Natural Products Database")
    print("-" * 70)
    np_database = load_natural_products_database(np_database_path, np_smiles_column)
    
    # Step 2: Generate network
    print("\n" + "-" * 70)
    print("Step 2: Generating Reaction Network")
    print("-" * 70)
    
    forward_network = enzymatic.generate_network(
        job_name=job_name,
        starters=user_starters,
        gen=generation,
        max_atoms=max_atoms,
        direction="forward",
        ruleset="JN3604IMT"
    )
    
    smiles_list = [mol.uid for mol in forward_network.mols]
    print(f"Generated {len(smiles_list)} molecules")
    
    # Step 3: Find natural product matches
    print("\n" + "-" * 70)
    print("Step 3: Matching Against Natural Products")
    print("-" * 70)
    
    np_matches = find_natural_product_matches(smiles_list, np_database)
    
    print(f"Found {len(np_matches)} natural product matches!")
    
    if np_matches:
        print("\nMatched Natural Products:")
        for i, (smi, info) in enumerate(list(np_matches.items())[:20]):
            mol = Chem.MolFromSmiles(smi)
            formula = rdMolDescriptors.CalcMolFormula(mol) if mol else "N/A"
            name = info.get('name', 'Unknown')
            print(f"  {i+1}. {smi}")
            print(f"      Formula: {formula}, Name: {name}")
        if len(np_matches) > 20:
            print(f"  ... and {len(np_matches) - 20} more matches")
    
    # Step 4: Save all molecules
    print("\n" + "-" * 70)
    print("Step 4: Saving Results")
    print("-" * 70)
    
    # Save all molecules with NP annotation
    all_molecules_path = Path(f"{job_name}_all_molecules.csv")
    with open(all_molecules_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['SMILES', 'Is_Starter', 'Is_NaturalProduct', 'NP_Name', 
                         'MolFormula', 'MolWeight', 'NumHeavyAtoms'])
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                mol_weight = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
                num_heavy = mol.GetNumHeavyAtoms()
            else:
                formula, mol_weight, num_heavy = "N/A", 0, 0
            
            is_np = smi in np_matches
            np_name = np_matches[smi].get('name', '') if is_np else ''
            
            writer.writerow([smi, smi in user_starters, is_np, np_name,
                           formula, mol_weight, num_heavy])
    
    print(f"Saved all molecules to: {all_molecules_path}")
    
    # Save NP matches separately
    if np_matches:
        np_matches_path = Path(f"{job_name}_np_matches.csv")
        with open(np_matches_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Generated_SMILES', 'NP_Original_SMILES', 'NP_Name', 
                            'InChIKey', 'MolFormula', 'MolWeight'])
            for smi, info in np_matches.items():
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    formula = rdMolDescriptors.CalcMolFormula(mol)
                    mol_weight = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
                else:
                    formula, mol_weight = "N/A", 0
                
                writer.writerow([smi, info.get('original_smiles', ''), 
                               info.get('name', ''), info.get('inchikey', ''),
                               formula, mol_weight])
        
        print(f"Saved NP matches to: {np_matches_path}")
    
    # Step 5: Generate pathways to natural products
    print("\n" + "-" * 70)
    print("Step 5: Generating Pathways to Natural Products")
    print("-" * 70)
    
    if np_matches:
        np_targets = set(np_matches.keys())
        print(f"Finding pathways to {len(np_targets)} natural product targets...")
        
        try:
            post_processing.one_step(
                networks={forward_network},
                total_generations=generation,
                starters=user_starters,
                helpers=helpers,
                target=np_targets,
                job_name=job_name,
            )
            print(f"Pathway files generated with prefix: {job_name}")
        except Exception as e:
            print(f"Warning: Pathway generation failed: {e}")
    else:
        print("No natural product matches found - skipping pathway generation")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Starter molecule: {starter_smiles}")
    print(f"Generations: {generation}")
    print(f"Total molecules generated: {len(smiles_list)}")
    print(f"Natural product matches: {len(np_matches)}")
    print(f"\nOutput files:")
    print(f"  - {job_name}_all_molecules.csv")
    if np_matches:
        print(f"  - {job_name}_np_matches.csv")
        print(f"  - {job_name}_* (pathway files)")
    
    return {
        'smiles_list': smiles_list,
        'np_matches': np_matches,
        'network': forward_network
    }


# MAIN
if __name__ == "__main__":
    
    start_time = time.time()
    
    # Configuration
    starter_smiles = 'CCC=C1OC(=O)[C@@H](C)[C@H]1O'
    job_name = "NP_DF"
    generation = 2
    
    # Path to natural products database
    # Common options:
    # - COCONUT: https://coconut.naturalproducts.net/download
    # - NPAtlas: https://www.npatlas.org/download
    # - Your custom database
    np_database_path = "../natural_products_data/naturalProductDF.csv"
    
    max_atoms = {
        'C': 12,
        'N': 5,
        'O': 8,
        'S': 3
    }
    
    results = run_doranet_with_np_matching(
        starter_smiles=starter_smiles,
        np_database_path=np_database_path,
        job_name=job_name,
        generation=generation,
        max_atoms=max_atoms,
        np_smiles_column=None  # Auto-detect, or specify like 'sugar_free_smiles'
    )

    print(f"Time: {time.time() - start_time:.2f}")