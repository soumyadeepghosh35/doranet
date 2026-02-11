#!/usr/bin/env python

'''
### Compare against a natural products database and use matches as targets
'''

#!/usr/bin/env python

import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
import csv
import time

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


def mol_to_canonical_smiles(mol) -> str:
    """Convert RDKit mol object to canonical SMILES."""
    try:
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


def mol_to_inchikey(mol) -> str:
    """Convert RDKit mol object to InChIKey."""
    try:
        if mol:
            inchi = MolToInchi(mol)
            if inchi:
                return InchiToInchiKey(inchi)
    except:
        pass
    return None


def get_mol_property(mol, property_names: list) -> str:
    """Get first available property from molecule."""
    if not mol:
        return None
    
    props = mol.GetPropsAsDict()
    
    for prop_name in property_names:
        if prop_name in props:
            val = props[prop_name]
            if val and str(val).strip():
                return str(val).strip()
    
    return None


def load_natural_products_from_sdf(sdf_path: str, database_name: str = "unknown") -> dict:
    """
    Load natural products database from SDF file.
    
    Args:
        sdf_path: Path to SDF file
        database_name: Name of database ('npatlas', 'coconut', or 'unknown')
    
    Returns:
        dict mapping canonical SMILES/InChIKey -> molecule info
    """
    print(f"Loading natural products database from: {sdf_path}")
    print(f"Database type: {database_name}")
    
    # Property names for different databases
    name_properties = {
        'npatlas': ['compound_names', 'npaid', 'compound_name', 'name', 'ID'],
        'coconut': ['name', 'coconut_id', 'Name', 'ID', 'compound_name'],
        'unknown': ['name', 'Name', 'compound_name', 'compound_names', 'ID', 'id']
    }
    
    id_properties = {
        'npatlas': ['npaid', 'npatlas_id', 'ID', 'id'],
        'coconut': ['coconut_id', 'ID', 'id'],
        'unknown': ['ID', 'id', 'compound_id']
    }
    
    # Additional metadata properties
    source_properties = {
        'npatlas': ['origin_type', 'origin_species', 'genus', 'species'],
        'coconut': ['source', 'organism', 'textTaxa'],
        'unknown': ['source', 'origin', 'organism']
    }
    
    db_key = database_name.lower() if database_name.lower() in name_properties else 'unknown'
    
    # Load SDF file
    supplier = Chem.SDMolSupplier(sdf_path)
    
    if supplier is None:
        raise ValueError(f"Could not load SDF file: {sdf_path}")
    
    np_database = {}
    valid_count = 0
    error_count = 0
    
    # First pass: check available properties
    sample_mol = None
    for mol in supplier:
        if mol is not None:
            sample_mol = mol
            break
    
    if sample_mol:
        available_props = list(sample_mol.GetPropsAsDict().keys())
        print(f"  Available properties: {available_props[:15]}{'...' if len(available_props) > 15 else ''}")
    
    # Reset supplier
    supplier = Chem.SDMolSupplier(sdf_path)
    
    print("  Processing molecules...")
    
    for idx, mol in enumerate(supplier):
        if mol is None:
            error_count += 1
            continue
        
        # Get canonical SMILES directly from molecule
        canonical = mol_to_canonical_smiles(mol)
        if not canonical:
            error_count += 1
            continue
        
        # Get InChIKey
        inchikey = mol_to_inchikey(mol)
        
        # Get name
        name = get_mol_property(mol, name_properties[db_key])
        
        # Get ID
        mol_id = get_mol_property(mol, id_properties[db_key])
        
        # Get source/organism info
        source = get_mol_property(mol, source_properties[db_key])
        
        # Get molecular formula
        try:
            formula = rdMolDescriptors.CalcMolFormula(mol)
        except:
            formula = None
        
        valid_count += 1
        
        entry = {
            'canonical_smiles': canonical,
            'inchikey': inchikey,
            'name': name,
            'id': mol_id,
            'source': source,
            'formula': formula,
            'database': database_name
        }
        
        # Index by canonical SMILES
        np_database[canonical] = entry
        
        # Also index by InChIKey for robust matching
        if inchikey:
            np_database[inchikey] = entry
        
        # Progress update
        if (idx + 1) % 50000 == 0:
            print(f"    Processed {idx + 1} molecules...")
    
    print(f"  Loaded {valid_count} valid natural products")
    print(f"  Skipped {error_count} invalid entries")
    
    return np_database


def load_multiple_np_databases(sdf_paths: list) -> dict:
    """
    Load multiple natural products databases.
    
    Args:
        sdf_paths: List of tuples (sdf_path, database_name)
                   e.g., [('NPAtlas_np_data.sdf', 'npatlas'), 
                          ('coconut_np_data.sdf', 'coconut')]
    
    Returns:
        Combined dictionary of all natural products
    """
    combined_database = {}
    
    for sdf_path, db_name in sdf_paths:
        print(f"\n{'='*50}")
        db = load_natural_products_from_sdf(sdf_path, db_name)
        
        # Merge databases
        for key, entry in db.items():
            if key not in combined_database:
                combined_database[key] = entry
            else:
                # Entry exists, append database info
                existing = combined_database[key]
                if existing['database'] != entry['database']:
                    existing['database'] = f"{existing['database']}, {entry['database']}"
    
    print(f"\n{'='*50}")
    print(f"Total unique compounds: {len(combined_database) // 2}")  # Divided by 2 since we index by both SMILES and InChIKey
    
    return combined_database


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
    np_databases: list,
    job_name: str,
    generation: int = 2,
    max_atoms: dict = None,
    helpers: set = None
):
    """
    Run DORAnet network generation, match against natural products,
    and generate pathways to matches.
    
    Args:
        starter_smiles: Starting molecule SMILES
        np_databases: List of tuples [(sdf_path, db_name), ...]
        job_name: Job name for output files
        generation: Number of generations
        max_atoms: Max atom constraints
        helpers: Helper molecules
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
    
    # Step 1: Load natural products databases
    print("\n" + "-" * 70)
    print("Step 1: Loading Natural Products Databases")
    print("-" * 70)
    
    np_database = load_multiple_np_databases(np_databases)
    
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
            db = info.get('database', 'Unknown')
            np_id = info.get('id', 'N/A')
            print(f"  {i+1}. {smi}")
            print(f"      Name: {name}, ID: {np_id}, Database: {db}")
        if len(np_matches) > 20:
            print(f"  ... and {len(np_matches) - 20} more matches")
    
    # Step 4: Save results
    print("\n" + "-" * 70)
    print("Step 4: Saving Results")
    print("-" * 70)
    
    # Save all molecules
    all_molecules_path = Path(f"{job_name}_all_molecules.csv")
    with open(all_molecules_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['SMILES', 'Is_Starter', 'Is_NaturalProduct', 'NP_Name', 
                        'NP_ID', 'NP_Database', 'NP_Source', 'MolFormula', 
                        'MolWeight', 'NumHeavyAtoms'])
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                mol_weight = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
                num_heavy = mol.GetNumHeavyAtoms()
            else:
                formula, mol_weight, num_heavy = "N/A", 0, 0
            
            is_np = smi in np_matches
            if is_np:
                info = np_matches[smi]
                np_name = info.get('name', '')
                np_id = info.get('id', '')
                np_db = info.get('database', '')
                np_source = info.get('source', '')
            else:
                np_name, np_id, np_db, np_source = '', '', '', ''
            
            writer.writerow([smi, smi in user_starters, is_np, np_name, np_id,
                           np_db, np_source, formula, mol_weight, num_heavy])
    
    print(f"Saved all molecules to: {all_molecules_path}")
    
    # Save NP matches separately
    if np_matches:
        np_matches_path = Path(f"{job_name}_np_matches.csv")
        with open(np_matches_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Generated_SMILES', 'NP_Name', 'NP_ID', 'NP_Database',
                           'NP_Source', 'InChIKey', 'MolFormula', 'MolWeight'])
            for smi, info in np_matches.items():
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    formula = rdMolDescriptors.CalcMolFormula(mol)
                    mol_weight = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
                else:
                    formula, mol_weight = "N/A", 0
                
                writer.writerow([smi, info.get('name', ''), info.get('id', ''),
                               info.get('database', ''), info.get('source', ''),
                               info.get('inchikey', ''), formula, mol_weight])
        
        print(f"Saved NP matches to: {np_matches_path}")
    
    # Step 5: Generate pathways
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
    job_name = "NP_SDF"
    generation = 2
    
    # Natural Products Databases (SDF format)
    # List of tuples: (path_to_sdf, database_name)
    np_databases = [
        ('../natural_products_data/NPAtlas_np_data.sdf', 'npatlas'),
        ('../natural_products_data/coconut_np_data.sdf', 'coconut')
    ]
    
    max_atoms = {
        'C': 12,
        'N': 5,
        'O': 8,
        'S': 3
    }
    
    results = run_doranet_with_np_matching(
        starter_smiles=starter_smiles,
        np_databases=np_databases,
        job_name=job_name,
        generation=generation,
        max_atoms=max_atoms
    )
    print(f"Time: {time.time() - start_time:.2f}")
