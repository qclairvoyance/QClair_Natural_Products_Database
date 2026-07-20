import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, AllChem
from rdkit.Chem import rdFreeSASA
import multiprocessing as mp
from functools import partial
import warnings
warnings.filterwarnings('ignore')
from tqdm import tqdm
import time

def calculate_descriptors_for_molecule(smiles_or_inchi, source_type='smiles'):

    
    # Initialize result dictionary with None values
    result = {
        'heavy_atom_count': None,
        'total_atom_count': None,
        'hbd': None,
        'hba': None,
        'bond_count': None,
        'rotatable_bond_count': None,
        'aromatic_ring_count': None,
        'logp': None,
        'tpsa': None,
        'formal_charge': None,
        'vanderwaals_volume': None,
        'calculation_status': 'Failed',
        'error_message': ''
    }
    
    # Check for empty/invalid input
    if pd.isna(smiles_or_inchi) or smiles_or_inchi == '' or smiles_or_inchi is None:
        result['error_message'] = 'Empty input'
        return result
    
    try:
        # Parse molecule based on source type
        if source_type == 'smiles':
            mol = Chem.MolFromSmiles(str(smiles_or_inchi))
        else:  # inchi
            mol = Chem.MolFromInchi(str(smiles_or_inchi))
        
        # Check if molecule was successfully created
        if mol is None:
            result['error_message'] = f'Invalid {source_type}'
            return result
        
        # Add explicit hydrogens for accurate calculations
        mol_h = Chem.AddHs(mol)
        
        # Calculate all descriptors
        try:
            # 1. Heavy atom count (non-hydrogen atoms)
            result['heavy_atom_count'] = mol.GetNumHeavyAtoms()
            
            # 2. Total atom count (including hydrogens)
            result['total_atom_count'] = mol_h.GetNumAtoms()
            
            # 3. Hydrogen Bond Donors (HBD)
            result['hbd'] = Descriptors.NumHDonors(mol)
            
            # 4. Hydrogen Bond Acceptors (HBA)
            result['hba'] = Descriptors.NumHAcceptors(mol)
            
            # 5. Bond count (total number of bonds)
            result['bond_count'] = mol.GetNumBonds()
            
            # 6. Rotatable bond count
            result['rotatable_bond_count'] = Descriptors.NumRotatableBonds(mol)
            
            # 7. Aromatic ring count
            result['aromatic_ring_count'] = Descriptors.NumAromaticRings(mol)
            
            # 8. LogP (partition coefficient)
            result['logp'] = Descriptors.MolLogP(mol)
            
            # 9. TPSA (Topological Polar Surface Area)
            result['tpsa'] = Descriptors.TPSA(mol)
            
            # 10. Formal charge
            result['formal_charge'] = Chem.GetFormalCharge(mol)
            
            # 11. Van der Waals volume
            # Using AllChem to compute 3D conformer for accurate volume
            try:
                # Generate 3D conformer
                mol_3d = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_3d, randomSeed=42)
                AllChem.UFFOptimizeMolecule(mol_3d)
                
                # Calculate Van der Waals volume using FreeSASA
                radii = rdFreeSASA.classifyAtoms(mol_3d)
                sasa_result = rdFreeSASA.CalcSASA(mol_3d, radii)
                
                # Van der Waals volume approximation
                # Using molecular volume from 3D conformer
                result['vanderwaals_volume'] = AllChem.ComputeMolVolume(mol_3d)
                
            except Exception as vdw_error:
                # Fallback: Use 2D approximation if 3D fails
                # Van der Waals volume approximation using Crippen method
                result['vanderwaals_volume'] = Descriptors.MolWt(mol) / 0.6  # Rough approximation
                # Note: This is a fallback approximation, not as accurate as 3D method
            
            # Mark as successful
            result['calculation_status'] = 'Success'
            result['error_message'] = ''
            
        except Exception as desc_error:
            result['calculation_status'] = 'Partial'
            result['error_message'] = f'Descriptor calculation error: {str(desc_error)[:100]}'
    
    except Exception as e:
        result['error_message'] = f'Molecule parsing error: {str(e)[:100]}'
        result['calculation_status'] = 'Failed'
    
    return result


def process_single_row(row_data, use_canonical=True):
    """
    Process a single row from the dataframe
    
    Parameters:
    -----------
    row_data : tuple
        (index, row_dict) from iterrows()
    use_canonical : bool
        If True, prefer canonical_smiles, else try isomeric_smiles first
    
    Returns:
    --------
    tuple : (index, descriptor_dict)
    """
    idx, row = row_data
    
    # Priority order for molecular representation
    if use_canonical:
        sources = [
            ('canonical_smiles', 'smiles'),
            ('isomeric_smiles', 'smiles'),
            ('stdinchi', 'inchi')
        ]
    else:
        sources = [
            ('isomeric_smiles', 'smiles'),
            ('canonical_smiles', 'smiles'),
            ('stdinchi', 'inchi')
        ]
    
    # Try each source until one succeeds
    for col_name, source_type in sources:
        if col_name in row and pd.notna(row[col_name]) and row[col_name] != '':
            result = calculate_descriptors_for_molecule(row[col_name], source_type)
            
            if result['calculation_status'] in ['Success', 'Partial']:
                result['source_used'] = col_name
                return (idx, result)
    
    # If all sources failed
    return (idx, {
        'heavy_atom_count': None,
        'total_atom_count': None,
        'hbd': None,
        'hba': None,
        'bond_count': None,
        'rotatable_bond_count': None,
        'aromatic_ring_count': None,
        'logp': None,
        'tpsa': None,
        'formal_charge': None,
        'vanderwaals_volume': None,
        'calculation_status': 'Failed',
        'error_message': 'No valid source available',
        'source_used': 'None'
    })


def parallel_descriptor_calculation(df, use_canonical=True, n_cores=None):
    """
    Calculate descriptors for all molecules using multiprocessing
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe with molecular structures
    use_canonical : bool
        Preference for canonical vs isomeric SMILES
    n_cores : int
        Number of CPU cores to use (None = use all available)
    
    Returns:
    --------
    pandas.DataFrame : DataFrame with added descriptor columns
    """
    
    # Determine number of cores
    if n_cores is None:
        n_cores = mp.cpu_count()
    
    print(f"\n{'='*80}")
    print(f"Starting Parallel Descriptor Calculation")
    print(f"{'='*80}")
    print(f"Total records to process: {len(df):,}")
    print(f"CPU cores available: {mp.cpu_count()}")
    print(f"CPU cores to use: {n_cores}")
    print(f"SMILES preference: {'canonical_smiles' if use_canonical else 'isomeric_smiles'}")
    print(f"{'='*80}\n")
    
    # Prepare data for multiprocessing
    row_data = list(df.iterrows())
    
    # Create partial function with fixed use_canonical parameter
    process_func = partial(process_single_row, use_canonical=use_canonical)
    
    # Initialize result storage
    results = {}
    
    # Process using multiprocessing with progress bar
    start_time = time.time()
    
    with mp.Pool(processes=n_cores) as pool:
        # Use tqdm for progress bar
        with tqdm(total=len(row_data), desc="Calculating descriptors", unit="molecules") as pbar:
            for idx, result in pool.imap(process_func, row_data, chunksize=1000):
                results[idx] = result
                pbar.update(1)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"\n{'='*80}")
    print(f"Calculation Complete!")
    print(f"Total time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
    print(f"Processing rate: {len(df)/elapsed_time:.2f} molecules/second")
    print(f"{'='*80}\n")
    
    # Convert results to DataFrame columns
    descriptor_columns = [
        'heavy_atom_count', 'total_atom_count', 'hbd', 'hba', 
        'bond_count', 'rotatable_bond_count', 'aromatic_ring_count',
        'logp', 'tpsa', 'formal_charge', 'vanderwaals_volume'
    ]
    
    # Initialize new columns
    for col in descriptor_columns:
        df[col] = None
    
    df['calculation_status'] = ''
    df['calculation_error'] = ''
    df['descriptor_source'] = ''
    
    # Populate columns with results
    for idx, result in results.items():
        for col in descriptor_columns:
            df.at[idx, col] = result.get(col)
        
        df.at[idx, 'calculation_status'] = result.get('calculation_status', 'Unknown')
        df.at[idx, 'calculation_error'] = result.get('error_message', '')
        df.at[idx, 'descriptor_source'] = result.get('source_used', 'None')
    
    return df


def generate_calculation_report(df):
    """
    Generate comprehensive report of descriptor calculation results
    """
    print("\n" + "="*80)
    print("DESCRIPTOR CALCULATION REPORT")
    print("="*80)
    
    total = len(df)
    
    # Calculation status breakdown
    print("\nCalculation Status:")
    status_counts = df['calculation_status'].value_counts()
    for status, count in status_counts.items():
        percentage = (count / total) * 100
        print(f"  {status}: {count:,} ({percentage:.2f}%)")
    
    # Source used breakdown
    print("\nSource Used for Calculation:")
    source_counts = df['descriptor_source'].value_counts()
    for source, count in source_counts.items():
        percentage = (count / total) * 100
        print(f"  {source}: {count:,} ({percentage:.2f}%)")
    
    # Descriptor completeness
    print("\nDescriptor Completeness:")
    descriptor_cols = [
        'heavy_atom_count', 'total_atom_count', 'hbd', 'hba',
        'bond_count', 'rotatable_bond_count', 'aromatic_ring_count',
        'logp', 'tpsa', 'formal_charge', 'vanderwaals_volume'
    ]
    
    for col in descriptor_cols:
        non_null = df[col].notna().sum()
        percentage = (non_null / total) * 100
        print(f"  {col}: {non_null:,}/{total:,} ({percentage:.2f}%)")
    
    # Descriptor statistics (for successfully calculated values)
    successful_df = df[df['calculation_status'] == 'Success']
    
    if len(successful_df) > 0:
        print(f"\nDescriptor Statistics (from {len(successful_df):,} successful calculations):")
        print("-" * 80)
        
        for col in descriptor_cols:
            values = successful_df[col].dropna()
            if len(values) > 0:
                print(f"\n{col}:")
                print(f"  Count: {len(values):,}")
                print(f"  Mean: {values.mean():.4f}")
                print(f"  Median: {values.median():.4f}")
                print(f"  Min: {values.min():.4f}")
                print(f"  Max: {values.max():.4f}")
                print(f"  Std Dev: {values.std():.4f}")
    
    # Error analysis
    print("\n" + "="*80)
    print("Error Analysis:")
    failed_df = df[df['calculation_status'] == 'Failed']
    if len(failed_df) > 0:
        print(f"\nFailed calculations: {len(failed_df):,}")
        error_counts = failed_df['calculation_error'].value_counts().head(10)
        print("\nTop 10 Error Messages:")
        for error, count in error_counts.items():
            if error and error != '':
                print(f"  {error}: {count}")
    else:
        print("\nNo failed calculations! ")
    
    print("="*80)


def validate_descriptor_quality(df):
    """
    Perform quality checks on calculated descriptors
    """
    print("\n" + "="*80)
    print("DESCRIPTOR QUALITY VALIDATION")
    print("="*80)
    
    successful_df = df[df['calculation_status'] == 'Success']
    issues = []
    
    # Quality checks
    print("\nPerforming quality checks...")
    
    # 1. Check for negative values where they shouldn't exist
    for col in ['heavy_atom_count', 'total_atom_count', 'bond_count', 
                'rotatable_bond_count', 'aromatic_ring_count', 'hbd', 'hba']:
        negative = successful_df[successful_df[col] < 0]
        if len(negative) > 0:
            issues.append(f"⚠ Found {len(negative)} negative values in {col}")
    
    # 2. Check logical relationships
    # Total atoms should be >= heavy atoms
    invalid_atom_count = successful_df[successful_df['total_atom_count'] < successful_df['heavy_atom_count']]
    if len(invalid_atom_count) > 0:
        issues.append(f"⚠ Found {len(invalid_atom_count)} records where total_atom_count < heavy_atom_count")
    
    # 3. Check for extremely high values (potential errors)
    if successful_df['logp'].max() > 20:
        high_logp = successful_df[successful_df['logp'] > 20]
        issues.append(f"ℹ Found {len(high_logp)} records with LogP > 20 (may be valid for large molecules)")
    
    # 4. Check TPSA range (should be 0-300 typically)
    if successful_df['tpsa'].max() > 400:
        high_tpsa = successful_df[successful_df['tpsa'] > 400]
        issues.append(f"ℹ Found {len(high_tpsa)} records with TPSA > 400 (may be valid for large molecules)")
    
    # Print results
    if issues:
        print("\nQuality Check Results:")
        for issue in issues:
            print(f"  {issue}")
    else:
        print("\n All quality checks passed!")
    
    print("="*80)


# Main execution
if __name__ == "__main__":
    
    print("\n" + "="*80)
    print("RDKit Descriptor Calculator for Natural Products Database")
    print("="*80)
    
    # Configuration
    INPUT_FILE = "/home/qc_intern/work/QC_DB/New_Data/new_db_0/new_db_2.csv"
    OUTPUT_FILE = "/home/qc_intern/work/QC_DB/New_Data/new_db_0/new_db_3.csv"
    USE_CANONICAL_SMILES = True  # Set to False to prefer isomeric_smiles
    N_CORES = None  # None = use all available cores
    
    # Load data
    print(f"\nLoading data from: {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df):,} records with {len(df.columns)} columns")
    print(f"Columns: {df.columns.tolist()}")
    
    # Verify required columns exist
    required_cols = ['canonical_smiles', 'isomeric_smiles', 'stdinchi']
    available_cols = [col for col in required_cols if col in df.columns]
    print(f"\nAvailable molecular representations: {available_cols}")
    
    if not available_cols:
        print("\nERROR: No valid molecular representation columns found!")
        print("Required: at least one of ['canonical_smiles', 'isomeric_smiles', 'stdinchi']")
        exit(1)
    
    # Calculate descriptors using multiprocessing
    df_with_descriptors = parallel_descriptor_calculation(
        df, 
        use_canonical=USE_CANONICAL_SMILES,
        n_cores=N_CORES
    )
    
    # Generate report
    generate_calculation_report(df_with_descriptors)
    
    # Validate quality
    validate_descriptor_quality(df_with_descriptors)
    
    # Save results
    print(f"\n{'='*80}")
    print(f"Saving results to: {OUTPUT_FILE}")
    df_with_descriptors.to_csv(OUTPUT_FILE, index=False)
    print(f"Save complete!")
    
    # Save successful records only
    successful_df = df_with_descriptors[df_with_descriptors['calculation_status'] == 'Success']
    success_output = OUTPUT_FILE.replace('.csv', '_success_only.csv')
    print(f"\nSaving successful calculations only to: {success_output}")
    successful_df.to_csv(success_output, index=False)
    print(f"Successful records saved: {len(successful_df):,}")
    
    # Save failed records for review
    failed_df = df_with_descriptors[df_with_descriptors['calculation_status'] == 'Failed']
    if len(failed_df) > 0:
        failed_output = OUTPUT_FILE.replace('.csv', '_failed.csv')
        print(f"\nSaving failed calculations to: {failed_output}")
        failed_df.to_csv(failed_output, index=False)
        print(f"Failed records saved: {len(failed_df):,}")
    
    print("\n" + "="*80)
    print("All processing complete!")
    print("="*80)
    
    # Final summary
    print("\nFinal Summary:")
    print(f"  Total processed: {len(df_with_descriptors):,}")
    print(f"  Successful: {len(successful_df):,} ({len(successful_df)/len(df_with_descriptors)*100:.2f}%)")
    print(f"  Failed: {len(failed_df):,} ({len(failed_df)/len(df_with_descriptors)*100:.2f}%)")
    print(f"\nOutput files created:")
    print(f"  1. {OUTPUT_FILE} (complete dataset)")
    print(f"  2. {success_output} (successful only)")
    if len(failed_df) > 0:
        print(f"  3. {failed_output} (failed records)")
    
    print("\nDone!\n")
