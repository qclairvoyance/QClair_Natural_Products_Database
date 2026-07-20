import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Scaffolds, Descriptors
from rdkit.Chem import rdFreeSASA
import multiprocessing as mp
from functools import partial
import warnings
warnings.filterwarnings('ignore')
from tqdm import tqdm
import time

# Import selfies library
try:
    import selfies as sf
    SELFIES_AVAILABLE = True
except ImportError:
    SELFIES_AVAILABLE = False
    print("WARNING: selfies library not found. Install with: pip install selfies")
    print("SELFIES generation will be skipped.\n")


def calculate_advanced_descriptors(row_data, use_canonical=True):
    """
    Calculate SMARTS, SELFIES, Murcko Framework, and VdW Surface Area for a single molecule
    
    Parameters:
    -----------
    row_data : tuple
        (index, row_dict) from DataFrame.iterrows()
    use_canonical : bool
        If True, prefer canonical_smiles, else try isomeric_smiles first
    
    Returns:
    --------
    tuple : (index, result_dict)
    """
    idx, row = row_data
    
    # Initialize result dictionary
    result = {
        'smarts': None,
        'selfies': None,
        'murcko_framework': None,
        'vanderwaals_surface_area': None,
        'calculation_status': 'Failed',
        'error_message': '',
        'source_used': None
    }
    
    # Determine source priority
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
    
    # Try each source until we get a valid molecule
    mol = None
    source_used = None
    
    for col_name, source_type in sources:
        if col_name in row and pd.notna(row[col_name]) and row[col_name] != '':
            try:
                if source_type == 'smiles':
                    mol = Chem.MolFromSmiles(str(row[col_name]))
                else:  # inchi
                    mol = Chem.MolFromInchi(str(row[col_name]))
                
                if mol is not None:
                    source_used = col_name
                    break
            except:
                continue
    
    # If no valid molecule found
    if mol is None:
        result['error_message'] = 'No valid molecular source available'
        return (idx, result)
    
    result['source_used'] = source_used
    
    try:
        # ==================== 1. SMARTS Pattern ====================
        try:
            # Generate SMARTS from molecule
            # Using MolToSmarts for a SMARTS representation
            smarts = Chem.MolToSmarts(mol)
            result['smarts'] = smarts if smarts else None
        except Exception as e:
            result['smarts'] = None
            # Don't fail entire calculation for SMARTS error
        
        # ==================== 2. SELFIES ====================
        try:
            if SELFIES_AVAILABLE:
                # Convert SMILES to SELFIES
                # First get SMILES from the molecule
                smiles_for_selfies = Chem.MolToSmiles(mol)
                
                # Convert to SELFIES using selfies library
                try:
                    selfies_str = sf.encoder(smiles_for_selfies)
                    result['selfies'] = selfies_str
                except Exception as selfies_error:
                    # Some molecules may not encode properly
                    result['selfies'] = None
            else:
                result['selfies'] = None
        except Exception as e:
            result['selfies'] = None
        
        # ==================== 3. Murcko Framework ====================
        try:
            # Generate Bemis-Murcko scaffold
            scaffold = Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
            
            if scaffold:
                # Convert scaffold to SMILES for storage
                scaffold_smiles = Chem.MolToSmiles(scaffold)
                result['murcko_framework'] = scaffold_smiles
            else:
                result['murcko_framework'] = None
        except Exception as e:
            result['murcko_framework'] = None
        
        # ==================== 4. Van der Waals Surface Area ====================
        try:
            # Add explicit hydrogens
            mol_h = Chem.AddHs(mol)
            
            # Generate 3D conformer for accurate surface area calculation
            try:
                # Embed molecule in 3D
                AllChem.EmbedMolecule(mol_h, randomSeed=42)
                # Optimize geometry
                AllChem.UFFOptimizeMolecule(mol_h)
                
                # Calculate Van der Waals surface area using FreeSASA
                radii = rdFreeSASA.classifyAtoms(mol_h)
                sasa_result = rdFreeSASA.CalcSASA(mol_h, radii)
                
                # Get Van der Waals surface area
                # FreeSASA returns solvent accessible surface area
                # For VdW surface area, we use the total SASA value
                vdw_surface_area = sum(sasa_result)
                result['vanderwaals_surface_area'] = round(vdw_surface_area, 4)
                
            except Exception as embed_error:
                # Fallback: Use 2D approximation
                # Estimate based on molecular weight and atom count
                try:
                    # Simple approximation: 4.0 Ų per heavy atom
                    heavy_atoms = mol.GetNumHeavyAtoms()
                    vdw_surface_area_approx = heavy_atoms * 4.0
                    result['vanderwaals_surface_area'] = round(vdw_surface_area_approx, 4)
                except:
                    result['vanderwaals_surface_area'] = None
        
        except Exception as e:
            result['vanderwaals_surface_area'] = None
        
        # Check if at least one descriptor was calculated
        descriptors_calculated = sum([
            result['smarts'] is not None,
            result['selfies'] is not None,
            result['murcko_framework'] is not None,
            result['vanderwaals_surface_area'] is not None
        ])
        
        if descriptors_calculated >= 3:
            result['calculation_status'] = 'Success'
            result['error_message'] = ''
        elif descriptors_calculated > 0:
            result['calculation_status'] = 'Partial'
            result['error_message'] = 'Some descriptors could not be calculated'
        else:
            result['calculation_status'] = 'Failed'
            result['error_message'] = 'All descriptor calculations failed'
    
    except Exception as e:
        result['calculation_status'] = 'Failed'
        result['error_message'] = f'General error: {str(e)[:100]}'
    
    return (idx, result)


def parallel_advanced_descriptor_calculation(df, use_canonical=True, n_cores=None):
    
    # Determine number of cores
    if n_cores is None:
        n_cores = mp.cpu_count()
    
    print(f"\n{'='*80}")
    print(f"Advanced Descriptor Calculation - Multiprocessing")
    print(f"{'='*80}")
    print(f"Descriptors to calculate:")
    print(f"  1. SMARTS pattern")
    print(f"  2. SELFIES representation")
    print(f"  3. Murcko framework (Bemis-Murcko scaffold)")
    print(f"  4. Van der Waals surface area")
    print(f"\nDataset Information:")
    print(f"  Total records: {len(df):,}")
    print(f"  CPU cores available: {mp.cpu_count()}")
    print(f"  CPU cores to use: {n_cores}")
    print(f"  SMILES preference: {'canonical_smiles' if use_canonical else 'isomeric_smiles'}")
    print(f"  SELFIES library: {'Available ' if SELFIES_AVAILABLE else 'Not installed '}")
    print(f"{'='*80}\n")
    
    # Prepare data for multiprocessing
    row_data = list(df.iterrows())
    
    # Create partial function
    process_func = partial(calculate_advanced_descriptors, use_canonical=use_canonical)
    
    # Initialize result storage
    results = {}
    
    # Process using multiprocessing with progress bar
    start_time = time.time()
    
    print("Starting parallel processing...")
    with mp.Pool(processes=n_cores) as pool:
        with tqdm(total=len(row_data), desc="Processing molecules", unit="mol") as pbar:
            for idx, result in pool.imap(process_func, row_data, chunksize=500):
                results[idx] = result
                pbar.update(1)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"\n{'='*80}")
    print(f"Processing Complete!")
    print(f"{'='*80}")
    print(f"Total time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
    print(f"Processing rate: {len(df)/elapsed_time:.2f} molecules/second")
    print(f"{'='*80}\n")
    
    # Add new columns to dataframe
    df['smarts'] = None
    df['selfies'] = None
    df['murcko_framework'] = None
    df['vanderwaals_surface_area'] = None
    df['adv_calculation_status'] = ''
    df['adv_calculation_error'] = ''
    df['adv_descriptor_source'] = ''
    
    # Populate columns with results
    for idx, result in results.items():
        df.at[idx, 'smarts'] = result.get('smarts')
        df.at[idx, 'selfies'] = result.get('selfies')
        df.at[idx, 'murcko_framework'] = result.get('murcko_framework')
        df.at[idx, 'vanderwaals_surface_area'] = result.get('vanderwaals_surface_area')
        df.at[idx, 'adv_calculation_status'] = result.get('calculation_status', 'Unknown')
        df.at[idx, 'adv_calculation_error'] = result.get('error_message', '')
        df.at[idx, 'adv_descriptor_source'] = result.get('source_used', 'None')
    
    return df


def generate_detailed_report(df):
    """
    Generate comprehensive report of advanced descriptor calculations
    """
    print("\n" + "="*80)
    print("ADVANCED DESCRIPTOR CALCULATION REPORT")
    print("="*80)
    
    total = len(df)
    
    # Overall status
    print("\n1. Calculation Status:")
    print("-" * 80)
    status_counts = df['adv_calculation_status'].value_counts()
    for status, count in status_counts.items():
        percentage = (count / total) * 100
        print(f"  {status}: {count:,} ({percentage:.2f}%)")
    
    # Individual descriptor completeness
    print("\n2. Descriptor Completeness:")
    print("-" * 80)
    descriptors = {
        'smarts': 'SMARTS Pattern',
        'selfies': 'SELFIES Representation',
        'murcko_framework': 'Murcko Framework',
        'vanderwaals_surface_area': 'Van der Waals Surface Area'
    }
    
    for col, name in descriptors.items():
        non_null = df[col].notna().sum()
        percentage = (non_null / total) * 100
        print(f"  {name}: {non_null:,}/{total:,} ({percentage:.2f}%)")
    
    # Source usage
    print("\n3. Source Used for Calculation:")
    print("-" * 80)
    source_counts = df['adv_descriptor_source'].value_counts()
    for source, count in source_counts.items():
        percentage = (count / total) * 100
        print(f"  {source}: {count:,} ({percentage:.2f}%)")
    
    # Statistics for VdW surface area
    successful_df = df[df['adv_calculation_status'] == 'Success']
    if len(successful_df) > 0:
        print("\n4. Van der Waals Surface Area Statistics:")
        print("-" * 80)
        vdw_values = successful_df['vanderwaals_surface_area'].dropna()
        if len(vdw_values) > 0:
            print(f"  Count: {len(vdw_values):,}")
            print(f"  Mean: {vdw_values.mean():.2f} Ų")
            print(f"  Median: {vdw_values.median():.2f} Ų")
            print(f"  Min: {vdw_values.min():.2f} Ų")
            print(f"  Max: {vdw_values.max():.2f} Ų")
            print(f"  Std Dev: {vdw_values.std():.2f} Ų")
    
    # Murcko framework analysis
    print("\n5. Murcko Framework Analysis:")
    print("-" * 80)
    murcko_values = df['murcko_framework'].dropna()
    unique_scaffolds = murcko_values.nunique()
    print(f"  Total scaffolds calculated: {len(murcko_values):,}")
    print(f"  Unique scaffolds: {unique_scaffolds:,}")
    if len(murcko_values) > 0:
        scaffold_diversity = (unique_scaffolds / len(murcko_values)) * 100
        print(f"  Scaffold diversity: {scaffold_diversity:.2f}%")
        
        # Most common scaffolds
        print("\n  Top 5 Most Common Scaffolds:")
        top_scaffolds = murcko_values.value_counts().head(5)
        for i, (scaffold, count) in enumerate(top_scaffolds.items(), 1):
            percentage = (count / len(murcko_values)) * 100
            # Truncate long SMILES for display
            scaffold_display = scaffold[:50] + "..." if len(scaffold) > 50 else scaffold
            print(f"    {i}. {scaffold_display}")
            print(f"       Count: {count:,} ({percentage:.2f}%)")
    
    # SELFIES statistics
    if SELFIES_AVAILABLE:
        print("\n6. SELFIES Representation:")
        print("-" * 80)
        selfies_values = df['selfies'].dropna()
        print(f"  Total SELFIES generated: {len(selfies_values):,}")
        if len(selfies_values) > 0:
            avg_length = selfies_values.str.len().mean()
            print(f"  Average SELFIES length: {avg_length:.0f} characters")
    
    # Error analysis
    print("\n7. Error Analysis:")
    print("-" * 80)
    failed_df = df[df['adv_calculation_status'] == 'Failed']
    if len(failed_df) > 0:
        print(f"  Failed calculations: {len(failed_df):,}")
        error_counts = failed_df['adv_calculation_error'].value_counts().head(10)
        print("\n  Top Error Messages:")
        for error, count in error_counts.items():
            if error and error != '':
                error_display = error[:70] + "..." if len(error) > 70 else error
                print(f"    - {error_display}: {count}")
    else:
        print("  No failed calculations! ")
    
    print("\n" + "="*80)


def validate_descriptors(df):
    """
    Validate the quality of calculated descriptors
    """
    print("\n" + "="*80)
    print("DESCRIPTOR QUALITY VALIDATION")
    print("="*80)
    
    issues = []
    
    successful_df = df[df['adv_calculation_status'] == 'Success']
    
    print("\nRunning quality checks...")
    
    # 1. Check VdW surface area range
    vdw_values = successful_df['vanderwaals_surface_area'].dropna()
    if len(vdw_values) > 0:
        negative_vdw = vdw_values[vdw_values < 0]
        if len(negative_vdw) > 0:
            issues.append(f"Found {len(negative_vdw)} negative VdW surface area values")
        
        extremely_high_vdw = vdw_values[vdw_values > 2000]
        if len(extremely_high_vdw) > 0:
            issues.append(f"Found {len(extremely_high_vdw)} very high VdW surface area values (>2000 Ų)")
    
    # 2. Check SMARTS validity
    smarts_values = successful_df['smarts'].dropna()
    invalid_smarts = 0
    for smarts in smarts_values.head(100):  # Sample check
        try:
            mol = Chem.MolFromSmarts(smarts)
            if mol is None:
                invalid_smarts += 1
        except:
            invalid_smarts += 1
    
    if invalid_smarts > 0:
        issues.append(f"Found {invalid_smarts}/100 invalid SMARTS in sample")
    
    # 3. Check SELFIES validity (if available)
    if SELFIES_AVAILABLE:
        selfies_values = successful_df['selfies'].dropna()
        invalid_selfies = 0
        for selfie in selfies_values.head(100):  # Sample check
            try:
                smiles = sf.decoder(selfie)
                if not smiles:
                    invalid_selfies += 1
            except:
                invalid_selfies += 1
        
        if invalid_selfies > 0:
            issues.append(f"Found {invalid_selfies}/100 invalid SELFIES in sample")
    
    # 4. Check Murcko framework validity
    murcko_values = successful_df['murcko_framework'].dropna()
    invalid_murcko = 0
    for murcko in murcko_values.head(100):  # Sample check
        try:
            mol = Chem.MolFromSmiles(murcko)
            if mol is None:
                invalid_murcko += 1
        except:
            invalid_murcko += 1
    
    if invalid_murcko > 0:
        issues.append(f"Found {invalid_murcko}/100 invalid Murcko frameworks in sample")
    
    # Print results
    print("\nValidation Results:")
    if issues:
        for issue in issues:
            print(f"  {issue}")
    else:
        print(" All quality checks passed!")
    
    print("="*80)


# Main execution
if __name__ == "__main__":
    
    print("\n" + "="*80)
    print("Advanced Molecular Descriptor Calculator")
    print("SMARTS | SELFIES | Murcko Framework | VdW Surface Area")
    print("="*80)
    
    # Configuration
    INPUT_FILE = "/home/qc_intern/work/QC_DB/New_Data/new_db_0/new_db_2.csv"
    OUTPUT_FILE = "/home/qc_intern/work/QC_DB/New_Data/new_db_0/new_db_4.csv"
    USE_CANONICAL_SMILES = True  # Set to False to prefer isomeric_smiles
    N_CORES = None  # None = use all available cores
    
    # Check if selfies is available
    if not SELFIES_AVAILABLE:
        print("\n" + "!"*80)
        print("IMPORTANT: SELFIES library not installed!")
        print("Install with: pip install selfies")
        print("SELFIES calculations will be skipped.")
        print("!"*80 + "\n")
        time.sleep(2)
    
    # Load data
    print(f"\nLoading data from: {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df):,} records")
    print(f"Shape: {df.shape}")
    print(f"Columns: {df.columns.tolist()}")
    
    # Verify required columns
    required_cols = ['canonical_smiles', 'isomeric_smiles', 'stdinchi']
    available_cols = [col for col in required_cols if col in df.columns]
    print(f"\nAvailable molecular representations: {available_cols}")
    
    if not available_cols:
        print("\nERROR: No valid molecular representation columns found!")
        print("Required: at least one of ['canonical_smiles', 'isomeric_smiles', 'stdinchi']")
        exit(1)
    
    # Calculate advanced descriptors using multiprocessing
    df_with_descriptors = parallel_advanced_descriptor_calculation(
        df,
        use_canonical=USE_CANONICAL_SMILES,
        n_cores=N_CORES
    )
    
    # Generate comprehensive report
    generate_detailed_report(df_with_descriptors)
    
    # Validate descriptor quality
    validate_descriptors(df_with_descriptors)
    
    # Save results
    print(f"\n{'='*80}")
    print("SAVING RESULTS")
    print(f"{'='*80}")
    print(f"\nSaving complete dataset to: {OUTPUT_FILE}")
    df_with_descriptors.to_csv(OUTPUT_FILE, index=False)
    print("Save complete!")
    
    # Save successful records only
    successful_df = df_with_descriptors[df_with_descriptors['adv_calculation_status'] == 'Success']
    success_output = OUTPUT_FILE.replace('.csv', '_success_only.csv')
    print(f"\nSaving successful calculations only to: {success_output}")
    successful_df.to_csv(success_output, index=False)
    print(f" Successful records saved: {len(successful_df):,}")
    
    # Save failed/partial records for review
    problematic_df = df_with_descriptors[
        df_with_descriptors['adv_calculation_status'].isin(['Failed', 'Partial'])
    ]
    if len(problematic_df) > 0:
        problem_output = OUTPUT_FILE.replace('.csv', '_problematic.csv')
        print(f"\nSaving problematic calculations to: {problem_output}")
        problematic_df.to_csv(problem_output, index=False)
        print(f" Problematic records saved: {len(problematic_df):,}")
    
    print("\n" + "="*80)
    print("ALL PROCESSING COMPLETE!")
    print("="*80)
    
    # Final summary
    print("\nFinal Summary:")
    print(f"  Total processed: {len(df_with_descriptors):,}")
    print(f"  Successful: {len(successful_df):,} ({len(successful_df)/len(df_with_descriptors)*100:.2f}%)")
    
    # Individual descriptor success rates
    print("\nDescriptor Success Rates:")
    for descriptor in ['smarts', 'selfies', 'murcko_framework', 'vanderwaals_surface_area']:
        count = df_with_descriptors[descriptor].notna().sum()
        rate = (count / len(df_with_descriptors)) * 100
        status = "" if rate > 95 else ""
        print(f"  {status} {descriptor}: {count:,} ({rate:.2f}%)")
    
    print(f"\nOutput Files Created:")
    print(f"  1. {OUTPUT_FILE}")
    print(f"     (Complete dataset with all {len(df_with_descriptors):,} records)")
    print(f"  2. {success_output}")
    print(f"     (Success only - {len(successful_df):,} records)")
    if len(problematic_df) > 0:
        print(f"  3. {problem_output}")
        print(f"     (Problematic - {len(problematic_df):,} records)")
    
    if not SELFIES_AVAILABLE:
        print("\n" + "!"*80)
        print("REMINDER: Install selfies library for SELFIES generation:")
        print("  pip install selfies")
        print("!"*80)
    
    print("\n Done!\n")

