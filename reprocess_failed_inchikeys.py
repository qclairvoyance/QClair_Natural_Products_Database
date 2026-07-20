import pandas as pd
import os
from datetime import datetime
from pubchem_cid_scraper import PubChemCIDExtractor
import time
import argparse


class FailedInChIKeyReprocessor:
    """Reprocess failed InChIKey entries"""
    
    def __init__(self, final_results_file, output_dir='output_retry'):
        """
        Initialize reprocessor
        
        Args:
            final_results_file (str): Path to final_results.csv
            output_dir (str): Directory for retry output files
        """
        self.final_results_file = final_results_file
        self.output_dir = output_dir
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        print(f"Loading final results from {final_results_file}...")
        self.df = pd.read_csv(final_results_file)
        
        print(f"Total records: {len(self.df):,}")
        
        # Initialize extractor
        self.extractor = PubChemCIDExtractor()
    
    def extract_failed_inchikeys(self):
        """Extract InChIKeys that failed or were not found"""
        # Find records where CID is missing/null
        failed_df = self.df[self.df['cid'].isna()].copy()
        
        print(f"\nFailed/Missing entries found: {len(failed_df):,}")
        
        if len(failed_df) == 0:
            print("No failed entries to reprocess!")
            return None
        
        # Extract unique InChIKeys from failed entries
        if 'inchikey' in failed_df.columns:
            failed_inchikeys = failed_df['inchikey'].dropna().unique()
        else:
            # Try to find the inchikey column
            inchikey_cols = [col for col in failed_df.columns if 'inchikey' in col.lower()]
            if inchikey_cols:
                failed_inchikeys = failed_df[inchikey_cols[0]].dropna().unique()
            else:
                raise ValueError("Could not find InChIKey column in the data")
        
        print(f"Unique failed InChIKeys: {len(failed_inchikeys):,}")
        
        # Save failed InChIKeys to a file for reference
        failed_list_file = os.path.join(self.output_dir, 'failed_inchikeys_list.csv')
        pd.DataFrame({'inchikey': failed_inchikeys}).to_csv(failed_list_file, index=False)
        print(f"Failed InChIKeys list saved to: {failed_list_file}")
        
        return failed_inchikeys
    
    def reprocess_failed(self, delay=0.5, save_interval=100, batch_size=5000):
        
        # Get failed InChIKeys
        failed_inchikeys = self.extract_failed_inchikeys()
        
        if failed_inchikeys is None or len(failed_inchikeys) == 0:
            return None
        
        # Check for already reprocessed entries
        checkpoint_file = os.path.join(self.output_dir, 'retry_checkpoint.csv')
        already_processed = set()
        
        if os.path.exists(checkpoint_file):
            try:
                checkpoint_df = pd.read_csv(checkpoint_file)
                if 'inchikey' in checkpoint_df.columns:
                    already_processed = set(checkpoint_df['inchikey'].dropna())
                    print(f"\nFound {len(already_processed):,} already reprocessed InChIKeys")
            except Exception as e:
                print(f"Warning: Could not read checkpoint: {e}")
        
        # Get remaining InChIKeys to process
        remaining_inchikeys = [k for k in failed_inchikeys if k not in already_processed]
        
        if len(remaining_inchikeys) == 0:
            print("\nAll failed InChIKeys already reprocessed!")
            return self.combine_retry_results()
        
        print(f"Remaining to reprocess: {len(remaining_inchikeys):,}")
        
        # Process with retry logic
        total_batches = (len(remaining_inchikeys) + batch_size - 1) // batch_size
        start_time = datetime.now()
        
        all_results = []
        if os.path.exists(checkpoint_file):
            all_results = pd.read_csv(checkpoint_file).to_dict('records')
        
        print(f"\n{'='*70}")
        print(f"STARTING REPROCESSING OF FAILED INCHIKEYS")
        print(f"{'='*70}")
        
        for batch_num in range(total_batches):
            batch_start = batch_num * batch_size
            batch_end = min((batch_num + 1) * batch_size, len(remaining_inchikeys))
            batch_keys = remaining_inchikeys[batch_start:batch_end]
            
            print(f"\nBatch {batch_num + 1}/{total_batches}")
            print(f"Processing {len(batch_keys):,} InChIKeys")
            
            for i, inchikey in enumerate(batch_keys, 1):
                # Retry with longer timeout or different approach
                result = self.extractor.get_cid_from_inchikey(inchikey)
                all_results.append(result)
                
                if i % 100 == 0:
                    print(f"  Processed {i:,}/{len(batch_keys):,} in this batch")
                
                # Save progress
                if i % save_interval == 0 or i == len(batch_keys):
                    results_df = pd.DataFrame(all_results)
                    results_df.to_csv(checkpoint_file, index=False)
                
                time.sleep(delay)
            
            # Calculate statistics for this batch
            batch_results = all_results[-len(batch_keys):]
            batch_df = pd.DataFrame(batch_results)
            success_count = len(batch_df[batch_df['status'] == 'success'])
            success_rate = (success_count / len(batch_df)) * 100
            
            print(f"  Batch success: {success_count}/{len(batch_df)} ({success_rate:.1f}%)")
        
        # Final save
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(checkpoint_file, index=False)
        
        elapsed = (datetime.now() - start_time).total_seconds()
        print(f"\n{'='*70}")
        print(f"REPROCESSING COMPLETE!")
        print(f"Time elapsed: {elapsed/3600:.2f} hours")
        print(f"{'='*70}")
        
        return self.combine_retry_results()
    
    def combine_retry_results(self):
        """Combine retry results with original data"""
        checkpoint_file = os.path.join(self.output_dir, 'retry_checkpoint.csv')
        
        if not os.path.exists(checkpoint_file):
            print("No retry checkpoint file found!")
            return None
        
        # Read retry results
        retry_df = pd.read_csv(checkpoint_file)
        
        print("\n" + "="*70)
        print("RETRY RESULTS SUMMARY")
        print("="*70)
        
        total_retried = len(retry_df)
        newly_successful = len(retry_df[retry_df['cid'].notna()])
        still_failed = len(retry_df[retry_df['cid'].isna()])
        
        print(f"Total retried: {total_retried:,}")
        print(f"Newly successful: {newly_successful:,}")
        print(f"Still failed: {still_failed:,}")
        print(f"Retry success rate: {(newly_successful/total_retried*100):.2f}%")
        
        # Update original results with retry data
        print("\nUpdating original results with retry data...")
        
        # Create a mapping of inchikey to new CID
        retry_mapping = retry_df[retry_df['cid'].notna()][['inchikey', 'cid']].set_index('inchikey')['cid'].to_dict()
        
        # Load original final results
        final_df = pd.read_csv(self.final_results_file)
        
        # Update CIDs for retried entries
        updated_count = 0
        for inchikey, new_cid in retry_mapping.items():
            mask = (final_df['inchikey'] == inchikey) & (final_df['cid'].isna())
            if mask.any():
                final_df.loc[mask, 'cid'] = new_cid
                final_df.loc[mask, 'status'] = 'success_retry'
                updated_count += 1
        
        print(f"Updated {updated_count:,} records in final results")
        
        # Save updated final results
        updated_file = os.path.join(self.output_dir, 'final_results_updated.csv')
        final_df.to_csv(updated_file, index=False)
        
        print("\n" + "="*70)
        print("UPDATED FINAL RESULTS")
        print("="*70)
        print(f"Total records: {len(final_df):,}")
        print(f"Successfully retrieved CIDs: {final_df['cid'].notna().sum():,}")
        print(f"Failed/Not found: {final_df['cid'].isna().sum():,}")
        print(f"Overall success rate: {(final_df['cid'].notna().sum()/len(final_df)*100):.2f}%")
        print(f"\nUpdated results saved to: {updated_file}")
        
        # Save still-failed entries for review
        still_failed_df = final_df[final_df['cid'].isna()]
        if len(still_failed_df) > 0:
            still_failed_file = os.path.join(self.output_dir, 'still_failed_entries.csv')
            still_failed_df.to_csv(still_failed_file, index=False)
            print(f"Still failed entries saved to: {still_failed_file}")
        
        return final_df


def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(
        description='Reprocess failed/missing InChIKeys from final results'
    )
    
    parser.add_argument('--final-results', 
                       default='/home/qc_intern/work/PubChem/output/final_results.csv',
                       help='Path to final_results.csv file')
    parser.add_argument('--output-dir', default='output_retry',
                       help='Output directory for retry results (default: output_retry)')
    parser.add_argument('--batch-size', type=int, default=5000,
                       help='Batch size for retry (default: 5000)')
    parser.add_argument('--delay', type=float, default=0.5,
                       help='Delay between requests in seconds (default: 0.5)')
    parser.add_argument('--save-interval', type=int, default=100,
                       help='Save progress every N records (default: 100)')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.final_results):
        print(f"Error: File not found: {args.final_results}")
        print("\nPlease provide the correct path to final_results.csv")
        return
    
    # Create reprocessor
    reprocessor = FailedInChIKeyReprocessor(
        final_results_file=args.final_results,
        output_dir=args.output_dir
    )
    
    # Reprocess failed entries
    reprocessor.reprocess_failed(
        delay=args.delay,
        save_interval=args.save_interval,
        batch_size=args.batch_size
    )


if __name__ == "__main__":
    print("""
    REPROCESS FAILED INCHIKEYS
    
    This script will:
    1. Read your final_results.csv
    2. Extract all failed/missing entries (where CID is null)
    3. Retry fetching PubChem CIDs for these entries
    4. Update the final results with newly found CIDs
    5. Save updated results to final_results_updated.csv
    
    USAGE EXAMPLES:
    
    # Basic usage (uses default path)
    python reprocess_failed_inchikeys.py
    
    # Specify custom final results file
    python reprocess_failed_inchikeys.py --final-results /path/to/final_results.csv
    
    # Custom output directory
    python reprocess_failed_inchikeys.py --output-dir retry_results
    
    # Faster processing (shorter delay)
    python reprocess_failed_inchikeys.py --delay 0.3
    
    # Smaller batches
    python reprocess_failed_inchikeys.py --batch-size 2000
    
    # For programmatic use:
    from reprocess_failed_inchikeys import FailedInChIKeyReprocessor
    
    reprocessor = FailedInChIKeyReprocessor('output/final_results.csv', 'output_retry')
    results = reprocessor.reprocess_failed(delay=0.5, batch_size=5000)
    """)
    
    import sys
    if len(sys.argv) == 1:
        print("\nStarting with default settings...")
        print("Press Ctrl+C to cancel, or wait 3 seconds to continue...")
        try:
            time.sleep(3)
            main()
        except KeyboardInterrupt:
            print("\n\nCancelled. Use --help for usage information.")
    else:
        main()

