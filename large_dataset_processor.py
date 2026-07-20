import pandas as pd
import os
from datetime import datetime
from pubchem_cid_scraper import PubChemCIDExtractor
import time
import argparse


class LargeDatasetProcessor:
    """Optimized processor for large InChIKey datasets"""
    
    def __init__(self, input_file, inchikey_column='inchikey', output_dir='output'):
    
        self.input_file = input_file
        self.inchikey_column = inchikey_column
        self.output_dir = output_dir
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Load data
        print(f"Loading data from {input_file}...")
        self.df = pd.read_csv(input_file)
        
        if inchikey_column not in self.df.columns:
            raise ValueError(f"Column '{inchikey_column}' not found in CSV")
        
        self.total_records = len(self.df)
        self.unique_inchikeys = self.df[inchikey_column].dropna().unique()
        
        print(f"Total records: {self.total_records:,}")
        print(f"Unique InChIKeys: {len(self.unique_inchikeys):,}")
        
        # Initialize extractor
        self.extractor = PubChemCIDExtractor()
    
    def get_processed_inchikeys(self):
        """Get list of already processed InChIKeys from checkpoint files"""
        processed = set()
        
        checkpoint_files = [f for f in os.listdir(self.output_dir) 
                          if f.startswith('checkpoint_') and f.endswith('.csv')]
        
        for file in checkpoint_files:
            try:
                df_checkpoint = pd.read_csv(os.path.join(self.output_dir, file))
                if 'inchikey' in df_checkpoint.columns:
                    processed.update(df_checkpoint['inchikey'].dropna().unique())
            except Exception as e:
                print(f"Warning: Could not read {file}: {e}")
        
        return processed
    
    def process_in_batches(self, batch_size=10000, delay=0.5, resume=True, save_interval=100):
        """
        Process InChIKeys in batches with continuous checkpoint saving
        
        Args:
            batch_size (int): Number of InChIKeys per batch
            delay (float): Delay between requests (seconds)
            resume (bool): Resume from checkpoints if available
            save_interval (int): Save progress every N records within batch
        """
        # Get already processed InChIKeys if resume is enabled
        processed_keys = set()
        if resume:
            processed_keys = self.get_processed_inchikeys()
            print(f"\nFound {len(processed_keys):,} already processed InChIKeys")
        
        # Get remaining InChIKeys
        remaining_keys = [k for k in self.unique_inchikeys if k not in processed_keys]
        
        if not remaining_keys:
            print("\nAll InChIKeys already processed!")
            return self.combine_results()
        
        print(f"Remaining to process: {len(remaining_keys):,}")
        
        # Calculate total batches needed
        total_inchikeys = len(self.unique_inchikeys)
        already_processed = len(processed_keys)
        current_batch_num = (already_processed // batch_size) + 1
        total_batches = (total_inchikeys + batch_size - 1) // batch_size
        
        print(f"Resuming from batch {current_batch_num}/{total_batches}")
        
        # Process in batches
        start_time = datetime.now()
        remaining_batches = total_batches - current_batch_num + 1
        
        for batch_idx in range(remaining_batches):
            batch_num = current_batch_num + batch_idx
            batch_start_idx = batch_idx * batch_size
            batch_end_idx = min((batch_idx + 1) * batch_size, len(remaining_keys))
            batch_keys = remaining_keys[batch_start_idx:batch_end_idx]
            
            checkpoint_file = os.path.join(
                self.output_dir,
                f'checkpoint_batch_{batch_num:04d}.csv'
            )
            
            print(f"\n{'='*70}")
            print(f"Batch {batch_num}/{total_batches}")
            print(f"Processing {len(batch_keys):,} InChIKeys in this batch")
            print(f"Checkpoint file: {checkpoint_file}")
            print(f"{'='*70}")
            
            # Check if this batch already has partial progress
            batch_results = []
            start_from = 0
            
            if resume and os.path.exists(checkpoint_file):
                try:
                    existing_df = pd.read_csv(checkpoint_file)
                    if 'inchikey' in existing_df.columns:
                        existing_keys = set(existing_df['inchikey'].dropna())
                        # Filter out already processed keys from this batch
                        batch_keys = [k for k in batch_keys if k not in existing_keys]
                        batch_results = existing_df.to_dict('records')
                        start_from = len(batch_results)
                        print(f"  Resuming batch {batch_num}: {start_from} already done, {len(batch_keys)} remaining")
                except Exception as e:
                    print(f"  Warning: Could not read existing checkpoint {checkpoint_file}: {e}")
                    print(f"  Starting batch from beginning")
            
            # Process remaining keys in this batch
            for i, inchikey in enumerate(batch_keys, 1):
                result = self.extractor.get_cid_from_inchikey(inchikey)
                batch_results.append(result)
                
                current_count = start_from + i
                
                if i % 100 == 0:
                    print(f"  Processed {current_count:,}/{start_from + len(batch_keys):,} in this batch")
                
                # Save progress incrementally
                if i % save_interval == 0 or i == len(batch_keys):
                    batch_df = pd.DataFrame(batch_results)
                    batch_df.to_csv(checkpoint_file, index=False)
                
                time.sleep(delay)
            
            # Final save of batch results
            batch_df = pd.DataFrame(batch_results)
            batch_df.to_csv(checkpoint_file, index=False)
            
            # Calculate statistics
            success_count = len(batch_df[batch_df['status'] == 'success'])
            success_rate = (success_count / len(batch_df)) * 100 if len(batch_df) > 0 else 0
            
            # Calculate time estimates
            elapsed = (datetime.now() - start_time).total_seconds()
            avg_time_per_batch = elapsed / (batch_idx + 1)
            batches_remaining = remaining_batches - (batch_idx + 1)
            eta_seconds = avg_time_per_batch * batches_remaining
            eta = datetime.now() + pd.Timedelta(seconds=eta_seconds)
            
            total_processed = already_processed + ((batch_idx + 1) * batch_size)
            
            print(f"\n  Batch Results:")
            print(f"    Successful: {success_count}/{len(batch_df)} ({success_rate:.1f}%)")
            print(f"    Saved to: {checkpoint_file}")
            print(f"\n  Progress:")
            print(f"    Completed: {batch_num}/{total_batches} batches")
            print(f"    Total processed so far: {min(total_processed, total_inchikeys):,}/{total_inchikeys:,}")
            print(f"    Elapsed time: {elapsed/3600:.2f} hours")
            if batches_remaining > 0:
                print(f"    ETA: {eta.strftime('%Y-%m-%d %H:%M:%S')}")
        
        print("\n" + "="*70)
        print("BATCH PROCESSING COMPLETE!")
        print("="*70)
        
        # Combine all results
        return self.combine_results()
    
    def combine_results(self):
        """Combine all checkpoint files into final result"""
        print("\nCombining all checkpoint files...")
        
        checkpoint_files = sorted([
            os.path.join(self.output_dir, f)
            for f in os.listdir(self.output_dir)
            if f.startswith('checkpoint_') and f.endswith('.csv')
        ])
        
        if not checkpoint_files:
            print("No checkpoint files found!")
            return None
        
        # Read and combine all checkpoints
        all_results = []
        for file in checkpoint_files:
            df = pd.read_csv(file)
            all_results.append(df)
        
        results_df = pd.concat(all_results, ignore_index=True)
        
        # Remove duplicates (keep last occurrence)
        results_df = results_df.drop_duplicates(subset=['inchikey'], keep='last')
        
        # Merge with original data
        final_df = self.df.merge(
            results_df,
            left_on=self.inchikey_column,
            right_on='inchikey',
            how='left'
        )
        
        # Save final results
        output_file = os.path.join(self.output_dir, 'final_results.csv')
        final_df.to_csv(output_file, index=False)
        
        # Print summary
        print("\n" + "="*70)
        print("FINAL RESULTS SUMMARY")
        print("="*70)
        print(f"Total records: {len(final_df):,}")
        print(f"Successfully retrieved CIDs: {final_df['cid'].notna().sum():,}")
        print(f"Failed/Not found: {final_df['cid'].isna().sum():,}")
        print(f"Success rate: {(final_df['cid'].notna().sum()/len(final_df)*100):.2f}%")
        print(f"\nFinal results saved to: {output_file}")
        
        # Save summary statistics
        summary = {
            'total_records': len(final_df),
            'successful': final_df['cid'].notna().sum(),
            'failed': final_df['cid'].isna().sum(),
            'success_rate': (final_df['cid'].notna().sum()/len(final_df)*100),
            'timestamp': datetime.now().isoformat()
        }
        
        summary_file = os.path.join(self.output_dir, 'summary.txt')
        with open(summary_file, 'w') as f:
            for key, value in summary.items():
                f.write(f"{key}: {value}\n")
        
        return final_df


def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(
        description='Process large dataset of InChIKeys to extract PubChem CIDs'
    )
    
    parser.add_argument('input_file', help='Input CSV file with InChIKeys')
    parser.add_argument('--column', default='inchikey',
                       help='Column name containing InChIKeys (default: inchikey)')
    parser.add_argument('--output-dir', default='output',
                       help='Output directory (default: output)')
    parser.add_argument('--batch-size', type=int, default=10000,
                       help='Batch size (default: 10000)')
    parser.add_argument('--delay', type=float, default=0.5,
                       help='Delay between requests in seconds (default: 0.5)')
    parser.add_argument('--save-interval', type=int, default=100,
                       help='Save progress every N records within batch (default: 100)')
    parser.add_argument('--no-resume', action='store_true',
                       help='Start fresh, ignore existing checkpoints')
    
    args = parser.parse_args()
    
    # Create processor
    processor = LargeDatasetProcessor(
        input_file=args.input_file,
        inchikey_column=args.column,
        output_dir=args.output_dir
    )
    
    # Process data
    processor.process_in_batches(
        batch_size=args.batch_size,
        delay=args.delay,
        resume=not args.no_resume,
        save_interval=args.save_interval
    )


if __name__ == "__main__":

    
    # If run without arguments, show usage
    import sys
    if len(sys.argv) == 1:
        print("\nNo input file provided. Use --help for usage information.")
    else:
        main()

