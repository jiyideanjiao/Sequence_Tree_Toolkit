import pandas as pd
import sys
import os

def process_ogg_data_pandas(input_csv, output_csv, threshold=0.20, output_dir=None):
    try:
        # Read the CSV into a pandas DataFrame
        df = pd.read_csv(input_csv)
        
        # Check for required columns
        required_columns = {'ogg_id', 'species', 'length', 'aa', 'gap'}
        if not required_columns.issubset(df.columns):
            missing = required_columns - set(df.columns)
            print("Error: Missing columns in input CSV: {}".format(', '.join(missing)))
            sys.exit(1)
        
        # Calculate gap_rate
        df['gap_rate'] = df.apply(lambda row: round(row['gap'] / row['aa'], 4) if row['aa'] != 0 else 'NA', axis=1)
        
        # Filter rows with gap_rate < threshold and gap_rate is not 'NA'
        filtered_df = df[(df['gap_rate'] < threshold) & (df['gap_rate'] != 'NA')]
        
        # Save to output CSV
        filtered_df.to_csv(output_csv, index=False)
        print("Gap rate calculation and filtering completed. Output written to '{}'.".format(output_csv))
        
        # Generate ogg_id text files
        if output_dir is None:
            output_dir = os.path.dirname(os.path.abspath(output_csv))
        else:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        
        # Group by ogg_id and write species lists
        grouped = filtered_df.groupby('ogg_id')['species']
        for ogg_id, species_series in grouped:
            txt_filename = os.path.join(output_dir, "{}.txt".format(ogg_id))
            species_series.to_csv(txt_filename, index=False, header=False)
            print("Generated file: {}".format(txt_filename))
    
    except FileNotFoundError:
        print("Error: The file '{}' does not exist.".format(input_csv))
    except PermissionError:
        print("Error: Permission denied when accessing '{}' or '{}'.".format(input_csv, output_csv))
    except Exception as e:
        print("An unexpected error occurred: {}".format(e))

def main():
    if len(sys.argv) not in [3, 4, 5]:
        print("Usage: python filter_trimal_OG.py <input_csv> <output_csv> [threshold] [output_dir]")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_csv = sys.argv[2]
    
    # Set default values
    threshold = 0.20
    output_dir = None
    
    if len(sys.argv) >= 4:
        try:
            threshold = float(sys.argv[3])
            if not (0 <= threshold <= 1):
                print("Error: Threshold must be between 0 and 1.")
                sys.exit(1)
        except ValueError:
            print("Error: Threshold must be a numeric value.")
            sys.exit(1)
    
    if len(sys.argv) == 5:
        output_dir = sys.argv[4]
    
    process_ogg_data_pandas(input_csv, output_csv, threshold, output_dir)

if __name__ == "__main__":
    main()
