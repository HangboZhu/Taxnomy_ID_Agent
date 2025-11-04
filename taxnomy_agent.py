from zai import ZhipuAiClient
import requests
from dotenv import load_dotenv
import os
import time
from ete3 import NCBITaxa
import traceback
import pandas as pd
from tqdm import tqdm
import numpy as np
import re
import argparse  # Added for argparse

# --- API Related And Data Cleaning Function ---

def clean_name(name: str) -> str:
    """(Helper) Clean name, remove * and strip whitespace"""
    if pd.isna(name):
        return np.nan
    cleaned_name = str(name).replace('*', '').strip()
    return cleaned_name

def is_invalid(name: str) -> bool:
    """(Helper) Check if name is invalid (nan, '', '-', 'nan')"""
    if pd.isna(name):
        return True
    return str(name).strip().lower() in ['', '-', 'nan']

def common_name_to_latin(zhipu_api_key: str, common_name: str) -> str:
    """(Helper) Use Zhipu GLM to convert a single common name to a Latin name"""
    try:
        client = ZhipuAiClient(api_key=zhipu_api_key)
        prompt = f"""Strictly follow 2 rules, no extra output:
1. Convert the species common name to Latin name, formatted as "Genus (capitalized first letter) + space + species epithet (lowercase)" (e.g., "domestic cat" → "Felis catus");
2. If the common name is unrecognizable, return only "unrecognizable".

Species common name: {common_name}"""
        
        response = client.chat.completions.create(
            model="glm-4.5",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,
            max_tokens=50,
            thinking={"type": "disabled"}
        )
        latin_name = response.choices[0].message.content.strip()
        return latin_name if latin_name else "unrecognizable"
    except Exception as e:
        return f"GLM call failed: {str(e)}"

def latin_to_common_name(zhipu_api_key: str, latin_name: str) -> str:
    """(NEW Helper) Use Zhipu GLM to convert a single Latin name to a common name"""
    try:
        client = ZhipuAiClient(api_key=zhipu_api_key)
        prompt = f"""Strictly follow 2 rules, no extra output:
1. Convert the species Latin name to its most common English name (e.g., "Felis catus" → "domestic cat");
2. If the Latin name is unrecognizable, return only "unrecognizable".

Latin name: {latin_name}"""
        
        response = client.chat.completions.create(
            model="glm-4.5",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1,
            max_tokens=50,
            thinking={"type": "disabled"}
        )
        common_name = response.choices[0].message.content.strip()
        return common_name if common_name else "unrecognizable"
    except Exception as e:
        return f"GLM call failed: {str(e)}"

# --- Batch Process Function ---

def batch_common_to_latin(zhipu_api_key: str, common_names: list[str]) -> dict[str, str]:
    """(Batch) Batch convert list of common names to Latin names"""
    if not common_names: return {}
    print(f"Batch converting {len(common_names)} common names...")
    latin_map = {}
    for name in tqdm(common_names, desc="Converting common names"):
        retries = 3
        for i in range(retries):
            latin_name = common_name_to_latin(zhipu_api_key, str(name))
            if "GLM call failed" not in latin_name:
                latin_map[name] = latin_name
                break
            time.sleep(2 * (i + 1))
        else:
            print(f"Error: GLM call failed for '{name}' after {retries} attempts.")
            latin_map[name] = "GLM call failed"
    return latin_map

def batch_latin_to_common(zhipu_api_key: str, latin_names: list[str]) -> dict[str, str]:
    """(Batch) Batch convert list of Latin names to common names"""
    if not latin_names: return {}
    print(f"Batch reverse converting {len(latin_names)} Latin names...")
    common_map = {}
    for name in tqdm(latin_names, desc="Reverse converting Latin names"):
        retries = 3
        for i in range(retries):
            common_name = latin_to_common_name(zhipu_api_key, str(name))
            if "GLM call failed" not in common_name:
                common_map[name] = common_name
                break
            time.sleep(2 * (i + 1))
        else:
            print(f"Error: GLM call failed for '{name}' after {retries} attempts.")
            common_map[name] = "GLM call failed"
    return common_map

def batch_latin_to_taxid_ete3(latin_names: list[str], cache_dir: str = None) -> dict[str, str]:
    """(Batch) Batch query Taxonomy IDs for a list of Latin names (using ete3)"""
    if not latin_names:
        return {}
        
    print(f"Batch querying TaxIDs for {len(latin_names)} Latin names...")
    print(f"Note: Initializing NCBI Taxonomy client...")
    
    taxid_map = {}
    try:
        ncbi = NCBITaxa(taxdump_file=cache_dir)
        raw_taxid_dict = ncbi.get_name_translator(latin_names)
        print("NCBI raw query data received.")

        for latin_name in latin_names:
            if latin_name not in raw_taxid_dict:
                taxid_map[latin_name] = f"Not found: '{latin_name}'"
            else:
                taxid_list = raw_taxid_dict[latin_name]
                taxid_map[latin_name] = str(taxid_list[0])
        
        print("TaxID query complete.")
        return taxid_map
    except Exception as e:
        error_msg = str(e).lower()
        error_str = f"Error: {str(e)}"
        if any(keyword in error_msg for keyword in ["cache", "taxdump", "sqlite"]):
            error_str = f"Cache error: {str(e)}"
        elif any(keyword in error_msg for keyword in ["connection", "network", "download"]):
            error_str = f"Network error: {str(e)}"
        print(f"A critical error occurred during TaxID query: {error_str}")
        for name in latin_names:
            if name not in taxid_map:
                taxid_map[name] = error_str
        return taxid_map

# --- Refactored CSV Processor ---

def process_taxonomy_csv(input_path: str, output_path: str, zhipu_api_key: str, ete3_cache_dir: str = None):
    """
    Reads a CSV, executes primary and fallback workflows, batch updates taxonomy info,
    and saves to a new CSV.
    """
    print("=" * 50)
    print(f"Starting processing for file: {input_path}")
    print("=" * 50)

    try:
        try:
            df = pd.read_csv(input_path)
        except UnicodeDecodeError:
            print("UTF-8 decoding failed, trying 'latin1' encoding...")
            df = pd.read_csv(input_path, encoding='latin1')
        except Exception as e:
            print(f"Error loading CSV: {e}")
            return

        df.columns = df.columns.str.strip()
        
        common_col = "Common Name"
        latin_col = "Latin name"
        taxid_col = "Taxonomy ID"

        for col in [common_col, latin_col, taxid_col]:
            if col not in df.columns:
                print(f"Warning: Column '{col}' not found. It will be created.")
                df[col] = np.nan

        print("Step 1/5: Cleaning Common Name and Latin Name columns...")
        df[common_col] = df[common_col].apply(clean_name)
        df[latin_col] = df[latin_col].apply(clean_name)

        common_is_invalid = df[common_col].apply(is_invalid)
        latin_is_invalid = df[latin_col].apply(is_invalid)

        mask_primary = ~common_is_invalid
        mask_fallback = common_is_invalid & ~latin_is_invalid

        print("Initializing temporary buffers...")
        df['temp_common'] = pd.Series(dtype='object')
        df['temp_latin'] = pd.Series(dtype='object')
        df['temp_taxid'] = pd.Series(dtype='object')

        print(f"\nStep 2/5: Executing primary workflow (based on Common Name)... {mask_primary.sum()} rows")
        if mask_primary.sum() > 0:
            unique_common_names = df.loc[mask_primary, common_col].unique()
            
            latin_map = batch_common_to_latin(zhipu_api_key, unique_common_names.tolist())
            
            unique_latin_names = list(set(
                name for name in latin_map.values() 
                if not is_invalid(name) and "failed" not in name and "unrecognizable" not in name
            ))
            taxid_map = batch_latin_to_taxid_ete3(unique_latin_names, ete3_cache_dir)
            
            df.loc[mask_primary, 'temp_latin'] = df.loc[mask_primary, common_col].map(latin_map)
            df.loc[mask_primary, 'temp_taxid'] = df.loc[mask_primary, 'temp_latin'].map(taxid_map)

        print(f"\nStep 3/5: Executing fallback workflow (based on Latin Name)... {mask_fallback.sum()} rows")
        if mask_fallback.sum() > 0:
            unique_latin_fallback = df.loc[mask_fallback, latin_col].unique().tolist()
            
            print("Fallback(a): Batch reverse querying Common Name (for completion)...")
            reverse_map = batch_latin_to_common(zhipu_api_key, unique_latin_fallback) # latin -> common
            
            # Attempt to complete common name regardless of subsequent success
            df.loc[mask_fallback, 'temp_common'] = df.loc[mask_fallback, latin_col].map(reverse_map)

            print("Fallback(b): Attempting direct TaxID query...")
            taxid_map_attempt1 = batch_latin_to_taxid_ete3(unique_latin_fallback, ete3_cache_dir)
            
            failed_latin_names = [
                name for name, taxid in taxid_map_attempt1.items()
                if "Not found" in str(taxid) or "Error" in str(taxid)
            ]
            
            # Map results from Attempt 1
            df.loc[mask_fallback, 'temp_latin'] = df.loc[mask_fallback, latin_col] # Latin name is the original one
            df.loc[mask_fallback, 'temp_taxid'] = df.loc[mask_fallback, 'temp_latin'].map(taxid_map_attempt1)

            if failed_latin_names:
                print(f"Fallback(c): Starting Attempt 2 for {len(failed_latin_names)} failed Latin names...")
                
                # We get the new common names from the populated 'temp_common' column
                mask_attempt2 = mask_fallback & df[latin_col].isin(failed_latin_names)

                # Get the new common names corresponding to these rows
                new_common_names = df.loc[mask_attempt2, 'temp_common'].unique().tolist()
                
                # Filter out invalid values
                new_common_names = [
                    name for name in new_common_names
                    if not is_invalid(name) and "failed" not in name and "unrecognizable" not in name
                ]

                forward_map = batch_common_to_latin(zhipu_api_key, new_common_names) # common -> new_latin
                
                new_latin_names = list(set(
                    name for name in forward_map.values()
                    if not is_invalid(name) and "failed" not in name and "unrecognizable" not in name
                ))
                taxid_map_attempt2 = batch_latin_to_taxid_ete3(new_latin_names, ete3_cache_dir)
                
                # Map results from Attempt 2 (only mapping latin and taxid, not temp_common)
                df.loc[mask_attempt2, 'temp_latin'] = df.loc[mask_attempt2, 'temp_common'].map(forward_map)
                df.loc[mask_attempt2, 'temp_taxid'] = df.loc[mask_attempt2, 'temp_latin'].map(taxid_map_attempt2)

        print("\nStep 4/5: Merging all workflow results...")
        
        df['temp_taxid'] = df['temp_taxid'].fillna(df['temp_latin'].map({
            "unrecognizable": "unrecognizable",
            "GLM call failed": "GLM call failed"
        }))
        
        # This merge logic is now correct
        # 1. Complete Common Name (not handled by primary workflow, handled by fallback[a])
        df[common_col] = df['temp_common'].fillna(df[common_col])
        
        # 2. Update Latin Name (from primary or fallback workflow)
        df[latin_col] = df['temp_latin'].fillna(df[latin_col])
        
        # 3. Update TaxID (from primary or fallback workflow)
        df[taxid_col] = df['temp_taxid'].fillna(df[taxid_col])

        df = df.drop(columns=['temp_common', 'temp_latin', 'temp_taxid'])

        df.to_csv(output_path, index=False)
        
        print("\nStep 5/5: Processing complete.")
        print(f"Updated CSV saved to: {output_path}")
        print("=" * 50)

    except Exception as e:
        print("\n" + "=" * 50)
        print(f"An unexpected error occurred during CSV processing:")
        traceback.print_exc()
        print("=" * 50)

# --- Main execution entry point ---

if __name__ == "__main__":
    
    # --- Setup Argparse ---
    parser = argparse.ArgumentParser(description="Process taxonomy data in a CSV file.")
    parser.add_argument(
        '-i', '--input', 
        help='Path to the input CSV file.', 
        default="./Host_Range_output.csv"
    )
    parser.add_argument(
        '-o', '--output', 
        help='Path to save the updated output CSV file.', 
        default="./Host_Range_output_update.csv"
    )
    parser.add_argument(
        '-d', '--cachedir', 
        help='Path to the ETE3 cache file (e.g., ./NCBI_taxnomy_db_dir/taxdump.tar.gz).', 
        default="./NCBI_taxnomy_db_dir/taxdump.tar.gz"
    )
    args = parser.parse_args()

    # --- Load API Key ---
    load_dotenv()
    YOUR_ZHIPU_API_KEY = os.getenv("ZHIPU_API_KEY")
    
    if not YOUR_ZHIPU_API_KEY:
        print("Error: ZHIPU_API_KEY not found in environment variables.")
        print("Please create a .env file and add ZHIPU_API_KEY='your_key'")
    else:
        # --- Use args from parser ---
        INPUT_CSV_PATH = args.input
        OUTPUT_CSV_PATH = args.output
        CUSTOM_ETE3_CACHE_DIR = args.cachedir

        print(f"Input file: {INPUT_CSV_PATH}")
        print(f"Output file: {OUTPUT_CSV_PATH}")
        if CUSTOM_ETE3_CACHE_DIR:
            print(f"ETE3 Cache file: {CUSTOM_ETE3_CACHE_DIR}")
        
        if not os.path.exists(INPUT_CSV_PATH):
            print(f"Error: Input file not found: {INPUT_CSV_PATH}")
            print("Please create the file and populate it with data before running.")
        else:
            process_taxonomy_csv(
                input_path=INPUT_CSV_PATH,
                output_path=OUTPUT_CSV_PATH,
                zhipu_api_key=YOUR_ZHIPU_API_KEY,
                ete3_cache_dir=CUSTOM_ETE3_CACHE_DIR
            )