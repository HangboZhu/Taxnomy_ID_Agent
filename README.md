# Taxonomy Agent

A Large Language Model-powered batch processing tool for species taxonomy information, designed to convert common species names to Latin scientific names and retrieve NCBI Taxonomy IDs.

## Features

- 🤖 **Intelligent Name Conversion**: Uses Zhipu GLM-4.5 model for bidirectional conversion between common and Latin names
- 📊 **Batch Data Processing**: Supports efficient batch processing of CSV files with large species datasets
- 🔄 **Dual Workflow System**: Primary workflow + fallback workflow ensuring high success rates
- 📚 **NCBI Database Integration**: Queries authoritative taxonomy IDs through ETE3 library
- 🛡️ **Robust Error Handling**: Comprehensive retry mechanisms and exception handling

## How It Works

### Core Algorithm Flow

1. **Primary Workflow**: Common name → Latin name → NCBI TaxID
2. **Fallback Workflow**: Latin name → Common name → Latin name → NCBI TaxID

### Intelligent Decision Logic

The project employs a rule-based decision system rather than relying on LLM for decision making:
- Prioritizes common name conversion (primary workflow)
- Enables Latin name-based fallback mechanism when common names are invalid or conversion fails
- Ensures data accuracy and completeness through multi-layer validation

## Installation

```bash
# Install dependencies using uv (recommended)
uv install

# Or using pip
pip install -e .
```

## Environment Configuration

Create a `.env` file and add your Zhipu API key:

```env
ZHIPU_API_KEY=your_zhipu_api_key_here
```

## Usage

### Basic Usage

```bash
# Use default file paths
python taxnomy_agent.py

# Use custom file paths
python taxnomy_agent.py -i input.csv -o output.csv -d /path/to/taxdump.tar.gz

# View all options
python taxnomy_agent.py --help
```

### Parameters

- `-i, --input`: Input CSV file path (default: `./Host_Range_output.csv`)
- `-o, --output`: Output CSV file path (default: `./Host_Range_output_update.csv`)
- `-d, --cachedir`: ETE3 cache file path (optional, default: `./NCBI_taxnomy_db_dir/taxdump.tar.gz`)

### Input Data Format

CSV files should contain the following columns (script will automatically create missing columns):
- `Common Name (Host)`: Species common names
- `Latin name (Host)`: Latin scientific names
- `Taxonomy ID (Host)`: NCBI Taxonomy IDs

## NCBI Taxonomy Database

The system uses the NCBI Taxonomy database for authoritative species information. You can manually download the database from:

**Download URL**: https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

The first run of ETE3 will automatically download ~100MB of taxonomy data to `~/.etetoolkit/` or use the custom cache directory specified with `-d` parameter.

## Core Functions

### Name Conversion Functions
- `common_name_to_latin()`: Convert common names to Latin names
- `latin_to_common_name()`: Convert Latin names to common names

### Batch Processing Functions
- `batch_common_to_latin()`: Batch common name conversion
- `batch_latin_to_taxid_ete3()`: Batch TaxID querying

### Data Processing Function
- `process_taxonomy_csv()`: Main CSV processing workflow

## Project Structure

```
taxnomy_agent/
├── taxnomy_agent.py          # Main script file
├── final.ipynb              # Development and debugging notebook
├── data/                    # Data folder
│   ├── sample_Host_Range_output.csv
│   └── output_test_taxonomy_updated.csv
├── NCBI_taxnomy_db_dir/     # NCBI taxonomy database cache
├── pyproject.toml           # Project configuration
├── uv.lock                  # Dependency lock file
└── .env                     # Environment variables (create manually)
```

## Technical Features

- **Smart Retry Mechanism**: Automatic retry on API call failures (up to 3 attempts)
- **Data Cleaning**: Automatic handling of special characters and null values
- **Encoding Compatibility**: Support for UTF-8 and latin1 encodings
- **Progress Display**: Real-time progress tracking with tqdm
- **Cache Optimization**: Local NCBI database caching after first run

## Important Notes

- First ETE3 run downloads ~100MB NCBI taxonomy database to `~/.etetoolkit/`
- Ensure stable internet connection for Zhipu API and NCBI database access
- Recommend testing with small datasets first before processing large-scale data

## License

This project is open source. Please refer to the project configuration files for license details.