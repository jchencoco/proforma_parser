# ProForma Parser

Parses ProForma strings, maps proteoforms to UniProt sequences, and calculates theoretical masses. Results are saved as CSV and Excel files.

## Usage

```text
python -m pip install -r requirements.txt
python profoma_parser.py
```

Set `INPUT_FILE`, `OUTPUT_FILE`, and column names in `profoma_parser.py`. Input requires `proForma` and `protein_accessions` columns.

The published compendium and its outputs are in `data/`. Internet access may be needed to retrieve UniProt and modification data.
