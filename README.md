# COMAD: COMmunity Assembly Dynamics of microbes

## Installation 

```bash
pip install -e.
```

## Sample run starting with biom file
```bash
comad neufit_biom --biom comad/tests/data/sample_biom --output_filename github_example --output_filepath comad/tests/data/testing_output
```

## Sample run starting with data and taxonomy files
```bash
comad neufit --_data_filename comad/tests/data/sample_data.csv --_taxonomy_filename comad/tests/data/sample_taxonomy.csv --output_filename github_example --output_folder_path comad/tests/data/testing_output/github_example
```
