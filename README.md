# Drosophila piRNA germ-plasm QC companion

A reproducible Python workflow for synthetic piRNA abundance QC, simplified ping-pong signature summaries, Aub localization imaging summaries, and germ-cell phenotype integration.

The repository uses synthetic demonstration data. It is not a reproduction of the PNAS paper and does not claim to match the Tomari Lab's internal analysis pipeline.

[Click here to view the Live Interactive Web App Demo]( https://drosophilapirnagermplasmapp-gomrqnjdehruzsiquwvkek.streamlit.app/)

## Local Windows CMD test

```cmd
python -m pip install -r requirements.txt
python -m pip install -e .
python run.py demo --data-dir data --config configs\demo.yaml --output-dir results\demo
python -m pytest -q
```

## Analyze your own exported tables

```cmd
python run.py analyze --metadata data\sample_metadata.csv --pirna data\pirna_counts.csv --transposons data\transposon_expression.csv --imaging data\imaging_summary.csv --germ-cells data\germ_cell_counts.csv --config configs\demo.yaml --output-dir results\analysis
```

## Docker

```cmd
docker build -t drosophila-pirna-germplasm .
docker run --rm -v "%cd%\results:/app/results" drosophila-pirna-germplasm
```

## Required tables

- `sample_metadata.csv`: `sample_id,genotype,stage,batch,replicate,tissue,library_type`
- `pirna_counts.csv`: `sample_id,pirna_id,transposon_family,cluster,sequence,length,count`
- `transposon_expression.csv`: `sample_id,transposon_family,expression_tpm`
- `imaging_summary.csv`: `image_id,sample_id,genotype,batch,replicate,aub_posterior_enrichment,germ_plasm_area_um2,aub_granule_count,background_intensity`
- `germ_cell_counts.csv`: `experiment_id,genotype,batch,replicate,total_embryos,mean_germ_cells,count_low_germ_cells`

## Scientific limitation

The ping-pong score is a simplified synthetic proxy. Real piRNA analysis requires strand-aware read mapping, genome/transposon annotations, overlap profiles, and lab-specific small-RNA preprocessing decisions.


## Streamlit interface

```cmd
python -m pip install -r requirements.txt
python -m pip install -e .
python run.py generate-demo --output-dir data
python -m streamlit run streamlit_app.py
```

Open `http://localhost:8501`.

Docker:

```cmd
docker build -t drosophila-pirna-germplasm-app .
docker run --rm -p 8501:8501 drosophila-pirna-germplasm-app
```
