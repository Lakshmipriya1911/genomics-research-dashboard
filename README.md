# Genomics Research & NGS Quality Intelligence

Independent portfolio prototype for a Bioinformatics Analyst I interview. Not affiliated with Cincinnati Children's.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

The public research expression file is bundled at `data/GSE232336_RNA-seq_TPM.txt.gz`, so **no upload is required**.

## Data provenance

- Real fixed research data: [NCBI GEO GSE232336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232336), processed RNA-seq TPM table.
- 57,820 gene rows, 10 sample columns, HaCaT and TE-7 cell lines, unstimulated and IL-13-treated conditions.
- The public research data are **not** a real-time clinical feed. The QC and variant modules contain only reproducible **synthetic** demonstration data.
- Research visualizations and descriptive contrasts are exploratory and not clinical findings or validated differential-expression analysis.

## Suggested interview walkthrough

1. Overview: study groups and CCL26 expression in treated versus untreated cell lines.
2. Research explorer: choose CCL26 and other genes, compare TPM, view log2 heatmap, export descriptive contrasts.
3. NGS QC: adjust illustrative thresholds and describe what would be required to validate them.
4. Variant review: explain difference between a UI for review and an actual validated variant calling pipeline.
5. Methods: explain provenance, reproducibility, limitations, and how to make the system production-ready.

## Project structure

- `app.py` Streamlit app
- `requirements.txt` Python dependencies
- `data/GSE232336_RNA-seq_TPM.txt.gz` real public processed RNA-seq matrix

## Interview upgrade (v2)
- Research Explorer: replicate-level scatter, heatmap, and descriptive IL-13 response bar chart from actual GSE232336 TPM data.
- NGS QC: Pass/Review/Fail demonstration labels, configurable thresholds, flag reasons, plots, and CSV audit export.
- NGS pipeline workflow: selectable synthetic run with seven stages, expandable descriptions, expected artifacts, and audit export. It does not connect to Nextflow, SLURM, or patient systems.

### Upgrade an existing local installation
1. Stop Streamlit in Terminal with Ctrl+C.
2. Back up your current project folder (especially any edits to app.py).
3. Extract the upgraded ZIP and open the new `cincy_genomics_dashboard` folder.
4. In Terminal: `cd ` followed by dragging the new folder into Terminal, then press Enter.
5. Run `python3 -m venv .venv`, `source .venv/bin/activate`, `python -m pip install -r requirements.txt`, then `python -m streamlit run app.py`.

Research data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232336


## Advanced RNA-seq analysis

Use **Advanced RNA-seq analysis** in the sidebar for searchable gene-symbol replicate plots, group summary statistics, and full-matrix descriptive IL-13 response rankings. Gene symbols with duplicate source rows are summed for display; stable gene IDs and raw counts are recommended for rigorous statistical analysis. The log2 ratios use mean TPM plus a +1 pseudocount and are **not** differential-expression significance results.


## RNA-seq expression QC & sample relationships
New sidebar page with PCA (NumPy SVD on log2(TPM+1), centered per gene), sample correlation heatmap, within-condition replicate pair correlations with adjustable illustrative flags, and per-sample expression distributions. These are exploratory checks on processed TPM, not raw sequencing QC or validated clinical metrics.
