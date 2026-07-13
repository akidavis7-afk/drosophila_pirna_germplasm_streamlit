from pathlib import Path
import pandas as pd
import yaml

from .core import AnalysisConfig, validate, pirna_summaries, pingpong, imaging_summaries, germ_summaries, integrate
from .plots import plot_length, plot_family_heatmap, plot_bar, plot_integrated

def load_config(path: Path) -> AnalysisConfig:
    raw=yaml.safe_load(path.read_text(encoding="utf-8")) if path.exists() else {}
    return AnalysisConfig(**(raw or {}))

def analyze_frames(metadata, pirna, transposons, imaging, germ_cells, config: AnalysisConfig):
    metadata,pirna,transposons,imaging,germ_cells=validate(metadata,pirna,transposons,imaging,germ_cells)
    cpm,sample_qc,family,cluster,length,bias=pirna_summaries(pirna,metadata)
    family_pp,geno_pp=pingpong(cpm,metadata)
    imaging_qc,aub,gp=imaging_summaries(imaging,config)
    germ_reps,germ_summary=germ_summaries(germ_cells)
    integrated,corr,ranked=integrate(sample_qc,geno_pp,aub,gp,germ_summary,transposons,metadata)
    return {
        "metadata":metadata,"pirna_cpm":cpm,"pirna_sample_qc":sample_qc,"pirna_family_summary":family,
        "pirna_cluster_summary":cluster,"pirna_length_distribution":length,"pirna_bias_summary":bias,
        "pingpong_family_signature":family_pp,"pingpong_score_by_genotype":geno_pp,
        "imaging_qc":imaging_qc,"aub_localization_summary":aub,"germ_plasm_imaging_summary":gp,
        "germ_cell_replicates":germ_reps,"germ_cell_phenotype_summary":germ_summary,
        "integrated_genotype_summary":integrated,"integrated_correlations":corr,"genotype_rescue_ranking":ranked,
    }

def save_results(results, output_dir: Path):
    output_dir.mkdir(parents=True,exist_ok=True)
    for key,df in results.items():
        if hasattr(df,"to_csv") and key!="metadata":
            df.to_csv(output_dir/f"{key}.csv",index=False)
    plot_length(results["pirna_length_distribution"], output_dir/"pirna_length_distribution.png")
    plot_family_heatmap(results["pirna_family_summary"], output_dir/"pirna_family_heatmap.png")
    plot_bar(results["pingpong_score_by_genotype"],"genotype","mean_pingpong_score","Ping-pong signature by genotype","Mean synthetic ping-pong score",output_dir/"pingpong_signature.png")
    aub_by=results["aub_localization_summary"].groupby("genotype",as_index=False)["mean_aub_posterior_enrichment"].mean()
    plot_bar(aub_by,"genotype","mean_aub_posterior_enrichment","Aub posterior enrichment","Aub posterior enrichment",output_dir/"aub_localization_by_genotype.png")
    plot_bar(results["germ_cell_phenotype_summary"],"genotype","mean_germ_cells","Germ-cell formation by genotype","Mean germ-cell count",output_dir/"germ_cell_counts_by_genotype.png")
    plot_integrated(results["integrated_genotype_summary"], output_dir/"integrated_summary_scatter.png")
    report=[
        "# Drosophila piRNA / germ-plasm QC companion report","",
        "The bundled data are synthetic. This is a topic-aligned engineering demonstration, not a reproduction of the PNAS paper or the Tomari Lab's internal pipeline.","",
        "## Integrated genotype summary","",results["integrated_genotype_summary"].to_markdown(index=False),"",
        "## Correlations","",results["integrated_correlations"].to_markdown(index=False),"",
        "## Limitation","",
        "The ping-pong score is a simplified synthetic proxy. Real piRNA analysis requires strand-aware read mapping, genome/transposon annotations, overlap profiles, and laboratory-specific preprocessing."
    ]
    (output_dir/"report.md").write_text("\n".join(report),encoding="utf-8")

def run_from_paths(metadata_path: Path,pirna_path: Path,transposon_path: Path,imaging_path: Path,germ_cells_path: Path,config_path: Path,output_dir: Path):
    results=analyze_frames(
        pd.read_csv(metadata_path),pd.read_csv(pirna_path),pd.read_csv(transposon_path),pd.read_csv(imaging_path),pd.read_csv(germ_cells_path),load_config(config_path)
    )
    save_results(results,output_dir)
