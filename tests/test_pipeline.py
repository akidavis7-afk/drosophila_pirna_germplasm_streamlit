from pathlib import Path
import pandas as pd
from pirna_germplasm_qc.synthetic import generate_demo_data
from pirna_germplasm_qc.pipeline import analyze_frames, save_results
from pirna_germplasm_qc.core import AnalysisConfig

def test_demo_pipeline(tmp_path: Path):
    data=tmp_path/"data"; out=tmp_path/"results"; generate_demo_data(data, seed=123)
    results=analyze_frames(pd.read_csv(data/"sample_metadata.csv"),pd.read_csv(data/"pirna_counts.csv"),pd.read_csv(data/"transposon_expression.csv"),pd.read_csv(data/"imaging_summary.csv"),pd.read_csv(data/"germ_cell_counts.csv"),AnalysisConfig())
    save_results(results,out)
    for name in ["pirna_sample_qc.csv","pingpong_score_by_genotype.csv","integrated_genotype_summary.csv","pirna_length_distribution.png","integrated_summary_scatter.png","report.md"]:
        assert (out/name).exists()
    assert set(results["integrated_genotype_summary"]["genotype"])=={"WT","tpp_mutant","tpp_rescue"}
