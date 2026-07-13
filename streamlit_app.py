from __future__ import annotations

from pathlib import Path
import sys
ROOT_DIR = Path(__file__).resolve().parent
SRC_DIR = ROOT_DIR / "src"
if SRC_DIR.exists() and str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from io import BytesIO
import zipfile
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from pirna_germplasm_qc.core import AnalysisConfig
from pirna_germplasm_qc.pipeline import analyze_frames
from pirna_germplasm_qc.synthetic import generate_demo_data

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

def ensure_demo():
    if not (DATA_DIR / "sample_metadata.csv").exists():
        generate_demo_data(DATA_DIR)

def load_inputs(use_demo, uploads):
    if use_demo:
        ensure_demo()
        return (
            pd.read_csv(DATA_DIR/"sample_metadata.csv"),
            pd.read_csv(DATA_DIR/"pirna_counts.csv"),
            pd.read_csv(DATA_DIR/"transposon_expression.csv"),
            pd.read_csv(DATA_DIR/"imaging_summary.csv"),
            pd.read_csv(DATA_DIR/"germ_cell_counts.csv"),
        )
    missing=[k for k,v in uploads.items() if v is None]
    if missing:
        raise ValueError("Upload all five CSV files: "+", ".join(missing))
    return tuple(pd.read_csv(uploads[k]) for k in ["metadata","pirna","transposons","imaging","germ_cells"])

def zip_results(results):
    keep=["pirna_sample_qc","pirna_family_summary","pirna_cluster_summary","pirna_length_distribution","pirna_bias_summary","pingpong_family_signature","pingpong_score_by_genotype","imaging_qc","aub_localization_summary","germ_plasm_imaging_summary","germ_cell_phenotype_summary","integrated_genotype_summary","integrated_correlations","genotype_rescue_ranking"]
    buf=BytesIO()
    with zipfile.ZipFile(buf,"w",zipfile.ZIP_DEFLATED) as z:
        for k in keep:
            z.writestr(f"{k}.csv",results[k].to_csv(index=False))
    buf.seek(0)
    return buf.getvalue()

def fig_length(df):
    fig,ax=plt.subplots(figsize=(8,5))
    for genotype,g in df.groupby("genotype"):
        g=g.sort_values("length"); ax.plot(g["length"],g["fraction"],marker="o",label=genotype)
    ax.set_xlabel("piRNA length (nt)"); ax.set_ylabel("Fraction"); ax.set_title("piRNA length distribution")
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout(); return fig

def fig_bar(df, x, y, title, ylabel):
    fig,ax=plt.subplots(figsize=(7,5))
    ax.bar(df[x],df[y]); ax.set_title(title); ax.set_ylabel(ylabel); fig.tight_layout(); return fig

def fig_integrated(df):
    fig,ax=plt.subplots(figsize=(8,5))
    sizes=35+8*df["mean_germ_cells"].fillna(0)
    ax.scatter(df["mean_pingpong_score"],df["mean_aub_posterior_enrichment"],s=sizes,alpha=.75)
    for _,r in df.iterrows():
        ax.annotate(r["genotype"],(r["mean_pingpong_score"],r["mean_aub_posterior_enrichment"]),xytext=(5,5),textcoords="offset points")
    ax.set_xlabel("Mean synthetic ping-pong score"); ax.set_ylabel("Mean Aub posterior enrichment")
    ax.set_title("Integrated piRNA → Aub → germ-cell summary"); ax.grid(alpha=.25); fig.tight_layout(); return fig

st.set_page_config(page_title="Drosophila piRNA germ-plasm companion", page_icon="🧬", layout="wide")
st.title("🧬 Drosophila piRNA / germ-plasm QC companion")
st.caption("Synthetic-data demonstration for piRNA abundance, ping-pong proxy, Aub localization, and germ-cell phenotype integration.")
st.info("This app uses synthetic demonstration data. It is not a reproduction of the PNAS paper or the Tomari Lab's internal workflow.")

with st.sidebar:
    st.header("Settings")
    background_threshold=st.slider("Imaging background threshold",0.05,1.0,0.30,0.01)
    fdr_threshold=st.slider("FDR review threshold",0.01,0.50,0.10,0.01)

tabs=st.tabs(["1. Data","2. Run analysis","3. piRNA QC","4. Aub / germ plasm","5. Germ-cell phenotype","6. Integrated summary","About"])

with tabs[0]:
    source=st.radio("Data source",["Bundled synthetic demonstration","Upload five CSV files"],horizontal=True)
    use_demo=source=="Bundled synthetic demonstration"
    uploads={"metadata":None,"pirna":None,"transposons":None,"imaging":None,"germ_cells":None}
    if use_demo:
        ensure_demo(); st.success("Synthetic demonstration data are ready.")
        st.dataframe(pd.read_csv(DATA_DIR/"sample_metadata.csv"), width="stretch", hide_index=True)
    else:
        c1,c2=st.columns(2)
        uploads["metadata"]=c1.file_uploader("sample_metadata.csv",type="csv")
        uploads["pirna"]=c2.file_uploader("pirna_counts.csv",type="csv")
        uploads["transposons"]=c1.file_uploader("transposon_expression.csv",type="csv")
        uploads["imaging"]=c2.file_uploader("imaging_summary.csv",type="csv")
        uploads["germ_cells"]=c1.file_uploader("germ_cell_counts.csv",type="csv")

with tabs[1]:
    if st.button("Run piRNA / germ-plasm analysis",type="primary",width="stretch"):
        try:
            frames=load_inputs(use_demo,uploads)
            st.session_state["results"]=analyze_frames(*frames, AnalysisConfig(background_threshold, fdr_threshold))
            st.success("Analysis completed.")
        except Exception as exc:
            st.error(str(exc))
    if "results" in st.session_state:
        r=st.session_state["results"]
        a,b,c,d=st.columns(4)
        a.metric("Samples", r["metadata"]["sample_id"].nunique())
        b.metric("piRNA rows", len(r["pirna_cpm"]))
        c.metric("Images passing QC", int(r["imaging_qc"]["imaging_qc_pass"].sum()))
        d.metric("Phenotype experiments", r["germ_cell_replicates"]["experiment_id"].nunique())

with tabs[2]:
    if "results" not in st.session_state: st.warning("Run the analysis first.")
    else:
        r=st.session_state["results"]
        fig=fig_length(r["pirna_length_distribution"]); st.pyplot(fig,width="stretch"); plt.close(fig)
        st.dataframe(r["pirna_sample_qc"],width="stretch",hide_index=True)
        st.dataframe(r["pingpong_score_by_genotype"],width="stretch",hide_index=True)

with tabs[3]:
    if "results" not in st.session_state: st.warning("Run the analysis first.")
    else:
        r=st.session_state["results"]
        aub=r["aub_localization_summary"].groupby("genotype",as_index=False)["mean_aub_posterior_enrichment"].mean()
        fig=fig_bar(aub,"genotype","mean_aub_posterior_enrichment","Aub posterior enrichment","Aub posterior enrichment")
        st.pyplot(fig,width="stretch"); plt.close(fig)
        st.dataframe(r["aub_localization_summary"],width="stretch",hide_index=True)
        st.dataframe(r["germ_plasm_imaging_summary"],width="stretch",hide_index=True)

with tabs[4]:
    if "results" not in st.session_state: st.warning("Run the analysis first.")
    else:
        r=st.session_state["results"]
        fig=fig_bar(r["germ_cell_phenotype_summary"],"genotype","mean_germ_cells","Germ-cell formation by genotype","Mean germ-cell count")
        st.pyplot(fig,width="stretch"); plt.close(fig)
        st.dataframe(r["germ_cell_phenotype_summary"],width="stretch",hide_index=True)

with tabs[5]:
    if "results" not in st.session_state: st.warning("Run the analysis first.")
    else:
        r=st.session_state["results"]
        fig=fig_integrated(r["integrated_genotype_summary"]); st.pyplot(fig,width="stretch"); plt.close(fig)
        st.dataframe(r["integrated_genotype_summary"],width="stretch",hide_index=True)
        st.dataframe(r["integrated_correlations"],width="stretch",hide_index=True)
        st.download_button("Download result tables as ZIP",data=zip_results(r),file_name="drosophila_pirna_germplasm_results.zip",mime="application/zip",width="stretch")

with tabs[6]:
    st.markdown("""
### Scope
This app integrates processed, exported tables. It does not trim adapters, map small-RNA reads, calculate real strand-aware 10-nt overlap profiles, or segment microscopy images.

### Important caveat
The bundled ping-pong metric is a synthetic demonstration proxy. Real piRNA analysis requires genome/transposon annotation, strand-aware mapping, and laboratory-specific small-RNA preprocessing.
""")
