from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

@dataclass(frozen=True)
class AnalysisConfig:
    background_threshold: float = 0.30
    fdr_threshold: float = 0.10

def _require(df, cols, name):
    missing=set(cols)-set(df.columns)
    if missing: raise ValueError(f"{name} missing columns: {sorted(missing)}")
    if df.empty: raise ValueError(f"{name} is empty")

def validate(metadata, pirna, transposons, imaging, germ):
    _require(metadata, ["sample_id","genotype","stage","batch","replicate","tissue","library_type"], "metadata")
    _require(pirna, ["sample_id","pirna_id","transposon_family","cluster","sequence","length","count"], "pirna_counts")
    _require(transposons, ["sample_id","transposon_family","expression_tpm"], "transposon_expression")
    _require(imaging, ["image_id","sample_id","genotype","batch","replicate","aub_posterior_enrichment","germ_plasm_area_um2","aub_granule_count","background_intensity"], "imaging")
    _require(germ, ["experiment_id","genotype","batch","replicate","total_embryos","mean_germ_cells","count_low_germ_cells"], "germ_cells")
    if metadata["sample_id"].duplicated().any(): raise ValueError("sample_id must be unique")
    samples=set(metadata["sample_id"].astype(str))
    for name,df in [("pirna",pirna),("transposons",transposons),("imaging",imaging)]:
        unknown=set(df["sample_id"].astype(str))-samples
        if unknown: raise ValueError(f"{name} has unknown sample IDs: {sorted(unknown)[:5]}")
    if pirna[["sample_id","pirna_id"]].duplicated().any(): raise ValueError("duplicate sample/pirna rows")
    for col in ["length","count"]: pirna[col]=pd.to_numeric(pirna[col], errors="raise")
    transposons["expression_tpm"]=pd.to_numeric(transposons["expression_tpm"], errors="raise")
    for col in ["aub_posterior_enrichment","germ_plasm_area_um2","aub_granule_count","background_intensity"]:
        imaging[col]=pd.to_numeric(imaging[col], errors="raise")
    for col in ["total_embryos","mean_germ_cells","count_low_germ_cells"]:
        germ[col]=pd.to_numeric(germ[col], errors="raise")
    if (pirna["count"]<0).any(): raise ValueError("piRNA counts must be non-negative")
    if (transposons["expression_tpm"]<0).any(): raise ValueError("transposon TPM must be non-negative")
    bad=(germ["total_embryos"]<=0)|(germ["count_low_germ_cells"]<0)|(germ["count_low_germ_cells"]>germ["total_embryos"])
    if bad.any(): raise ValueError("invalid germ-cell count rows")
    return metadata.copy(), pirna.copy(), transposons.copy(), imaging.copy(), germ.copy()

def pirna_cpm(pirna):
    df=pirna.copy()
    lib=df.groupby("sample_id")["count"].transform("sum")
    df["cpm"]=df["count"]/lib.replace(0,np.nan)*1_000_000
    df["log2_cpm"]=np.log2(df["cpm"].fillna(0)+1)
    df["has_1u"]=df["sequence"].str.upper().str[0].eq("U")
    df["has_10a"]=df["sequence"].str.upper().str.len().ge(10)&df["sequence"].str.upper().str[9].eq("A")
    return df

def pirna_summaries(pirna, metadata):
    cpm=pirna_cpm(pirna)
    sample=(cpm.groupby("sample_id",as_index=False)
        .agg(total_reads=("count","sum"),detected_pirnas=("count",lambda s:int((s>0).sum())),mean_length=("length","mean"),one_u_fraction=("has_1u","mean"),ten_a_fraction=("has_10a","mean")))
    sample=metadata.merge(sample,on="sample_id",how="left")
    x=cpm.merge(metadata[["sample_id","genotype","batch"]],on="sample_id")
    family=(x.groupby(["genotype","batch","transposon_family"],as_index=False)
        .agg(total_count=("count","sum"),mean_cpm=("cpm","mean"),median_log2_cpm=("log2_cpm","median"),detected_pirnas=("pirna_id","nunique")))
    cluster=(x.groupby(["genotype","batch","cluster"],as_index=False)
        .agg(total_count=("count","sum"),mean_cpm=("cpm","mean"),detected_pirnas=("pirna_id","nunique")))
    length=x.groupby(["genotype","length"],as_index=False)["count"].sum()
    length["fraction"]=length["count"]/length.groupby("genotype")["count"].transform("sum").replace(0,np.nan)
    bias=(x.groupby(["genotype","batch","transposon_family"],as_index=False)
        .agg(weighted_1u_fraction=("has_1u","mean"),weighted_10a_fraction=("has_10a","mean")))
    return cpm,sample,family,cluster,length,bias

def pingpong(cpm, metadata):
    x=cpm.merge(metadata[["sample_id","genotype","batch"]],on="sample_id")
    rows=[]
    for keys,g in x.groupby(["genotype","batch","transposon_family"],sort=True):
        w=g["count"].clip(lower=0).to_numpy(float)
        motif=np.average(g["has_1u"] & g["has_10a"], weights=w) if w.sum()>0 else np.nan
        score=float(motif*np.log2(g["count"].sum()+1)) if np.isfinite(motif) else np.nan
        rows.append({"genotype":keys[0],"batch":keys[1],"transposon_family":keys[2],"motif_fraction":motif,"pingpong_score":score})
    family=pd.DataFrame(rows)
    geno=(family.groupby("genotype",as_index=False)
        .agg(mean_pingpong_score=("pingpong_score","mean"),median_pingpong_score=("pingpong_score","median"),families=("transposon_family","nunique")))
    return family, geno

def imaging_summaries(imaging, config):
    df=imaging.copy()
    df["background_qc_fail"]=df["background_intensity"]>config.background_threshold
    df["metric_qc_fail"]=(df["aub_posterior_enrichment"]<=0)|(df["germ_plasm_area_um2"]<0)|(df["aub_granule_count"]<0)
    df["imaging_qc_pass"]=~(df["background_qc_fail"]|df["metric_qc_fail"])
    passed=df[df["imaging_qc_pass"]]
    aub=(passed.groupby(["genotype","batch"],as_index=False)
        .agg(images=("image_id","nunique"),mean_aub_posterior_enrichment=("aub_posterior_enrichment","mean"),sem_aub_posterior_enrichment=("aub_posterior_enrichment","sem"),mean_aub_granule_count=("aub_granule_count","mean")))
    gp=(passed.groupby(["genotype","batch"],as_index=False)
        .agg(images=("image_id","nunique"),mean_germ_plasm_area_um2=("germ_plasm_area_um2","mean"),sem_germ_plasm_area_um2=("germ_plasm_area_um2","sem"),median_background_intensity=("background_intensity","median")))
    return df,aub,gp

def wilson(count,total,z=1.96):
    p=count/total
    den=1+z*z/total
    center=(p+z*z/(2*total))/den
    margin=z*np.sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return max(0,center-margin), min(1,center+margin)

def germ_summaries(germ):
    reps=germ.copy()
    reps["low_germ_cell_fraction"]=reps["count_low_germ_cells"]/reps["total_embryos"]
    intervals=reps.apply(lambda r:wilson(r["count_low_germ_cells"],r["total_embryos"]),axis=1)
    reps["low_fraction_wilson_low"]=[v[0] for v in intervals]
    reps["low_fraction_wilson_high"]=[v[1] for v in intervals]
    summary=(reps.groupby("genotype",as_index=False)
        .agg(experiments=("experiment_id","nunique"),total_embryos=("total_embryos","sum"),mean_germ_cells=("mean_germ_cells","mean"),sem_germ_cells=("mean_germ_cells","sem"),low_germ_cell_fraction=("low_germ_cell_fraction","mean")))
    return reps,summary

def integrate(sample_qc, pingpong_geno, aub, gp, germ_summary, transposons, metadata):
    p=sample_qc.groupby("genotype",as_index=False).agg(mean_total_pirna_reads=("total_reads","mean"),mean_detected_pirnas=("detected_pirnas","mean"),mean_1u_fraction=("one_u_fraction","mean"),mean_10a_fraction=("ten_a_fraction","mean"))
    a=aub.groupby("genotype",as_index=False).agg(mean_aub_posterior_enrichment=("mean_aub_posterior_enrichment","mean"),mean_aub_granule_count=("mean_aub_granule_count","mean"))
    g=gp.groupby("genotype",as_index=False).agg(mean_germ_plasm_area_um2=("mean_germ_plasm_area_um2","mean"))
    t=transposons.merge(metadata[["sample_id","genotype"]],on="sample_id").groupby("genotype",as_index=False).agg(mean_transposon_tpm=("expression_tpm","mean"))
    out=(p.merge(pingpong_geno[["genotype","mean_pingpong_score"]],on="genotype",how="left")
        .merge(a,on="genotype",how="left").merge(g,on="genotype",how="left")
        .merge(germ_summary[["genotype","mean_germ_cells","low_germ_cell_fraction"]],on="genotype",how="left")
        .merge(t,on="genotype",how="left"))
    wt=out[out["genotype"]=="WT"]["mean_germ_cells"]
    mut=out[out["genotype"].str.contains("mutant",case=False,na=False)]["mean_germ_cells"]
    if len(wt) and len(mut) and float(wt.iloc[0]-mut.iloc[0])!=0:
        out["rescue_index"]=(out["mean_germ_cells"]-float(mut.iloc[0]))/float(wt.iloc[0]-mut.iloc[0])
    else:
        out["rescue_index"]=np.nan
    corr_rows=[]
    for label,xcol,ycol in [
        ("pingpong_score_vs_aub_posterior_enrichment","mean_pingpong_score","mean_aub_posterior_enrichment"),
        ("aub_posterior_enrichment_vs_germ_cell_count","mean_aub_posterior_enrichment","mean_germ_cells"),
    ]:
        sub=out[[xcol,ycol]].dropna()
        if len(sub)>=3: rho,pv=spearmanr(sub[xcol],sub[ycol])
        else: rho,pv=np.nan,np.nan
        corr_rows.append({"comparison":label,"n":len(sub),"spearman_rho":rho,"p_value":pv})
    ranked=out.copy()
    ranked["rank_by_germ_cells"]=ranked["mean_germ_cells"].rank(ascending=False,method="dense")
    return out.sort_values("mean_germ_cells",ascending=False), pd.DataFrame(corr_rows), ranked.sort_values("rank_by_germ_cells")
