from pathlib import Path
import numpy as np
import pandas as pd

FAMILIES=["roo","gypsy","copia","HeT-A","TART","I-element"]
CLUSTERS=["42AB","flamenco","20A","cluster_7","cluster_11"]

def _seq(rng, length, motif_strength):
    bases=np.array(list("ACGU"))
    s=rng.choice(bases, size=length).tolist()
    if rng.random() < 0.72*motif_strength: s[0]="U"
    if length >= 10 and rng.random() < 0.62*motif_strength: s[9]="A"
    return "".join(s)

def generate_demo_data(output_dir: Path, seed: int=20260710) -> None:
    rng=np.random.default_rng(seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata=[]; pirna=[]; trans=[]; imaging=[]; germ=[]
    genotypes=[("WT",1.00,1.00,1.00,1.00),("tpp_mutant",0.42,0.45,0.40,0.36),("tpp_rescue",0.74,0.72,0.70,0.68)]
    ids=[f"piR_{i:04d}" for i in range(1,181)]
    info={}
    for i,pid in enumerate(ids):
        info[pid]=(FAMILIES[i%len(FAMILIES)], CLUSTERS[i%len(CLUSTERS)], int(rng.choice(np.arange(23,30))))
    for batch in ["B1","B2","B3"]:
        batch_shift=rng.lognormal(0,0.06)
        for genotype,pfac,motif,afac,gfac in genotypes:
            for rep in range(1,4):
                sid=f"{genotype}_{batch}_R{rep}"
                metadata.append({"sample_id":sid,"genotype":genotype,"stage":"ovary_stage_9","batch":batch,"replicate":f"R{rep}","tissue":"ovary","library_type":"small_rna"})
                for pid in ids:
                    fam,cluster,length=info[pid]
                    mean=90*(1+0.16*FAMILIES.index(fam))*batch_shift*pfac
                    if fam in ["roo","gypsy","copia"]: mean*=1.25
                    count=int(rng.negative_binomial(18,18/(18+max(mean,1))))
                    pirna.append({"sample_id":sid,"pirna_id":pid,"transposon_family":fam,"cluster":cluster,"sequence":_seq(rng,length,motif),"length":length,"count":count})
                for fam in FAMILIES:
                    tpm=rng.lognormal(np.log(8+1.2*FAMILIES.index(fam)),0.18)
                    # Mild family-specific changes only: reduced piRNA does not automatically mean massive de-repression.
                    if genotype=="tpp_mutant" and fam in ["gypsy","I-element"]: tpm*=1.18
                    if genotype=="tpp_rescue" and fam in ["gypsy","I-element"]: tpm*=1.07
                    trans.append({"sample_id":sid,"transposon_family":fam,"expression_tpm":tpm})
                for img in range(1,9):
                    imaging.append({"image_id":f"{sid}_IMG{img:02d}","sample_id":sid,"genotype":genotype,"batch":batch,"replicate":f"R{rep}","aub_posterior_enrichment":rng.normal(2.3*afac,0.18),"germ_plasm_area_um2":max(0,rng.normal(42*afac,4.5)),"aub_granule_count":max(0,rng.poisson(45*afac)),"background_intensity":rng.normal(0.18,0.035)})
                total=int(rng.integers(70,110))
                mean_gc=max(0,rng.normal(28*gfac,2.2))
                low_prob=np.clip(0.06+(1-gfac)*0.70,0.02,0.90)
                germ.append({"experiment_id":f"{sid}_GC","genotype":genotype,"batch":batch,"replicate":f"R{rep}","total_embryos":total,"mean_germ_cells":mean_gc,"count_low_germ_cells":int(rng.binomial(total,low_prob))})
    pd.DataFrame(metadata).to_csv(output_dir/"sample_metadata.csv",index=False)
    pd.DataFrame(pirna).to_csv(output_dir/"pirna_counts.csv",index=False)
    pd.DataFrame(trans).to_csv(output_dir/"transposon_expression.csv",index=False)
    pd.DataFrame(imaging).to_csv(output_dir/"imaging_summary.csv",index=False)
    pd.DataFrame(germ).to_csv(output_dir/"germ_cell_counts.csv",index=False)
