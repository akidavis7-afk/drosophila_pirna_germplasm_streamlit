import pandas as pd
from pirna_germplasm_qc.core import pirna_cpm

def test_cpm_sums_to_million():
    df=pd.DataFrame({"sample_id":["s1","s1"],"pirna_id":["a","b"],"transposon_family":["roo","gypsy"],"cluster":["42AB","42AB"],"sequence":["UCCCCCCCCA","ACCCCCCCCC"],"length":[25,25],"count":[10,30]})
    out=pirna_cpm(df)
    assert abs(out["cpm"].sum()-1_000_000)<1e-6
    assert out["has_1u"].iloc[0]
