import pandas as pd
import pytest
from pirna_germplasm_qc.core import validate

def test_missing_columns_error():
    with pytest.raises(ValueError):
        validate(pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame())
