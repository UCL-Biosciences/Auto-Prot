import sys
from pathlib import Path

# Add project root to sys.path
sys.path.append(str(Path(__file__).resolve().parents[1]))


import pytest
import pandas as pd
from utils.data_utils import get_subset

def test_get_subset_happy_path():
    df = pd.DataFrame({'value': [1, 2]}, index=['sample_A', 'sample_B'])
    result = get_subset(df, 'A')
    assert len(result) == 1
    assert 'sample_A' in result.index

def test_get_subset_no_matches():
    df = pd.DataFrame({'value': [1, 2]}, index=['sample_A', 'sample_B'])
    with pytest.raises(ValueError, match="No matches found for subset"):
        get_subset(df, 'C')
