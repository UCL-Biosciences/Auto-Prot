import pandas as pd
import pytest

from src.utils.data_io import load_data, make_outdir


def test_load_data_csv(tmp_path):
    path = tmp_path / "test.csv"
    path.write_text("a,b\n1,2\n3,4")
    df = load_data(str(path))
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["a", "b"]

def test_load_data_tsv(tmp_path):
    path = tmp_path / "test.tsv"
    path.write_text("a\tb\n1\t2\n3\t4")
    df = load_data(str(path))
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["a", "b"]

def test_load_data_invalid_extension(tmp_path):
    path = tmp_path / "test.xlsx"
    path.write_text("junk")
    with pytest.raises(ValueError, match="Unsupported file format"):
        load_data(str(path))

def test_load_data_strips_spaces(tmp_path):
    f = tmp_path / "test.csv"
    f.write_text(" sample_id , treatment \n S1 , Ctl \n S2 , Drug ")
    df = load_data(str(f))
    assert list(df.columns) == ["sample_id", "treatment"]
    assert df.iloc[0, 0] == "S1"
    assert df.iloc[1, 1] == "Drug"

def test_make_outdir_creates_dir(tmp_path):
    out = tmp_path / "myout"
    make_outdir(str(out))
    assert out.exists()
    assert (out / "data").exists()
    assert (out / "plots").exists()

def test_make_outdir_no_subdirs(tmp_path):
    out = tmp_path / "flat"
    make_outdir(str(out), make_subdirs=False)
    assert out.exists()
    assert not (out / "data").exists()
    assert not (out / "plots").exists()

def test_make_outdir_when_dir_already_exists(tmp_path):
    out = tmp_path / "existing"
    out.mkdir()
    make_outdir(str(out))  # Should not error
    assert out.exists()