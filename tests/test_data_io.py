import pandas as pd
import pytest

from src.utils.data_io import load_data, make_outdir


def test_load_data_csv(tmp_path):
    """Test loading a CSV file."""
    # Create a temporary CSV file
    path = tmp_path / "test.csv"
    # Write some sample data to the file
    path.write_text("a,b\n1,2\n3,4")
    df = load_data(str(path))
    # Check if the loaded data is a DataFrame and has the correct columns
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["a", "b"]


def test_load_data_tsv(tmp_path):
    """Test loading a TSV file."""
    path = tmp_path / "test.tsv"
    path.write_text("a\tb\n1\t2\n3\t4")
    df = load_data(str(path))
    # Check if the loaded data is a DataFrame and has the correct columns
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["a", "b"]


def test_load_data_invalid_extension(tmp_path):
    """Test loading a file with an unsupported extension."""
    # Create a temporary file with an unsupported extension
    path = tmp_path / "test.xlsx"
    # Write some sample data to the file
    path.write_text("junk")
    # Attempt to load the file and expect a ValueError
    with pytest.raises(ValueError, match="Unsupported file format"):
        load_data(str(path))


def test_load_data_strips_spaces(tmp_path):
    """Test loading a CSV file with spaces in headers."""
    # Create a temporary CSV file with spaces in headers
    f = tmp_path / "test.csv"
    f.write_text(" sample_id , treatment \n S1 , Ctl \n S2 , Drug ")
    # Load the data and check if spaces are stripped from headers
    df = load_data(str(f))
    # Check if the loaded data is a DataFrame and has the correct columns
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["sample_id", "treatment"]
    # Check if the data is loaded correctly
    assert df.iloc[0, 0] == "S1"
    assert df.iloc[1, 1] == "Drug"


def test_make_outdir_creates_dir(tmp_path):
    """Test that make_outdir creates the directory and subdirectories."""
    out = tmp_path / "myout"
    make_outdir(str(out))
    # Check if the directory and subdirectories are created
    assert out.exists()
    assert (out / "data").exists()
    assert (out / "plots").exists()


def test_make_outdir_no_subdirs(tmp_path):
    """Test that make_outdir creates a flat directory without subdirectories."""
    out = tmp_path / "flat"
    make_outdir(str(out), make_subdirs=False)
    # Check if the directory is created without subdirectories
    assert out.exists()
    assert not (out / "data").exists()
    assert not (out / "plots").exists()


def test_make_outdir_when_dir_already_exists(tmp_path):
    """Test that make_outdir does not error when the directory already exists."""
    out = tmp_path / "existing"
    out.mkdir()
    # Ensure the directory exists before calling make_outdir
    make_outdir(str(out))  # Should not error
    assert out.exists()
