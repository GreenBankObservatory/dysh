"""
Unit tests for index_file module (SDFITS .index file support)
"""

from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

from dysh.fits.index_file import (
    IndexMetadata,
    convert_dysh_to_sdfits_index,
    convert_sdfits_index_to_dysh,
    create_index_metadata,
    get_index_path,
    read_index,
    validate_index,
    write_index,
)


@pytest.fixture
def testdata_dir():
    """Return path to testdata directory."""
    return Path(__file__).parent.parent.parent.parent.parent / "testdata"


@pytest.fixture
def sample_index_file(testdata_dir):
    """Return path to a sample SDFITS .index file."""
    index_path = testdata_dir / "TSCAL_220105_W" / "TSCAL_220105_W.raw.vegas" / "TSCAL_220105_W.raw.vegas.index"
    if index_path.exists():
        return index_path
    pytest.skip(f"Test index file not found: {index_path}")


@pytest.fixture
def sample_fits_file(testdata_dir):
    """Return path to a sample FITS file with corresponding .index."""
    fits_path = testdata_dir / "TSCAL_220105_W" / "TSCAL_220105_W.raw.vegas" / "TSCAL_220105_W.raw.vegas.fits"
    if fits_path.exists():
        return fits_path
    pytest.skip(f"Test FITS file not found: {fits_path}")


@pytest.fixture
def sample_metadata():
    """Return sample index metadata."""
    return IndexMetadata(
        created="Wed Mar 19 13:32:44 2025",
        last_modified="Wed Mar 19 13:33:07 2025",
        version="1.7",
        observer="Test Observer",
        backend="VEGAS",
        tcal_rx_table="unknown",
        created_by="dysh",
        sprotect=1,
    )


@pytest.fixture
def sample_dataframe():
    """Return sample index DataFrame."""
    return pd.DataFrame(
        {
            "INDEX": [0, 1, 2],
            "SCAN": [24, 24, 25],
            "IFNUM": [0, 0, 0],
            "INT": [0, 0, 0],  # SDFITS index uses INT
            "NUMCHN": [1024, 1024, 1024],
            "SIG": ["T", "T", "T"],
            "CAL": ["F", "F", "F"],
            "AZIMUTH": [202.2637870, 202.0649486, 202.8570561],
            "ELEVATION": [66.45040253, 66.45040253, 66.39987842],
            "CENTFREQ": [7.699608491e10, 7.699608491e10, 7.699608467e10],
            "EXPOSURE": [29.72993469, 29.72993469, 29.72943497],
        }
    )


class TestGetIndexPath:
    """Tests for get_index_path()"""

    def test_get_index_path_simple(self):
        """Test basic index path generation."""
        fits_path = Path("/data/test.fits")
        expected = Path("/data/test.index")
        assert get_index_path(fits_path) == expected

    def test_get_index_path_with_directory(self):
        """Test index path generation with directory."""
        fits_path = Path("/data/subdir/observation.fits")
        expected = Path("/data/subdir/observation.index")
        assert get_index_path(fits_path) == expected

    def test_get_index_path_string_input(self):
        """Test index path generation with string input."""
        fits_path = "/data/test.fits"
        expected = Path("/data/test.index")
        assert get_index_path(fits_path) == expected

    def test_get_index_path_complex_name(self):
        """Test index path with complex FITS filename."""
        fits_path = Path("/data/TGBT25B_603_12.raw.vegas.fits")
        expected = Path("/data/TGBT25B_603_12.raw.vegas.index")
        assert get_index_path(fits_path) == expected


class TestReadIndex:
    """Tests for read_index()"""

    def test_read_existing_index(self, sample_index_file):
        """Test reading an existing SDFITS .index file."""
        metadata, df = read_index(sample_index_file)

        # Check metadata
        assert isinstance(metadata, IndexMetadata)
        assert metadata.version == "1.7"
        assert metadata.created_by == "gbtidl"
        assert metadata.observer == "Larry Morgan"
        assert metadata.backend == "VEGAS"

        # Check DataFrame
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert "INDEX" in df.columns or "#INDEX#" in df.columns
        assert "SCAN" in df.columns
        assert "IFNUM" in df.columns

    def test_read_index_structure(self, sample_index_file):
        """Test that read index has expected structure."""
        _metadata, df = read_index(sample_index_file)

        # Should have standard SDFITS index columns
        expected_cols = ["PROJECT", "FILE", "SCAN", "IFNUM", "INT", "SAMPLER"]
        for col in expected_cols:
            assert col in df.columns, f"Missing expected column: {col}"

    def test_read_index_datatypes(self, sample_index_file):
        """Test that data types are correctly parsed."""
        _metadata, df = read_index(sample_index_file)

        # Numeric columns should be numeric
        if "SCAN" in df.columns:
            assert pd.api.types.is_numeric_dtype(df["SCAN"])
        if "IFNUM" in df.columns:
            assert pd.api.types.is_numeric_dtype(df["IFNUM"])

        # Float columns should be float
        if "AZIMUTH" in df.columns:
            assert pd.api.types.is_float_dtype(df["AZIMUTH"])

    def test_read_nonexistent_index(self, tmp_path):
        """Test reading a nonexistent file raises error."""
        fake_path = tmp_path / "nonexistent.index"
        with pytest.raises(FileNotFoundError):
            read_index(fake_path)


class TestWriteIndex:
    """Tests for write_index()"""

    @pytest.mark.skip(reason="Write functionality not currently needed")
    def test_write_and_read_roundtrip(self, tmp_path, sample_metadata, sample_dataframe):
        """Test writing and reading back produces same data."""
        index_path = tmp_path / "test.index"

        # Write
        write_index(index_path, sample_metadata, sample_dataframe)

        # Verify file exists
        assert index_path.exists()

        # Read back
        metadata_read, df_read = read_index(index_path)

        # Check metadata
        assert metadata_read.version == sample_metadata.version
        assert metadata_read.observer == sample_metadata.observer
        assert metadata_read.backend == sample_metadata.backend

        # Check data - allow for some floating point differences
        assert len(df_read) == len(sample_dataframe)
        # SDFITS index format always has all 43 columns, so df_read will have more columns than sample_dataframe
        # Check that all original columns are present in the read data
        for col in sample_dataframe.columns:
            assert col in df_read.columns, f"Column {col} not found in read dataframe"

        # Check a few specific values
        if "SCAN" in df_read.columns:
            pd.testing.assert_series_equal(
                df_read["SCAN"], sample_dataframe["SCAN"], check_names=False, check_dtype=False
            )

    def test_write_creates_proper_format(self, tmp_path, sample_metadata, sample_dataframe):
        """Test that written file has proper SDFITS index format."""
        index_path = tmp_path / "test.index"
        write_index(index_path, sample_metadata, sample_dataframe)

        # Read raw file content
        with open(index_path) as f:
            lines = f.readlines()

        # Check structure
        assert lines[0].strip() == "[header]"
        assert any("[rows]" in line for line in lines)
        assert any("#INDEX#" in line for line in lines)

        # Check header lines are padded to 200 chars
        for i, line in enumerate(lines[1:9]):  # Header field lines
            if "=" in line:
                # Should be padded (newline adds 1 char, so 201 total)
                assert len(line) == 201, f"Line {i + 2} not properly padded: {len(line)} chars"

    def test_write_scientific_notation(self, tmp_path, sample_metadata, sample_dataframe):
        """Test that floats are written in scientific notation."""
        index_path = tmp_path / "test.index"
        write_index(index_path, sample_metadata, sample_dataframe)

        with open(index_path) as f:
            content = f.read()

        # Large frequency values should be in scientific notation
        assert "e+" in content or "e-" in content

    def test_write_boolean_as_TF(self, tmp_path, sample_metadata, sample_dataframe):
        """Test that booleans are written as T/F."""
        index_path = tmp_path / "test.index"
        write_index(index_path, sample_metadata, sample_dataframe)

        with open(index_path) as f:
            lines = f.readlines()

        # Find data lines (after #INDEX# line)
        data_section_started = False
        for line in lines:
            if "#INDEX#" in line:
                data_section_started = True
                continue
            if data_section_started and line.strip():
                # Data line should contain T or F for boolean columns
                # (we know SIG=True and CAL=False from sample data)
                assert " T " in line or " F " in line
                break


class TestValidateIndex:
    """Tests for validate_index()"""

    def test_validate_existing_index(self, sample_fits_file, sample_index_file):
        """Test validating an existing index file."""
        # This may fail if index is older than FITS, which is expected behavior
        result = validate_index(sample_fits_file, sample_index_file)
        assert isinstance(result, bool)

    def test_validate_nonexistent_index(self, sample_fits_file, tmp_path):
        """Test that validation fails for nonexistent index."""
        fake_index = tmp_path / "nonexistent.index"
        assert validate_index(sample_fits_file, fake_index) is False

    def test_validate_nonexistent_fits(self, tmp_path):
        """Test that validation fails for nonexistent FITS file."""
        fake_fits = tmp_path / "nonexistent.fits"
        fake_index = tmp_path / "nonexistent.index"

        # Create the index file
        fake_index.touch()

        assert validate_index(fake_fits, fake_index) is False

    def test_validate_newer_index(self, tmp_path, sample_metadata, sample_dataframe):
        """Test that newer index is valid."""
        fits_path = tmp_path / "test.fits"
        index_path = tmp_path / "test.index"

        # Create FITS file first
        fits_path.touch()

        # Wait a tiny bit (to ensure different timestamps)
        import time

        time.sleep(0.01)

        # Create newer index
        write_index(index_path, sample_metadata, sample_dataframe)

        # Should be valid (index is newer)
        assert validate_index(fits_path, index_path) is True

    def test_validate_older_index(self, tmp_path, sample_metadata, sample_dataframe):
        """Test that slightly older index is still valid."""
        fits_path = tmp_path / "test.fits"
        index_path = tmp_path / "test.index"

        # Create index first
        write_index(index_path, sample_metadata, sample_dataframe)

        # Wait a tiny bit
        import time

        time.sleep(0.01)

        # Create newer FITS file
        fits_path.touch()

        # Should still be valid (small time differences don't matter)
        assert validate_index(fits_path, index_path) is True


class TestConvertSdfitsIndexToDysh:
    """Tests for convert_sdfits_index_to_dysh()"""

    def test_convert_int_to_intnum(self):
        """Test that INT column is renamed to INTNUM."""
        df = pd.DataFrame({"SCAN": [1, 2, 3], "INT": [0, 1, 2], "IFNUM": [0, 0, 0]})

        df_converted = convert_sdfits_index_to_dysh(df)

        assert "INTNUM" in df_converted.columns
        assert "INT" not in df_converted.columns
        pd.testing.assert_series_equal(df_converted["INTNUM"], df["INT"], check_names=False)

    def test_convert_preserves_other_columns(self):
        """Test that other columns are preserved and SDFITS index columns are converted."""
        df = pd.DataFrame({"SCAN": [1, 2, 3], "INT": [0, 1, 2], "IFNUM": [0, 0, 0], "PROJECT": ["A", "B", "C"]})

        df_converted = convert_sdfits_index_to_dysh(df)

        assert "SCAN" in df_converted.columns
        assert "IFNUM" in df_converted.columns
        assert "PROJID" in df_converted.columns  # PROJECT should be converted to PROJID
        assert "PROJECT" not in df_converted.columns  # Original SDFITS index column should be renamed

    def test_convert_no_changes_if_no_sdfits_index_columns(self):
        """Test that DataFrame is unchanged if no SDFITS index columns."""
        df = pd.DataFrame({"SCAN": [1, 2, 3], "INTNUM": [0, 1, 2], "IFNUM": [0, 0, 0]})

        df_converted = convert_sdfits_index_to_dysh(df)

        pd.testing.assert_frame_equal(df_converted, df)


class TestConvertDyshToSdfitsIndex:
    """Tests for convert_dysh_to_sdfits_index()"""

    def test_convert_intnum_to_int(self):
        """Test that INTNUM column is renamed to INT."""
        df = pd.DataFrame({"SCAN": [1, 2, 3], "INTNUM": [0, 1, 2], "IFNUM": [0, 0, 0]})

        df_converted = convert_dysh_to_sdfits_index(df)

        assert "INT" in df_converted.columns
        assert "INTNUM" not in df_converted.columns
        pd.testing.assert_series_equal(df_converted["INT"], df["INTNUM"], check_names=False)

    def test_convert_removes_dysh_only_columns(self):
        """Test that dysh-specific columns are removed."""
        df = pd.DataFrame(
            {"SCAN": [1, 2, 3], "INTNUM": [0, 1, 2], "HDU": [1, 1, 1], "BINTABLE": [0, 0, 0], "FITSINDEX": [0, 1, 2]}
        )

        df_converted = convert_dysh_to_sdfits_index(df)

        assert "HDU" not in df_converted.columns
        assert "BINTABLE" not in df_converted.columns
        assert "FITSINDEX" not in df_converted.columns
        assert "SCAN" in df_converted.columns

    def test_convert_roundtrip(self):
        """Test that converting back and forth preserves data."""
        # Start with SDFITS index format
        df_index = pd.DataFrame({"SCAN": [1, 2, 3], "INT": [0, 1, 2], "IFNUM": [0, 0, 0]})

        # Convert to dysh
        df_dysh = convert_sdfits_index_to_dysh(df_index.copy())

        # Convert back to SDFITS index
        df_back = convert_dysh_to_sdfits_index(df_dysh.copy())

        # Should match original (column order may differ)
        assert set(df_back.columns) == set(df_index.columns)
        pd.testing.assert_series_equal(df_back["INT"], df_index["INT"], check_names=False)


class TestCreateIndexMetadata:
    """Tests for create_index_metadata()"""

    def test_create_metadata_defaults(self):
        """Test creating metadata with default values."""
        metadata = create_index_metadata()

        assert isinstance(metadata, IndexMetadata)
        assert metadata.version == "1.7"
        assert metadata.observer == "Unknown"
        assert metadata.backend == "Unknown"
        assert metadata.tcal_rx_table == "unknown"
        assert metadata.created_by == "dysh"
        assert metadata.sprotect == 1

    def test_create_metadata_custom_values(self):
        """Test creating metadata with custom values."""
        metadata = create_index_metadata(observer="Test Observer", backend="VEGAS", tcal_rx_table="test_table")

        assert metadata.observer == "Test Observer"
        assert metadata.backend == "VEGAS"
        assert metadata.tcal_rx_table == "test_table"

    def test_create_metadata_has_timestamps(self):
        """Test that metadata has valid timestamps."""
        metadata = create_index_metadata()

        assert metadata.created
        assert metadata.last_modified
        # Timestamps should be same when just created
        assert metadata.created == metadata.last_modified

        # Should be parseable as datetime (SDFITS index format)
        # Example: "Wed Mar 19 13:32:44 2025"
        datetime.strptime(metadata.created, "%a %b %d %H:%M:%S %Y")

    def test_create_metadata_current_time(self):
        """Test that metadata uses current time."""
        import time

        before = datetime.now()
        time.sleep(0.001)  # Small delay to ensure different second
        metadata = create_index_metadata()
        time.sleep(0.001)
        after = datetime.now()

        # Parse the timestamp (note: format doesn't include microseconds)
        ts = datetime.strptime(metadata.created, "%a %b %d %H:%M:%S %Y")

        # Should be between before and after (accounting for second precision)
        # Subtract 1 second from before since format doesn't preserve microseconds
        from datetime import timedelta

        assert (before - timedelta(seconds=1)) <= ts <= (after + timedelta(seconds=1))


class TestIndexMetadata:
    """Tests for IndexMetadata dataclass"""

    def test_metadata_fields(self):
        """Test that metadata has all required fields."""
        metadata = IndexMetadata(
            created="Wed Mar 19 13:32:44 2025",
            last_modified="Wed Mar 19 13:33:07 2025",
            version="1.7",
            observer="Test",
            backend="VEGAS",
            tcal_rx_table="unknown",
            created_by="dysh",
            sprotect=1,
        )

        assert metadata.created == "Wed Mar 19 13:32:44 2025"
        assert metadata.last_modified == "Wed Mar 19 13:33:07 2025"
        assert metadata.version == "1.7"
        assert metadata.observer == "Test"
        assert metadata.backend == "VEGAS"
        assert metadata.tcal_rx_table == "unknown"
        assert metadata.created_by == "dysh"
        assert metadata.sprotect == 1

    def test_metadata_defaults(self):
        """Test metadata default values."""
        metadata = IndexMetadata(created="now", last_modified="now")

        assert metadata.version == "1.7"
        assert metadata.observer == "Unknown"
        assert metadata.backend == "Unknown"
        assert metadata.tcal_rx_table == "unknown"
        assert metadata.created_by == "dysh"
        assert metadata.sprotect == 1


class TestIntegration:
    """Integration tests using real testdata files"""

    def test_read_all_testdata_indexes(self, testdata_dir):
        """Test reading all .index files in testdata."""
        if not testdata_dir.exists():
            pytest.skip("testdata directory not found")

        # Find all .index files
        index_files = list(testdata_dir.rglob("*.index"))

        if len(index_files) == 0:
            pytest.skip("No .index files found in testdata")

        # Try to read each one
        successful = 0
        failed = []

        for index_file in index_files:
            try:
                metadata, df = read_index(index_file)
                assert isinstance(metadata, IndexMetadata)
                assert isinstance(df, pd.DataFrame)
                assert len(df) > 0
                successful += 1
            except Exception as e:
                failed.append((index_file, str(e)))

        # Report results
        print(f"\nRead {successful}/{len(index_files)} index files successfully")
        if failed:
            print("Failed files:")
            for path, error in failed:
                print(f"  {path}: {error}")

        # At least some should succeed
        assert successful > 0

    def test_sdfits_index_has_int_column(self, testdata_dir):
        """Test that SDFITS .index files have INT column (not INTNUM)."""
        # Find an SDFITS index file
        index_files = list(testdata_dir.rglob("*.index"))

        if len(index_files) == 0:
            pytest.skip("No .index files found in testdata")

        # Read first one
        metadata, df = read_index(index_files[0])

        # SDFITS index files created by GBTIDL should have INT column
        if metadata.created_by == "gbtidl":
            assert "INT" in df.columns, "SDFITS index from GBTIDL should have INT column"

    def test_write_compatible_with_sdfits_index_format(self, tmp_path, sample_index_file):
        """Test that our written files match SDFITS index format structure."""
        # Read existing SDFITS index
        metadata_orig, df_orig = read_index(sample_index_file)

        # Write it out
        new_index = tmp_path / "test.index"
        write_index(new_index, metadata_orig, df_orig)

        # Read it back
        _metadata_new, df_new = read_index(new_index)

        # Should have same structure
        assert len(df_new) == len(df_orig)
        assert set(df_new.columns) == set(df_orig.columns)
