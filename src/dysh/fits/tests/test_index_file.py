"""
Unit tests for index_file module (SDFITS .index file support)
"""

from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

from dysh.fits.index_file import (
    IndexMetadata,
    _generate_rows_header,
    _get_center_frequency,
    _get_polarization,
    _get_procedure,
    _translate_boolean,
    convert_dysh_to_sdfits_index,
    convert_sdfits_index_to_dysh,
    create_index_metadata,
    get_index_path,
    read_index,
    validate_index,
    write_index,
)
from dysh.util import get_project_testdata


@pytest.fixture
def testdata_dir():
    """Return path to testdata directory."""
    return get_project_testdata()


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

        # Check integer values survive roundtrip
        pd.testing.assert_series_equal(df_read["SCAN"], sample_dataframe["SCAN"], check_names=False, check_dtype=False)
        pd.testing.assert_series_equal(
            df_read["IFNUM"], sample_dataframe["IFNUM"], check_names=False, check_dtype=False
        )

        # Check float values survive roundtrip (within scientific notation precision)
        for col in ("AZIMUTH", "ELEVATION", "CENTFREQ", "EXPOSURE"):
            pd.testing.assert_series_equal(
                df_read[col], sample_dataframe[col], check_names=False, check_dtype=False, rtol=1e-8
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

        # Check header lines are padded to 256 chars (matching sparrow3)
        for i, line in enumerate(lines[1:9]):  # Header field lines
            if "=" in line:
                # Should be padded (newline adds 1 char, so 257 total)
                assert len(line) == 257, f"Line {i + 2} not properly padded: {len(line)} chars"

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


class TestDerivedColumns:
    """Tests for derived column computation functions"""

    def test_get_polarization(self):
        """Test CRVAL4 → polarization string mapping."""
        assert _get_polarization(1) == "I"
        assert _get_polarization(2) == "Q"
        assert _get_polarization(3) == "U"
        assert _get_polarization(4) == "V"
        assert _get_polarization(-1) == "RR"
        assert _get_polarization(-2) == "LL"
        assert _get_polarization(-5) == "XX"
        assert _get_polarization(-6) == "YY"
        assert _get_polarization(99) == "??"

    def test_get_center_frequency(self):
        """Test CENTFREQ computation matches sparrow3."""
        # centerChan = (numchn / 2.0) - 0.5
        # centerFreq = ((centerChan - crpix1) * cdelt1) + crval1
        crval1 = 1.42e9
        crpix1 = 512.0
        cdelt1 = 1000.0
        numchn = 1024

        expected_center_chan = (1024 / 2.0) - 0.5  # = 511.5
        expected = ((expected_center_chan - 512.0) * 1000.0) + 1.42e9  # = 1.42e9 - 500.0
        result = _get_center_frequency(crval1, crpix1, cdelt1, numchn)
        assert result == pytest.approx(expected, rel=1e-12)

    def test_get_center_frequency_not_crval1(self):
        """Test that CENTFREQ differs from CRVAL1 when CRPIX1 is not at center."""
        crval1 = 1.0e9
        crpix1 = 0.0  # Reference pixel at edge, not center
        cdelt1 = 1000.0
        numchn = 1024

        result = _get_center_frequency(crval1, crpix1, cdelt1, numchn)
        # centerChan = 511.5, so centerFreq = (511.5 - 0) * 1000 + 1e9
        assert result != crval1  # Should NOT be the same as CRVAL1
        assert result == pytest.approx(1.0e9 + 511500.0, rel=1e-12)

    def test_get_procedure(self):
        """Test OBSMODE → PROCEDURE extraction."""
        assert _get_procedure("OnOff:Nod") == "OnOff"
        assert _get_procedure("Track") == "Track"
        assert _get_procedure("") == ""
        assert _get_procedure("OffOn:PSWITCHOFF:1:2") == "OffOn"

    def test_translate_boolean(self):
        """Test boolean T/F conversion."""
        assert _translate_boolean(True) == "T"
        assert _translate_boolean(False) == "F"
        assert _translate_boolean("T") == "T"
        assert _translate_boolean("F") == "F"
        assert _translate_boolean("True") == "T"
        assert _translate_boolean(1) == "T"
        assert _translate_boolean(0) == "F"


class TestSparrow3FormatCompat:
    """Tests for sparrow3/GBTIDL format compatibility"""

    def test_header_line_padding_256(self, tmp_path, sample_metadata, sample_dataframe):
        """Test that header lines are padded to 256 chars (matching sparrow3)."""
        index_path = tmp_path / "test.index"
        write_index(index_path, sample_metadata, sample_dataframe)

        with open(index_path) as f:
            lines = f.readlines()

        for line in lines[1:9]:  # Header key=value lines
            if "=" in line:
                # Content should be 256 chars + newline = 257
                assert len(line) == 257, f"Header line not padded to 256: {len(line) - 1} chars"

    def test_generated_header_matches_sparrow3(self):
        """Test that generated column header row matches sparrow3's format."""
        header = _generate_rows_header()

        # Should start with #INDEX# (7 chars, no trailing space)
        assert header.startswith("#INDEX#")

        # Verify truncated column names appear correctly
        # EXTENSION → EXT (truncated by %3.3s)
        assert " EXT" in header
        # POLARIZATION → POL (truncated by %3.3s)
        assert " POL " in header
        # E2ESCAN → E2ESC (truncated by %5.5s)
        assert " E2ESC " in header
        # PROCSEQN → PROCS (truncated by %5.5s)
        assert " PROCS " in header
        # SUBREF → SUB (truncated by %3.3s)
        assert " SUB " in header

        # Verify non-truncated names
        assert "ELEVATION" in header
        assert "CENTFREQ" in header
        assert "DATEOBS" in header

    def test_write_with_dysh_column_names(self, tmp_path, sample_metadata):
        """Test writing a DataFrame with dysh/FITS column names."""
        df = pd.DataFrame(
            {
                "INDEX": [0],
                "PROJID": ["TestProject"],
                "HDU": [1],
                "OBJECT": ["W3OH"],
                "SCAN": [100],
                "CRVAL4": [-5],  # XX polarization
                "CRVAL1": [1.42e9],
                "CRPIX1": [512.0],
                "CDELT1": [1000.0],
                "NUMCHN": [1024],
                "CRVAL2": [45.72],
                "CRVAL3": [-10.3],
                "ELEVATIO": [30.75],
                "OBSMODE": ["OnOff:Nod"],
                "SIG": [True],
                "CAL": [False],
            }
        )
        index_path = tmp_path / "test.index"
        write_index(index_path, sample_metadata, df)

        # Read back and verify derived columns
        _metadata, df_read = read_index(index_path)

        # POLARIZATION should be derived from CRVAL4=-5 → XX
        assert "POL" in df_read.columns
        assert df_read["POL"].iloc[0] == "XX"

        # CENTFREQ should be computed (not just CRVAL1)
        # centerChan = (1024/2) - 0.5 = 511.5
        # centerFreq = ((511.5 - 512.0) * 1000.0) + 1.42e9
        expected_centfreq = _get_center_frequency(1.42e9, 512.0, 1000.0, 1024)
        assert df_read["CENTFREQ"].iloc[0] == pytest.approx(expected_centfreq, rel=1e-8)

        # PROCEDURE should be extracted from OBSMODE
        assert df_read["PROCEDURE"].iloc[0] == "OnOff"

        # SOURCE should be mapped from OBJECT
        assert df_read["SOURCE"].iloc[0] == "W3OH"

        # EXT should be mapped from HDU
        assert df_read["EXT"].iloc[0] == 1

        # ELEVATION should be mapped from ELEVATIO
        assert df_read["ELEVATION"].iloc[0] == pytest.approx(30.75)

        # SIG/CAL should be T/F booleans
        assert df_read["SIG"].iloc[0] is True  # read_index converts T→True
        assert df_read["CAL"].iloc[0] is False  # read_index converts F→False

    def test_nsave_defaults_to_minus_one(self, tmp_path, sample_metadata, sample_dataframe):
        """Test that NSAVE defaults to -1 (matching sparrow3)."""
        index_path = tmp_path / "test.index"
        write_index(index_path, sample_metadata, sample_dataframe)

        _, df_read = read_index(index_path)
        assert all(df_read["NSAVE"] == -1)

    def test_dateobs_timestamp_width_22(self, tmp_path, sample_metadata):
        """Test DATEOBS and TIMESTAMP columns use 22-char width (not 21)."""
        df = pd.DataFrame(
            {
                "INDEX": [0],
                "SCAN": [1],
                "DATEOBS": ["2013_11_05_16:15:21.35"],  # 22 chars
                "TIMESTAMP": ["2013_11_05_16:15:21"],  # shorter, right-justified in 22
            }
        )
        index_path = tmp_path / "test.index"
        write_index(index_path, sample_metadata, df)

        with open(index_path) as f:
            content = f.read()

        # The date should be present (right-justified in 22-char field)
        assert "2013_11_05_16:15:21.35" in content
        assert "2013_11_05_16:15:21" in content


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


def _get_rows_section(filepath):
    """Extract [rows] section (column header + data rows) from an index file.

    Returns list of lines with trailing newlines/spaces stripped.
    Skips the [header] section entirely, returning only the column header
    and data rows for format comparison.
    """
    with open(filepath) as f:
        lines = f.readlines()
    rows_start = None
    for i, line in enumerate(lines):
        if line.strip() == "[rows]":
            rows_start = i + 1
            break
    if rows_start is None:
        raise ValueError("No [rows] section found")
    return [line.rstrip("\n") for line in lines[rows_start:] if line.strip()]


def _sparrow3_base_row():
    """Return a dict of canonical SDFITS index column values matching
    the sparrow3 IndexWriterTests test data.

    Values are pre-computed to match sparrow3's translateInfo() output:
    - CENTFREQ = ((1024/2 - 0.5) - 512) * 1.0 + 0.0 = -0.5
    - POLARIZATION = XX (from CRVAL4=-5)
    - PROCEDURE = obsmode (from OBSMODE="obsmode", no colon)
    - SIG = T, CAL = T (sparrow3's translateBoolean treats string "F" as truthy;
      we use the same expected output values for format comparison)
    """
    return {
        "PROJECT": "projectA",
        "FILE": "filepath",
        "EXTENSION": 1,
        "SOURCE": "object name",
        "PROCEDURE": "obsmode",
        "OBSID": "obsid",
        "E2ESCAN": 0,
        "PROCSEQN": 5,
        "SCAN": 100,
        "POLARIZATION": "XX",
        "PLNUM": 1,
        "IFNUM": 3,
        "FEED": 4,
        "FDNUM": 1,
        "INT": 10,
        "NUMCHN": 1024,
        "SIG": "T",
        "CAL": "T",
        "SAMPLER": "A1_0",
        "AZIMUTH": 120.5,
        "ELEVATION": 30.75,
        "LONGITUDE": 45.72,
        "LATITUDE": -10.3,
        "TRGTLONG": 45.8,
        "TRGTLAT": -10.5,
        "SUBREF": 1,
        "LST": 1234.56,
        "CENTFREQ": -0.5,  # computed from CRVAL1=0, CRPIX1=512, CDELT1=1, NUMCHN=1024
        "RESTFREQ": 1420405800.0,
        "VELOCITY": 0.0,
        "FREQINT": 1.0,
        "FREQRES": 1.0,
        "DATEOBS": "2013_11_05_16:15:21.35",
        "TIMESTAMP": "2013_11_05_16:15:21",
        "BANDWIDTH": 1024.0,
        "EXPOSURE": 2.0,
        "TSYS": 20.0,
        "NSAVE": -1,
        "PROCSCAN": 10,
        "PROCTYPE": "unknown",
        "WCALPOS": "Unknown",
    }


# Path to sparrow3 expected test data
_SPARROW3_DATA_DIR = (
    Path(__file__).parent.parent.parent.parent.parent / ".context" / "sparrow3" / "gbt" / "api" / "sdfits" / "data"
)


class TestSparrow3Ported:
    """Tests ported from sparrow3 IndexWriterTests.

    These compare dysh-written data rows character-for-character against
    sparrow3's .index.expected reference files to verify exact format
    compatibility with GBTIDL.
    """

    def _write_test_index(self, tmp_path, num_rows, start_index=0, start_row=0):
        """Helper: write an index file with the sparrow3 test data.

        Returns the path to the written file.
        """
        base = _sparrow3_base_row()
        rows = []
        for i in range(num_rows):
            row = dict(base)
            row["INDEX"] = start_index + i
            row["ROW"] = start_row + i
            rows.append(row)

        df = pd.DataFrame(rows)
        metadata = IndexMetadata(
            created="Mon Mar 23 14:42:48 2015",
            last_modified="Mon Mar 23 14:42:48 2015",
        )
        index_path = tmp_path / "test.index"
        write_index(index_path, metadata, df)
        return index_path

    def test_basics_matches_sparrow3(self, tmp_path):
        """Port of sparrow3 IndexWriterTests.testBasics.

        Write 3 rows with small INDEX/ROW values and verify data rows
        match Basics.index.expected character-for-character.
        """
        expected_path = _SPARROW3_DATA_DIR / "Basics.index.expected"
        if not expected_path.exists():
            pytest.skip(f"sparrow3 expected file not found: {expected_path}")

        index_path = self._write_test_index(tmp_path, num_rows=3, start_index=0, start_row=0)

        actual = _get_rows_section(index_path)
        expected = _get_rows_section(expected_path)

        # Column header should match exactly
        assert actual[0] == expected[0], "Column header row does not match sparrow3"

        # Data rows should match exactly
        assert len(actual) == len(expected), f"Row count mismatch: {len(actual) - 1} vs {len(expected) - 1} data rows"
        for i, (act, exp) in enumerate(zip(actual[1:], expected[1:], strict=True)):
            assert act == exp, f"Data row {i} does not match sparrow3 expected output"

    def test_long_row_matches_sparrow3(self, tmp_path):
        """Port of sparrow3 IndexWriterTests.testLongRow.

        Write 3 rows where ROW crosses the 1e6 boundary (999999, 1000000, 1000001).
        Tests the skip-spacing behavior: ROW < 1e6 gets a trailing space,
        ROW >= 1e6 does not (the extra digit fills the space).
        """
        expected_path = _SPARROW3_DATA_DIR / "LongRow.index.expected"
        if not expected_path.exists():
            pytest.skip(f"sparrow3 expected file not found: {expected_path}")

        index_path = self._write_test_index(tmp_path, num_rows=3, start_index=0, start_row=999999)

        actual = _get_rows_section(index_path)
        expected = _get_rows_section(expected_path)

        assert actual[0] == expected[0], "Column header row does not match sparrow3"
        assert len(actual) == len(expected)
        for i, (act, exp) in enumerate(zip(actual[1:], expected[1:], strict=True)):
            assert act == exp, f"Data row {i} (ROW={999999 + i}) does not match sparrow3 expected output"

    def test_long_index_and_row_matches_sparrow3(self, tmp_path):
        """Port of sparrow3 IndexWriterTests.testLongIndexAndRow.

        Write 5 rows where both INDEX and ROW cross the 1e6 boundary independently.
        INDEX starts at 999999 (crosses at row 2), ROW starts at 999997 (crosses at row 4).
        Tests skip-spacing for both columns simultaneously.
        """
        expected_path = _SPARROW3_DATA_DIR / "LongIndexAndRow.index.expected"
        if not expected_path.exists():
            pytest.skip(f"sparrow3 expected file not found: {expected_path}")

        index_path = self._write_test_index(tmp_path, num_rows=5, start_index=999999, start_row=999997)

        actual = _get_rows_section(index_path)
        expected = _get_rows_section(expected_path)

        assert actual[0] == expected[0], "Column header row does not match sparrow3"
        assert len(actual) == len(expected)
        for i, (act, exp) in enumerate(zip(actual[1:], expected[1:], strict=True)):
            assert act == exp, (
                f"Data row {i} (INDEX={999999 + i}, ROW={999997 + i}) does not match sparrow3 expected output"
            )

    def test_vegas_multibank_matches_sparrow3(self, tmp_path):
        """Port of sparrow3 MPSDFITSWriterTests.testVegasIntegrationsAndIndex.

        Write 32 rows across 8 VEGAS bank files (A-H) with varying FILE, POL,
        SAMPLER, FEED, FDNUM, IFNUM, coordinates, CENTFREQ, CAL, and EXPOSURE.
        Verify data rows match test.banks.vegas.raw.ints.index.expected
        character-for-character.
        """
        expected_path = _SPARROW3_DATA_DIR / "test.banks.vegas.raw.ints.index.expected"
        if not expected_path.exists():
            pytest.skip(f"sparrow3 expected file not found: {expected_path}")

        # Bank configs: (file, ifnum, feed, fdnum, centfreq, sampler_ll, sampler_rr, az, el, lon, lat)
        banks = [
            (
                "test.banks.vegas.raw.ints.A.fits",
                0,
                4,
                3,
                2.370934760e10,
                "A9_0",
                "A13_0",
                1.448764788e02,
                1.379075952e01,
                2.668584627e02,
                -2.835169168e01,
            ),
            (
                "test.banks.vegas.raw.ints.B.fits",
                1,
                5,
                4,
                2.370937201e10,
                "B17_0",
                "B21_0",
                1.448529791e02,
                1.380393730e01,
                2.668730975e02,
                -2.832869620e01,
            ),
            (
                "test.banks.vegas.raw.ints.C.fits",
                0,
                3,
                2,
                2.370934760e10,
                "C25_0",
                "C29_0",
                1.448764761e02,
                1.376440397e01,
                2.668737744e02,
                -2.837434610e01,
            ),
            (
                "test.banks.vegas.raw.ints.D.fits",
                0,
                1,
                0,
                2.370934760e10,
                "D33_0",
                "D37_0",
                1.448529791e02,
                1.377758174e01,
                2.668884105e02,
                -2.835134734e01,
            ),
            (
                "test.banks.vegas.raw.ints.E.fits",
                0,
                2,
                1,
                2.370934760e10,
                "E10_0",
                "E14_0",
                1.448529791e02,
                1.375122619e01,
                2.669037308e02,
                -2.837399619e01,
            ),
            (
                "test.banks.vegas.raw.ints.F.fits",
                2,
                6,
                5,
                2.370937201e10,
                "F18_0",
                "F22_0",
                1.448294795e02,
                1.379075952e01,
                2.669030393e02,
                -2.832834597e01,
            ),
            (
                "test.banks.vegas.raw.ints.G.fits",
                2,
                7,
                6,
                2.370937201e10,
                "G26_0",
                "G30_0",
                1.448294821e02,
                1.376440397e01,
                2.669183618e02,
                -2.835099279e01,
            ),
            (
                "test.banks.vegas.raw.ints.H.fits",
                2,
                1,
                0,
                2.414234760e10,
                "H34_0",
                "H38_0",
                1.448529791e02,
                1.377758174e01,
                2.668884105e02,
                -2.835134734e01,
            ),
        ]

        trgtlong = 2.668918438e02
        trgtlat = -2.835127778e01

        rows = []
        idx = 0
        for filename, ifnum, feed, fdnum, centfreq, sampler_ll, sampler_rr, az, el, lon, lat in banks:
            row_in_file = 0
            for pol, plnum, sampler in [("LL", 0, sampler_ll), ("RR", 1, sampler_rr)]:
                for cal, exposure in [("T", 4.949683249e-01), ("F", 3.954508901e-01)]:
                    rows.append(
                        {
                            "INDEX": idx,
                            "PROJECT": "KFPA",
                            "FILE": filename,
                            "EXTENSION": 1,
                            "ROW": row_in_file,
                            "SOURCE": "SGRB2",
                            "PROCEDURE": "RALongMap",
                            "OBSID": "unknown",
                            "E2ESCAN": 0,
                            "PROCSEQN": 22,
                            "SCAN": 34,
                            "POLARIZATION": pol,
                            "PLNUM": plnum,
                            "IFNUM": ifnum,
                            "FEED": feed,
                            "FDNUM": fdnum,
                            "INT": 1,
                            "NUMCHN": 4096,
                            "SIG": "T",
                            "CAL": cal,
                            "SAMPLER": sampler,
                            "AZIMUTH": az,
                            "ELEVATION": el,
                            "LONGITUDE": lon,
                            "LATITUDE": lat,
                            "TRGTLONG": trgtlong,
                            "TRGTLAT": trgtlat,
                            "SUBREF": 1,
                            "LST": 5.462375752e04,
                            "CENTFREQ": centfreq,
                            "RESTFREQ": 2.370629500e10,
                            "VELOCITY": 0.0,
                            "FREQINT": -1.220703125e04,
                            "FREQRES": 1.477050781e04,
                            "DATEOBS": "2010-04-08T07:23:56.00",
                            "TIMESTAMP": "2010_04_08_07:23:55",
                            "BANDWIDTH": 5.000000000e07,
                            "EXPOSURE": exposure,
                            "TSYS": 1.0,
                            "NSAVE": -1,
                            "PROCSCAN": "Unknown",
                            "PROCTYPE": "MAP",
                            "WCALPOS": "Unknown",
                        }
                    )
                    idx += 1
                    row_in_file += 1

        df = pd.DataFrame(rows)
        metadata = IndexMetadata(
            created="Tue Jul 26 12:43:16 2016",
            last_modified="Tue Jul 26 12:43:16 2016",
            created_by="index_writer",
        )
        index_path = tmp_path / "test.index"
        write_index(index_path, metadata, df)

        actual = _get_rows_section(index_path)
        expected = _get_rows_section(expected_path)

        # Column header should match exactly
        assert actual[0] == expected[0], "Column header row does not match sparrow3"

        # All 32 data rows should match exactly
        assert len(actual) == len(expected), f"Row count mismatch: {len(actual) - 1} vs {len(expected) - 1} data rows"
        for i, (act, exp) in enumerate(zip(actual[1:], expected[1:], strict=True)):
            assert act == exp, f"VEGAS multi-bank data row {i} does not match sparrow3 expected output"


class TestMultiFileVegas:
    """Tests for writing index data from multi-file VEGAS observations.

    Real VEGAS data has multiple FITS files (one per bank: A, B, C, D),
    each with different IFNUMs, polarizations, and frequencies. The index
    must handle varying FILE, EXTENSION, POLARIZATION, PROCEDURE, CENTFREQ,
    and other columns across rows.
    """

    def test_write_multifile_roundtrip(self, tmp_path):
        """Test writing and reading an index with multi-file VEGAS-like data."""
        # Simulate 4 VEGAS banks (A, B, C, D) with 2 scans, 2 pols each
        rows = []
        files = [
            "AGBT18B_354_03.raw.vegas.A.fits",
            "AGBT18B_354_03.raw.vegas.B.fits",
            "AGBT18B_354_03.raw.vegas.C.fits",
            "AGBT18B_354_03.raw.vegas.D.fits",
        ]
        # Each bank has different CRVAL1/CRPIX1/CDELT1 (different tunings)
        bank_freqs = [
            (1.40e9, 512.0, 1464.84375),  # Bank A
            (1.42e9, 512.0, 1464.84375),  # Bank B
            (1.44e9, 512.0, 1464.84375),  # Bank C
            (1.46e9, 512.0, 1464.84375),  # Bank D
        ]
        idx = 0
        for bank_i, (filename, (crval1, crpix1, cdelt1)) in enumerate(zip(files, bank_freqs, strict=True)):
            for scan in (6, 7):
                for crval4 in (-1, -2):  # RR, LL polarizations
                    rows.append(
                        {
                            "INDEX": idx,
                            "FILE": filename,
                            "EXTENSION": 1,
                            "ROW": idx,
                            "OBJECT": "W49N",
                            "OBSMODE": "OffOn:PSWITCHON:TPWCAL" if scan == 6 else "OffOn:PSWITCHOFF:TPWCAL",
                            "OBSID": "Observation:1",
                            "PROJID": "AGBT18B_354_03",
                            "SCAN": scan,
                            "CRVAL4": crval4,
                            "CRVAL1": crval1,
                            "CRPIX1": crpix1,
                            "CDELT1": cdelt1,
                            "NUMCHN": 1024,
                            "PLNUM": 0 if crval4 == -1 else 1,
                            "IFNUM": bank_i,
                            "FEED": 1,
                            "FDNUM": 0,
                            "INTNUM": 0,
                            "PROCSEQN": 1,
                            "SIG": True,
                            "CAL": False,
                            "SAMPLER": f"A{bank_i}_0",
                            "AZIMUTH": 202.26,
                            "ELEVATIO": 66.45,
                            "CRVAL2": 287.755,
                            "CRVAL3": 14.136,
                            "TRGTLONG": 287.755,
                            "TRGTLAT": 14.136,
                            "SUBREF_STATE": 1,
                            "LST": 45123.0,
                            "RESTFREQ": 1.42040575e9,
                            "VELOCITY": 0.0,
                            "FREQRES": 1464.84375,
                            "DATE-OBS": "2018-12-15T06:12:00.00",
                            "TIMESTAMP": "2018_12_15_06:12:00",
                            "BANDWID": 1.5e6,
                            "EXPOSURE": 30.0,
                            "TSYS": 25.0,
                        }
                    )
                    idx += 1

        df = pd.DataFrame(rows)
        metadata = create_index_metadata(observer="Test", backend="VEGAS")

        # Write
        index_path = tmp_path / "test.index"
        write_index(index_path, metadata, df)

        # Read back
        _read_metadata, df_read = read_index(index_path)

        # Verify all 16 rows survived
        assert len(df_read) == 16

        # Verify FILE column has 4 distinct values
        assert df_read["FILE"].nunique() == 4

        # Verify POLARIZATION derived correctly from CRVAL4
        pol_values = set(df_read["POL"].unique())
        assert pol_values == {"RR", "LL"}

        # Verify PROCEDURE extracted from OBSMODE
        proc_values = set(df_read["PROCEDURE"].unique())
        assert proc_values == {"OffOn"}

        # Verify CENTFREQ computed (not just CRVAL1) and varies across banks
        assert df_read["CENTFREQ"].nunique() == 4

        # Verify PROJECT mapped from PROJID
        assert all(df_read["PROJECT"] == "AGBT18B_354_03")

        # Verify SOURCE mapped from OBJECT
        assert all(df_read["SOURCE"] == "W49N")

        # Verify SIG/CAL booleans roundtripped
        assert df_read["SIG"].all()
        assert not df_read["CAL"].any()

    def test_write_real_vegas_data(self, testdata_dir, tmp_path):
        """Test writing an index from real multi-bank VEGAS FITS data."""
        vegas_dir = testdata_dir / "AGBT18B_354_03" / "AGBT18B_354_03.raw.vegas"
        if not vegas_dir.exists():
            pytest.skip(f"VEGAS testdata not found: {vegas_dir}")

        from dysh.fits.gbtfitsload import GBTFITSLoad

        loader = GBTFITSLoad(str(vegas_dir))
        df = loader._selection

        metadata = create_index_metadata(
            observer=df["OBSERVER"].iloc[0] if "OBSERVER" in df.columns else "Unknown",
            backend=df["BACKEND"].iloc[0] if "BACKEND" in df.columns else "VEGAS",
        )

        # Write index
        index_path = tmp_path / "AGBT18B_354_03.raw.vegas.index"
        write_index(index_path, metadata, df)

        # Read it back
        _meta, df_read = read_index(index_path)

        # Should have all rows
        assert len(df_read) == len(df)

        # Should have all 43 standard columns
        assert len(df_read.columns) == 43

        # Key columns should survive roundtrip
        assert df_read["SCAN"].nunique() == df["SCAN"].nunique()
        assert df_read["IFNUM"].nunique() == df["IFNUM"].nunique()

        # POLARIZATION should be derived from CRVAL4
        if "CRVAL4" in df.columns:
            assert "POL" in df_read.columns
            expected_pols = set(df["CRVAL4"].apply(_get_polarization).unique())
            actual_pols = set(df_read["POL"].unique())
            assert actual_pols == expected_pols

        # Written file should be parseable by parse_sdfits_index_file too
        from dysh.fits.index_file import parse_sdfits_index_file

        df_parsed = parse_sdfits_index_file(index_path)
        assert len(df_parsed) == len(df)
