"""Memory profiling tests for flag array initialization.

Run with: uv run pytest src/dysh/fits/tests/test_flag_memory.py -xvs
"""

import tracemalloc

import numpy as np

from dysh import util
from dysh.fits import gbtfitsload
from dysh.fits.lazyflag import LazyFlagArray


def test_flag_memory_baseline():
    """Measure memory consumed by flag array initialization."""
    fits_path = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas"

    # Measure WITH flags
    tracemalloc.start()
    snap_before = tracemalloc.take_snapshot()
    sdf = gbtfitsload.GBTFITSLoad(fits_path)
    snap_after = tracemalloc.take_snapshot()

    stats = snap_after.compare_to(snap_before, "lineno")
    total_with = sum(s.size for s in stats if s.size > 0)

    # Report nrows x nchan per bintable
    for s in sdf._sdf:
        for bi in range(len(s._bintable)):
            nr, nc = s.nrows(bi), s.nchan(bi)
            expected_dense = 2 * nr * nc  # 2 bool arrays (old approach)
            print(f"  bintable {bi}: {nr}x{nc} = {expected_dense / 1024**2:.1f} MB (would be dense)")

    # Measure WITHOUT flags
    snap_before2 = tracemalloc.take_snapshot()
    sdf2 = gbtfitsload.GBTFITSLoad(fits_path, skipflags=True)  # noqa: F841
    snap_after2 = tracemalloc.take_snapshot()
    total_without = sum(s.size for s in snap_after2.compare_to(snap_before2, "lineno") if s.size > 0)

    flag_overhead = total_with - total_without
    print(f"\nWith flags:    {total_with / 1024**2:.1f} MB")
    print(f"Without flags: {total_without / 1024**2:.1f} MB")
    print(f"Flag overhead: {flag_overhead / 1024**2:.1f} MB")

    tracemalloc.stop()

    # With lazy flags, the overhead should be very small (< 1 MB for test data)
    # The old dense approach would allocate 2 * nrows * nchan bytes per bintable
    assert flag_overhead < 5 * 1024**2, f"Flag overhead {flag_overhead / 1024**2:.1f} MB is unexpectedly large"


def test_lazy_flag_no_large_alloc():
    """Verify LazyFlagArray doesn't allocate large arrays on init."""
    # Simulate a 77GB file: 1.2M rows x 16K channels
    # Old approach: 2 x 1.2M x 16K = ~40 GB
    # New approach: ~200 bytes
    arr = LazyFlagArray(1_200_000, 16384)
    assert arr.shape == (1_200_000, 16384)
    assert len(arr._modified) == 0

    # Accessing a small slice should only allocate for those rows
    tracemalloc.start()
    snap1 = tracemalloc.take_snapshot()
    result = arr[np.array([0, 1, 2])]  # 3 rows
    snap2 = tracemalloc.take_snapshot()
    tracemalloc.stop()

    assert result.shape == (3, 16384)
    alloc = sum(s.size for s in snap2.compare_to(snap1, "lineno") if s.size > 0)
    # 3 rows x 16K = 48KB + overhead, should be well under 1MB
    assert alloc < 1 * 1024**2, f"Slice allocation {alloc / 1024:.1f} KB is too large"


def test_or_rows_sparse_memory():
    """Flagging individual rows should only store those rows."""
    arr = LazyFlagArray(1_000_000, 16384)

    # Flag 100 rows
    for i in range(100):
        arr.or_rows([i * 1000], np.ones(16384, dtype=bool))

    assert len(arr._modified) == 100
    # 100 rows x 16K channels x 1 byte = 1.6 MB + dict overhead
    total_bytes = sum(v.nbytes for v in arr._modified.values())
    assert total_bytes < 2 * 1024**2, f"Sparse storage {total_bytes / 1024**2:.1f} MB is too large"
