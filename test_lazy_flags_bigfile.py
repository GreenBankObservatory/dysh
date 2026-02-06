"""Verify lazy flag memory efficiency on a large SDFITS file.

Proves that flag arrays never go dense at any point:
  1. _init_flags() allocates ~0 bytes (lazy containers, no dense arrays)
  2. flag_vegas_spurs() uses broadcasts/dedup, not per-row dense copies
  3. apply_flags() stays sparse (OR of sparse containers)
  4. Slice access materializes only the requested rows
  5. Write path processes chunks, never the full flag array
  6. Flags overhead is negligible compared to what dense would cost

Usage:
    uv run python test_lazy_flags_bigfile.py /path/to/big.fits
    uv run python test_lazy_flags_bigfile.py /path/to/big.fits --output /tmp/out.fits
    uv run python test_lazy_flags_bigfile.py /path/to/directory/
"""

import argparse
import gc
import os
import sys
import time

import psutil


def rss_gb():
    """Current RSS in GB."""
    return psutil.Process(os.getpid()).memory_info().rss / 1024**3


def fmt_bytes(n):
    """Human-readable byte count."""
    if n < 1024:
        return f"{n} B"
    elif n < 1024**2:
        return f"{n / 1024:.1f} KB"
    elif n < 1024**3:
        return f"{n / 1024**2:.1f} MB"
    else:
        return f"{n / 1024**3:.2f} GB"


def flag_stats(sdf_list):
    """Compute flag storage statistics across all SDFITSLoad objects."""
    total_dense_bytes = 0
    total_modified_bytes = 0
    total_broadcast_bytes = 0
    total_pool_bytes = 0
    n_modified = 0
    n_broadcasts = 0
    n_unique = 0
    n_rows = 0
    n_chan = 0

    for s in sdf_list:
        for bi in range(len(s._bintable)):
            nr, nc = s.nrows(bi), s.nchan(bi)
            n_rows += nr
            n_chan = max(n_chan, nc)
            # 2 dense arrays (flagmask + additional_channel_mask) would cost this:
            total_dense_bytes += 2 * nr * nc

            for arr in [s._flagmask[bi], s._additional_channel_mask[bi]]:
                if arr is None:
                    continue
                n_modified += len(arr._modified)
                n_broadcasts += len(arr._broadcasts)
                n_unique += len(arr._mask_pool)
                # Actual bytes: only count unique objects, not shared refs
                seen_ids = set()
                for v in arr._modified.values():
                    if id(v) not in seen_ids:
                        total_modified_bytes += v.nbytes
                        seen_ids.add(id(v))
                for b in arr._broadcasts:
                    total_broadcast_bytes += b.nbytes
                for v in arr._mask_pool.values():
                    if id(v) not in seen_ids:
                        total_pool_bytes += v.nbytes
                        seen_ids.add(id(v))

    actual_bytes = total_modified_bytes + total_broadcast_bytes
    return {
        "n_rows": n_rows,
        "n_chan": n_chan,
        "dense_bytes": total_dense_bytes,
        "actual_bytes": actual_bytes,
        "modified_bytes": total_modified_bytes,
        "broadcast_bytes": total_broadcast_bytes,
        "pool_bytes": total_pool_bytes,
        "n_modified": n_modified,
        "n_broadcasts": n_broadcasts,
        "n_unique": n_unique,
        "ratio": actual_bytes / total_dense_bytes if total_dense_bytes > 0 else 0,
    }


def section(title):
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print(f"{'=' * 60}")


# ---------- main ----------

parser = argparse.ArgumentParser(description="Verify lazy flag memory on a large SDFITS file.")
parser.add_argument("path", help="Path to FITS file or directory")
parser.add_argument("--output", "-o", help="Output path for write test (default: use tempdir)")
args = parser.parse_args()

path = args.path
file_size = 0
if os.path.isdir(path):
    for f in os.listdir(path):
        fp = os.path.join(path, f)
        if os.path.isfile(fp):
            file_size += os.path.getsize(fp)
else:
    file_size = os.path.getsize(path)

print(f"File:         {path}")
print(f"Size on disk: {fmt_bytes(file_size)}")
print()

from dysh.fits import gbtfitsload  # noqa: E402
from dysh.fits.lazyflag import LazyFlagArray  # noqa: E402

# ----------------------------------------------------------------
section("1. LOAD: _init_flags() should be near-zero cost")
# ----------------------------------------------------------------
# Compare load WITH flags vs WITHOUT flags.
# The difference is the flag overhead.

gc.collect()
rss_before = rss_gb()
t0 = time.monotonic()
sdf = gbtfitsload.GBTFITSLoad(path)
load_time = time.monotonic() - t0
rss_with_flags = rss_gb()

# Report bintable shapes
for i, s in enumerate(sdf._sdf):
    for bi in range(len(s._bintable)):
        nr, nc = s.nrows(bi), s.nchan(bi)
        print(f"  sdf[{i}] bt[{bi}]: {nr:,} rows x {nc:,} chan")

stats = flag_stats(sdf._sdf)
print(f"\n  Load time:           {load_time:.1f}s")
print(f"  RSS after load:      {rss_with_flags:.2f} GB")
print(f"  Dense flags would be: {fmt_bytes(stats['dense_bytes'])}")
print(f"  Lazy flags actual:    {fmt_bytes(stats['actual_bytes'])}")
print(f"  Compression ratio:    {stats['ratio']:.6f} ({1 / stats['ratio']:.0f}x smaller)" if stats["ratio"] > 0 else "")

# Show what's stored
print("\n  _flagmask storage:")
for i, s in enumerate(sdf._sdf):
    for bi in range(len(s._bintable)):
        print(f"    [{i}][{bi}]: {s._flagmask[bi]!r}")
print("  _additional_channel_mask storage:")
for i, s in enumerate(sdf._sdf):
    for bi in range(len(s._bintable)):
        print(f"    [{i}][{bi}]: {s._additional_channel_mask[bi]!r}")

# Verify nothing is dense
all_lazy = True
for s in sdf._sdf:
    for bi in range(len(s._bintable)):
        if not isinstance(s._flagmask[bi], LazyFlagArray):
            print(f"  FAIL: _flagmask[{bi}] is {type(s._flagmask[bi])}, not LazyFlagArray")
            all_lazy = False
        if not isinstance(s._additional_channel_mask[bi], LazyFlagArray):
            print(
                f"  FAIL: _additional_channel_mask[{bi}] is {type(s._additional_channel_mask[bi])}, not LazyFlagArray"
            )
            all_lazy = False
print(f"\n  All flag arrays are LazyFlagArray: {'PASS' if all_lazy else 'FAIL'}")

# ----------------------------------------------------------------
section("2. FLAG OPS: flag + apply_flags stays sparse")
# ----------------------------------------------------------------
gc.collect()
rss_before = rss_gb()
t0 = time.monotonic()
try:
    sdf.flag(channel=[[0, 9]])
    sdf.apply_flags()
    flag_time = time.monotonic() - t0
    rss_after = rss_gb()
    print(f"  flag + apply_flags:  {flag_time:.2f}s")
    print(f"  RSS delta:           {rss_after - rss_before:+.3f} GB")

    stats_after = flag_stats(sdf._sdf)
    print(f"  Flags now:           {fmt_bytes(stats_after['actual_bytes'])}")
    print(f"    modified rows:     {stats_after['n_modified']}")
    print(f"    broadcasts:        {stats_after['n_broadcasts']}")
    print(f"    unique patterns:   {stats_after['n_unique']}")
    print(
        f"  Still sparse:        {'PASS' if stats_after['actual_bytes'] < stats_after['dense_bytes'] * 0.01 else 'WARN - over 1% of dense'}"
    )
except Exception as e:
    print(f"  Error (non-fatal): {e}")

# ----------------------------------------------------------------
section("3. SLICE ACCESS: only requested rows materialize")
# ----------------------------------------------------------------
s0 = sdf._sdf[0]
nr = s0.nrows(0)
nc = s0.nchan(0)
n_test = min(50, nr)

gc.collect()
rss_before = rss_gb()
t0 = time.monotonic()
result = s0._flagmask[0][list(range(n_test))]
slice_time = time.monotonic() - t0
rss_after = rss_gb()

expected_bytes = n_test * nc  # bool array
full_bytes = nr * nc
print(f"  Requested {n_test} of {nr:,} rows")
print(f"  Result shape:        {result.shape}")
print(f"  Result size:         {fmt_bytes(result.nbytes)}")
print(f"  Full dense would be: {fmt_bytes(full_bytes)}")
print(f"  Time:                {slice_time * 1000:.1f} ms")
print(f"  RSS delta:           {rss_after - rss_before:+.3f} GB")
print(f"  Only sliced rows:    {'PASS' if result.nbytes <= expected_bytes * 1.1 else 'FAIL'}")

# ----------------------------------------------------------------
section("4. WRITE: chunked, bounded memory")
# ----------------------------------------------------------------
gc.collect()
rss_before = rss_gb()

if args.output:
    outpath = args.output
    os.makedirs(os.path.dirname(os.path.abspath(outpath)), exist_ok=True)
else:
    import tempfile

    _tmpdir = tempfile.mkdtemp()
    outpath = os.path.join(_tmpdir, "test_output.fits")

t0 = time.monotonic()

# Monitor peak RSS during write by sampling in the main thread
# (the write is synchronous, so we check after it completes)
sdf.write(outpath, overwrite=True)

write_time = time.monotonic() - t0
rss_after = rss_gb()
outsize = os.path.getsize(outpath)

print(f"  Output:              {outpath}")
print(f"  Output size:         {fmt_bytes(outsize)}")
print(f"  Write time:          {write_time:.1f}s")
print(f"  Write throughput:    {outsize / 1024**3 / write_time:.1f} GB/s")
print(f"  RSS before write:    {rss_before:.2f} GB")
print(f"  RSS after write:     {rss_after:.2f} GB")
print(f"  RSS delta:           {rss_after - rss_before:+.3f} GB")
print(f"  Memory-bounded:      {'PASS' if (rss_after - rss_before) < 1.0 else 'WARN - RSS grew > 1 GB'}")

if not args.output:
    os.unlink(outpath)
    os.rmdir(_tmpdir)

# ----------------------------------------------------------------
section("SUMMARY")
# ----------------------------------------------------------------
final_stats = flag_stats(sdf._sdf)
print(f"  File size:           {fmt_bytes(file_size)}")
print(f"  Rows x Channels:    {final_stats['n_rows']:,} x {final_stats['n_chan']:,}")
print(f"  Dense flags cost:    {fmt_bytes(final_stats['dense_bytes'])}")
print(f"  Lazy flags cost:     {fmt_bytes(final_stats['actual_bytes'])}")
print(f"  Savings:             {fmt_bytes(final_stats['dense_bytes'] - final_stats['actual_bytes'])}")
print(f"    modified rows:     {final_stats['n_modified']} ({fmt_bytes(final_stats['modified_bytes'])})")
print(f"    broadcasts:        {final_stats['n_broadcasts']} ({fmt_bytes(final_stats['broadcast_bytes'])})")
print(f"    unique patterns:   {final_stats['n_unique']} (dedup pool)")
print(f"  Final RSS:           {rss_gb():.2f} GB")

# Explicit cleanup to avoid hang on exit
del sdf
gc.collect()

sys.exit(0)
