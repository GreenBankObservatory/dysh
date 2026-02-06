"""Test lazy flag memory on a large SDFITS file.

Usage:
    uv run python test_lazy_flags_bigfile.py /path/to/big.fits
    uv run python test_lazy_flags_bigfile.py /path/to/directory/
"""

import os
import sys
import tracemalloc

import psutil


def mem_gb():
    return psutil.Process(os.getpid()).memory_info().rss / 1024**3


if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <fits_file_or_directory>")
    sys.exit(1)

path = sys.argv[1]
file_size = 0
if os.path.isdir(path):
    for f in os.listdir(path):
        fp = os.path.join(path, f)
        if os.path.isfile(fp):
            file_size += os.path.getsize(fp)
else:
    file_size = os.path.getsize(path)

print(f"File: {path}")
print(f"Size on disk: {file_size / 1024**3:.1f} GB")
print(f"RSS before import: {mem_gb():.2f} GB")
print()

from dysh.fits import gbtfitsload  # noqa: E402

print(f"RSS after import: {mem_gb():.2f} GB")

# --- Load with flags (lazy) ---
print("\n--- Loading with flags (lazy) ---")
tracemalloc.start()
rss_before = mem_gb()
snap1 = tracemalloc.take_snapshot()

sdf = gbtfitsload.GBTFITSLoad(path)

snap2 = tracemalloc.take_snapshot()
rss_after = mem_gb()
alloc_with = sum(s.size for s in snap2.compare_to(snap1, "lineno") if s.size > 0)
tracemalloc.stop()

print(f"RSS: {rss_before:.2f} -> {rss_after:.2f} GB (delta {rss_after - rss_before:.2f} GB)")
print(f"tracemalloc delta: {alloc_with / 1024**2:.1f} MB")

# Report bintable shapes
for i, s in enumerate(sdf._sdf):
    for bi in range(len(s._bintable)):
        nr, nc = s.nrows(bi), s.nchan(bi)
        dense_gb = 2 * nr * nc / 1024**3
        print(f"  sdf[{i}] bintable {bi}: {nr} rows x {nc} chan (dense would be {dense_gb:.1f} GB)")
        fm = s._flagmask[bi]
        print(f"    _flagmask: {fm!r}")
        print(f"    _additional_channel_mask: {s._additional_channel_mask[bi]!r}")

# --- Test flag operations ---
print("\n--- Testing flag operations ---")
rss_before = mem_gb()

# Try flagging some channels
try:
    sdf.flag(channel=[[0, 9]])
    sdf.apply_flags()
    print(f"flag + apply_flags: RSS {rss_before:.2f} -> {mem_gb():.2f} GB")
except Exception as e:
    print(f"flag/apply_flags error (non-fatal): {e}")

# --- Test small slice access ---
print("\n--- Testing slice access (50 rows) ---")
rss_before = mem_gb()
tracemalloc.start()
snap1 = tracemalloc.take_snapshot()

s0 = sdf._sdf[0]
result = s0._flagmask[0][list(range(min(50, s0.nrows(0))))]

snap2 = tracemalloc.take_snapshot()
slice_alloc = sum(s.size for s in snap2.compare_to(snap1, "lineno") if s.size > 0)
tracemalloc.stop()

print(f"Slice shape: {result.shape}")
print(f"Slice alloc: {slice_alloc / 1024:.1f} KB")
print(f"RSS: {rss_before:.2f} -> {mem_gb():.2f} GB")

# --- Summary ---
print("\n" + "=" * 50)
print("SUMMARY")
print("=" * 50)
print(f"File size:        {file_size / 1024**3:.1f} GB")
print(f"Final RSS:        {mem_gb():.2f} GB")
total_dense = sum(2 * s.nrows(bi) * s.nchan(bi) for s in sdf._sdf for bi in range(len(s._bintable)))
print(f"Dense flags would be: {total_dense / 1024**3:.1f} GB")
total_sparse = sum(
    sum(v.nbytes for v in s._flagmask[bi]._modified.values()) for s in sdf._sdf for bi in range(len(s._bintable))
)
print(f"Lazy flags actual:    {total_sparse / 1024**2:.1f} MB")
