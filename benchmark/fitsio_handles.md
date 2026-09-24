# fitsio handle reuse: rationale and trade-offs

Notes from planning the performance-redux series (recovering closed PR #1049).
Decision: use **persistent handles** (one open `fitsio.FITS` handle per SDFITS file,
cached on `SDFITSLoad`).

## Why reopening is expensive

Every selective read in `src/dysh/fits/sdfitsload.py` currently does
`with fitsio.FITS(self._filename) as f:` — open the file, parse the primary header,
scan and parse the bintable extension headers to locate HDUs, read the requested rows,
close. That setup work repeats on every call even though the file has not changed.
The hot paths call it frequently:

- `rawspectra()` reopens per calibration call (every `getps`/`gettp`/`getsigref`)
- `nchan` opens the whole file just to inspect the DATA column dtype
- `load_full_rows()` reopens on every lazy-metadata fetch

On a local SSD with a warm page cache each open costs ~milliseconds — a real but modest
win. On NFS the cost is much larger: NFS close-to-open consistency semantics force
attribute revalidation with the server on each `open()` (lookup/getattr round-trips),
and header blocks are re-fetched rather than served from cache. Since GBO users
typically reduce data off NFS mounts, this hits the primary audience. PR #1049's
benchmarks diverged by host for exactly this reason: the SSD laptop improved modestly,
while `celeste` (NFS mounts) showed 4–10x overall gains. A survey reduction looping
`getps` over scans reopens the same file dozens to hundreds of times.

## Downsides of holding handles open

- **Stale handles on file change.** If the file is rewritten or appended while a handle
  is open, reads can miss new rows or see inconsistent data. Matters most for GBO
  "online" mode, where SDFITS files grow during an observation — the handle cache needs
  invalidation on refresh/write paths.
- **NFS ESTALE.** If the file is replaced server-side, the held handle errors with
  ESTALE; a reopen-and-retry path is needed.
- **File descriptor count.** One fd per SDFITS file for the life of the `GBTFITSLoad`.
  A multi-bank VEGAS session is typically 8 files — far below the default 1024 ulimit —
  but a batch pipeline holding many loaders could accumulate fds.
- **Serialization.** Open C-level handles cannot be pickled or deep-copied; any
  pickle/multiprocessing use of `SDFITSLoad` needs `__getstate__`/`__setstate__` to drop
  handles and lazily reopen.
- **Lifecycle complexity.** Requires an explicit `close()` (plus `__del__` fallback) —
  more code than the current self-cleaning `with` blocks.

None of these are blockers — invalidate on write/online-refresh, reopen on error,
exclude handles from pickling — but they are why this change is low-medium risk rather
than free. The LRU-pool alternative only addresses the fd-count concern (the least
serious one) at the cost of more code, which is why persistent handles were chosen.
