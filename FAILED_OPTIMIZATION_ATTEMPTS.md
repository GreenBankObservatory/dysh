# Failed Optimization Attempts

This document records optimization attempts that were tried during the `bench` branch performance work and did not hold up. The goal is to preserve negative results so we do not repeat the same experiments without new evidence.

It covers:
- `hi_survey` science-path work
- startup/import salvage attempts
- lazy-loading / selection-path experiments
- benchmark-side debugging detours

It does not include changes that were kept and are now part of `bench` branch.

## Summary

The general pattern was:
- broad profiling identified real hotspot families
- some straightforward fixes worked and were kept
- many later attempts targeted real hotspots but lost on implementation details, hidden pandas costs, or correctness

The most important conclusion is that the remaining `hi_survey` cost is concentrated in calibration-path metadata/selection overhead, but several intuitive ways of reducing that overhead turned out to be worse in practice.

## 1. `getsigref` / `hi_survey` Attempts That Regressed

### 1.1 Match GBTIDL's default `avgref=False` behavior

Attempt:
- change `dysh.getsigref(ref=<scan int>)` to avoid pre-averaging the reference unless needed, mirroring GBTIDL's default branch

Why it looked promising:
- GBTIDL only pre-averages the reference in certain cases
- `hi_survey` has matched signal/reference integration counts, so using per-integration refs looked closer to the GBTIDL algorithm

What happened:
- local `hi_survey` script time regressed significantly
- the change was reverted

Takeaway:
- matching GBTIDL's branch structure is not automatically a win in `dysh`
- `dysh`'s current PS path is optimized around an averaged reference representation, so a more faithful algorithmic match still needs a cheap internal container to pay off

### 1.2 Lighter `use_wcs=False` spectrum construction in the `getsigref` path

Attempt:
- make `Spectrum.make_spectrum()` skip more target/observer/WCS-adjacent work when `use_wcs=False`

Why it looked promising:
- profiles showed nontrivial time in spectrum construction and coordinate-related setup
- `getsigref` and internal TP references often only need intermediate spectra

What happened:
- local `hi_survey` script time regressed
- the change was reverted

Takeaway:
- the expensive-looking coordinate path was not the dominant remaining cost at that point
- removing it from the wrong layer added overhead or broke fast assumptions elsewhere

### 1.3 NumPy indexing rewrite for TP/PS exposure and delta-freq setup

Attempt:
- reduce pandas overhead in scan initialization by replacing some row-selection work with NumPy-oriented indexing

Why it looked promising:
- `_calc_exposure()` and `_calc_delta_freq()` were visible in profiles

What happened:
- no meaningful improvement in the real script
- not kept

Takeaway:
- those helpers were not the right target
- their overhead was smaller than the surrounding metadata and selection churn

### 1.4 Direct internal TP-reference builder for `getsigref`

Attempt:
- bypass parts of `gettp().timeaverage()` for the TP reference used inside `getsigref(ref=<scan int>)`

Why it looked promising:
- nested `gettp()` inside `getsigref` remained visible in profiles
- GBTIDL uses much lighter reference containers internally

What happened:
- edge cases became messy quickly
- performance either did not improve or regressed
- the attempt was abandoned

Takeaway:
- there is still likely a real win here, but it needs a principled internal container design rather than an ad hoc shortcut

### 1.5 Skip coordinate precompute for internal TP references

Attempt:
- add `precompute_coordinates` plumbing through `gettp` / `TPScan` / `ScanBase` so the internal TP reference path in `getsigref` could skip observer/target precompute

Why it looked promising:
- target/observer precompute remained visible after earlier wins

What happened:
- the first version broke tests because of bad signature placement
- after fixing the API plumbing, the hot-call profile still regressed badly
- the whole experiment was reverted

Takeaway:
- the coordinate precompute path was not expensive enough to justify the extra branching and complexity

## 2. Lazy-Load / Selection Experiments That Regressed

### 2.1 In-place merged selection and flag updates during partial lazy loads

Attempt:
- avoid full rebuilds after partial lazy loads by mutating `Selection` / `Flag` state in place
- introduce helpers like `_merge_lazy_loaded_rows()` and `_empty_lazy_column()`

Why it looked promising:
- `_rebuild_merged_index()` and selection rebuild showed up clearly in hot-call profiles

What happened:
- pandas `__setitem__` storms dominated the profile
- hot `getsigref(...).timeaverage()` calls regressed badly
- the change was fully reverted

Takeaway:
- this was the clearest example of a real hotspot with the wrong implementation strategy
- mutating the merged structures in place cost more than rebuilding them

### 2.2 Lightweight selection rebuild helper

Attempt:
- add a `_lightweight_rebuild()` path in `selection.py` and use it for partial lazy-load refreshes

Why it looked promising:
- full `Selection` reconstruction looked expensive

What happened:
- tests passed, but the real hot-call profile regressed sharply
- the change was reverted

Takeaway:
- "lightweight" at the API level did not translate into cheaper actual pandas work

### 2.3 Trim vane-only metadata from generic calibration lazy loads

Attempt:
- reduce the generic calibration metadata subset by removing `TWARM` / `TAMBIENT`

Why it looked promising:
- smaller metadata subsets should mean less lazy-load overhead

What happened:
- Argus `vanecal()` broke
- even where it did not break, the hot path did not improve enough to justify the risk
- the change was reverted

Takeaway:
- calibration metadata dependencies are broader than they first appear
- the subset list must be treated as a correctness boundary, not a best-effort optimization

### 2.4 Combined signal/ref lazy load in `getsigref(ref=int)`

Attempt:
- combine signal and reference lazy loads to reduce repeated loader overhead

Why it looked promising:
- `_load_full_rows_if_needed()` was still visible in the hot path

What happened:
- `test_getsigref` regressed
- the 5-run script mean got much worse
- the experiment was reverted

Takeaway:
- combining loads changed selection/caching behavior in a way that lost more than it gained

### 2.5 Single-pass `pd.concat()` rebuild in `_rebuild_merged_index()`

Attempt:
- replace the incremental merged-index rebuild loop with a single `pd.concat(frames, copy=False)` over all underlying `_index` frames

Why it looked promising:
- `_rebuild_merged_index()` remained visible in the hot `getsigref` path
- the existing implementation does repeated concat work and copies the first frame eagerly

What happened:
- targeted tests passed, but the hot `getsigref(...).timeaverage()` cProfile total regressed slightly
- it also introduced a `Pandas4Warning` around the `copy` keyword
- the change was reverted

Takeaway:
- the obvious concat cleanup did not translate into a real win on this workload
- `_rebuild_merged_index()` is still a real cost center, but this specific rewrite is not worth keeping

## 3. Startup-Branch Ideas That Did Not Transfer Cleanly

### 3.1 Deferred IERS initialization on top of `bench`

Attempt:
- salvage the lazy IERS initialization pattern from `startup-bench-optimize`

Why it looked promising:
- `import dysh` spent noticeable time in Astropy IERS setup

What happened:
- import time improved a lot
- but `hi_survey` script time regressed slightly because the cost moved into runtime
- only CLI lazy imports were kept; deferred IERS was not

Takeaway:
- startup wins that simply defer work into the first real science path are not helpful for the benchmark we care about

### 3.2 Wholesale reuse of `startup-bench-optimize`

Attempt:
- compare or salvage the startup branch as a possible performance base

Why it looked promising:
- large total runtime in `dysh` was clearly startup-heavy

What happened:
- `startup-bench-optimize` was much worse than `bench` for `hi_survey`
- it also carried older, slower science-path code and headless/UI problems

Takeaway:
- only narrow startup-only pieces should ever be extracted from that branch
- the branch as a whole is not a useful base for `hi_survey`

## 4. Correctness Regressions Triggered by Optimization Work

These were not all pure "perf regressions", but they were negative outcomes from the optimization passes and are worth recording.

### 4.1 Argus lazy-load regressions

Observed failures:
- `KeyError: 'TWARM'`
- `AttributeError` / bad `RADESYS`
- stale `_fully_loaded_columns` state causing later feed failures

What caused them:
- metadata subset trimming was too aggressive
- returned rows were stale after post-load fixups
- global "fully loaded" state was treated as row-level truth

Outcome:
- fixed with targeted regressions and commit groups
- the failed variants should not be retried without tighter per-row invariants

### 4.2 GWCS slicing crash from `use_wcs=False` TP references

Observed failure:
- `AttributeError: 'SpectralGWCS' object has no attribute 'wcs'` in `Spectrum.__getitem__`

What caused it:
- internal TP references started coming through a `use_wcs=False` path
- slicing code assumed a FITS-style `.wcs` object

Outcome:
- this was a real bug in slicing assumptions
- it demonstrated that some performance work exposed latent correctness bugs outside the original hotspot

### 4.3 Broad metadata regressions in `bench`

Observed failures during larger test runs:
- `KeyError: 'FEEDXOFF'`
- `KeyError: 'BANDWID'`
- nod/calibration regressions
- notebook and plotting fallout

What this means:
- some intermediate branch states were not generally safe even if the narrow benchmark improved

Outcome:
- this reinforced that metadata subset tuning must be checked against broader calibration and plotting tests

## 5. Benchmark / Measurement Detours That Did Not Produce a Keepable Win

### 5.1 Chasing a mixed benchmark number before splitting startup vs script time

Issue:
- early comparisons used a single subprocess wall time

Why it misled:
- it mixed startup/import costs with actual reduction time

Outcome:
- the benchmark runner was later improved to report total vs script vs startup separately
- this was a measurement correction, not an optimization

### 5.2 Argus benchmark mismatch before script restructuring

Issue:
- the original `dysh` Argus script called `vanecal()` per feed in a way that did more setup work than the GBTIDL benchmark

Outcome:
- the benchmark was restructured so `tcal` setup happens once
- this exposed real relative performance better

Negative lesson:
- benchmark structure itself can fabricate a performance gap

### 5.3 Broken GBTIDL benchmark output and marker formatting

Issue:
- at one point Argus GBTIDL verification failed because helper procedures were not being compiled correctly and output markers were malformed

Outcome:
- fixed on the benchmark side

Negative lesson:
- do not trust missing markers or missing verify fields until raw stdout has been inspected

## 6. Existing Issue Exposed, Not Introduced: `hi_survey` Baseline Regions

This is included here because it looked like a performance-branch regression for a while, but it was actually a separate pre-existing correctness issue.

Observed problem:
- `hi_survey` RMS values disagreed strongly with GBTIDL even when the branch was otherwise fast

What we found:
- the frequency-space `baseline(include=...)` regions in the benchmark script behaved poorly on the descending axis
- converting those same regions to channel pairs produced dysh results much closer to GBTIDL

Outcome:
- documented separately in `HI_SURVEY_BASELINE_ISSUE.md`

Takeaway:
- do not attribute every benchmark mismatch to the latest performance change
- some mismatches are pre-existing script/library semantics gaps

## 7. Branch / Workflow Mistakes Worth Recording

### 7.1 Working in the wrong branch/worktree

This happened more than once during the session:
- work was accidentally performed on `startup-bench-optimize` or other non-target branches
- fixes and instrumentation had to be disentangled from branch drift

Takeaway:
- verify branch/worktree before substantial profiling or edits
- especially when benchmark branches intentionally diverge in behavior

### 7.2 Benchmark repo drift

At several points, local `dysh-gbtidl-bench` state did not match the latest committed benchmark scripts:
- missing stage labels
- stale debug changes
- local-only behavior differences

Takeaway:
- always confirm nested repo HEAD when comparing results across runs

## 8. What Not To Retry Blindly

Do not retry these without new profiling evidence or a different implementation plan:

- in-place merged selection/flag mutation for partial lazy loads
- "lightweight" selection rebuild paths that still lean on pandas mutation
- broad metadata subset trimming without end-to-end calibration coverage
- startup wins that simply defer cost into the first science-path call
- GBTIDL-algorithm parity changes without matching `dysh` internal data structures
- coordinate-precompute skipping in TP reference paths

## 9. What These Failures Still Taught Us

Even the failed attempts narrowed the search space:

- The remaining `hi_survey` science-path cost is still mostly calibration-path metadata and selection overhead.
- `_load_full_rows_if_needed()` and its rebuild side effects are real targets, but naive partial mutation strategies lose.
- Some apparent hotspots, especially target/observer precompute, are not worth special-casing at the current scale.
- Benchmark correctness and benchmark structure matter as much as raw timing when comparing against GBTIDL.
