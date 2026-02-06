"""Lazy flag arrays for memory-efficient SDFITS channel masking.

Instead of eagerly allocating dense (nrows, nchan) boolean arrays for flags,
these classes store only modified rows in a sparse dictionary overlay.
For a 77GB file (~1.2M rows x 16384 channels), this reduces flag memory
from ~40GB to a few MB (proportional to flagged rows only).
"""

import tempfile

import fitsio
import numpy as np

from ..log import logger

# Number of rows to process at a time when materializing large arrays
CHUNK_SIZE = 10_000

# Threshold in bytes above which to_dense() uses memmap instead of RAM
MEMMAP_THRESHOLD = 1 * 1024**3  # 1 GB


class LazyFlagArray:
    """A sparse lazy wrapper around a (nrows, nchan) boolean flag array.

    Stores only rows that have been modified via `or_rows()` in a dict.
    Optionally reads base flags from a FITS FLAGS column on demand.

    Parameters
    ----------
    nrows : int
        Number of rows.
    nchan : int
        Number of channels per row.
    filename : str or None
        Path to FITS file for on-demand FLAGS column reads.
    hdu_index : int or None
        HDU index (1-based) in the FITS file.
    has_flags_column : bool
        Whether the FITS file has a FLAGS column to read from.
    memmap_threshold : int
        Byte threshold above which to_dense() uses numpy memmap.
    """

    def __init__(
        self,
        nrows,
        nchan,
        filename=None,
        hdu_index=None,
        has_flags_column=False,
        memmap_threshold=MEMMAP_THRESHOLD,
    ):
        self._nrows = nrows
        self._nchan = nchan
        self._filename = filename
        self._hdu_index = hdu_index
        self._has_flags_column = has_flags_column
        self._memmap_threshold = memmap_threshold
        # Sparse overlay: row_index -> 1D bool array of length nchan
        self._modified = {}
        # Broadcast masks: 1D masks that apply to ALL rows (avoids per-row copies)
        self._broadcasts = []

    @property
    def shape(self):
        return (self._nrows, self._nchan)

    @property
    def ndim(self):
        return 2

    @property
    def dtype(self):
        return np.dtype(bool)

    def __len__(self):
        return self._nrows

    def _read_flags_rows(self, rows):
        """Read FLAGS column for specific rows from the FITS file.

        Parameters
        ----------
        rows : array-like
            Row indices to read.

        Returns
        -------
        np.ndarray
            Boolean array of shape (len(rows), nchan).
        """
        if not self._has_flags_column or self._filename is None:
            return np.zeros((len(rows), self._nchan), dtype=bool)

        rows_array = np.asarray(rows)
        with fitsio.FITS(self._filename) as f:
            flags_data = f[self._hdu_index].read(columns=["FLAGS"], rows=rows_array)["FLAGS"]
        return flags_data.astype(bool)

    def __getitem__(self, key):
        """Return dense ndarray for the requested rows.

        Supports integer, slice, and array indexing on the first axis,
        and optional second-axis indexing (e.g., [rows, channels]).

        Returns
        -------
        np.ndarray
            Dense boolean array for the requested slice.
        """
        # Handle tuple indexing: (row_key, col_key)
        if isinstance(key, tuple):
            row_key = key[0]
            col_key = key[1] if len(key) > 1 else slice(None)
        else:
            row_key = key
            col_key = slice(None)

        # Normalize row_key to an array of row indices
        if isinstance(row_key, (int, np.integer)):
            row_indices = np.array([row_key])
            squeeze_row = True
        elif isinstance(row_key, slice):
            row_indices = np.arange(*row_key.indices(self._nrows))
            squeeze_row = False
        else:
            row_indices = np.asarray(row_key)
            squeeze_row = False

        # Build dense result for these rows
        result = self._read_flags_rows(row_indices)

        # Apply broadcast masks (same mask for all rows)
        for bmask in self._broadcasts:
            result |= bmask

        # Apply sparse overlay (per-row modifications)
        for i, row_idx in enumerate(row_indices):
            row_idx_int = int(row_idx)
            if row_idx_int in self._modified:
                result[i] |= self._modified[row_idx_int]

        # Apply column indexing
        if not (isinstance(col_key, slice) and col_key == slice(None)):
            result = result[:, col_key]

        # Squeeze single-row dimension if integer indexing
        if squeeze_row:
            result = result[0]

        return result

    def or_rows(self, rows, mask):
        """Apply OR operation to specific rows.

        This is the write operation that replaces `array[rows] |= mask`.

        Parameters
        ----------
        rows : array-like
            Row indices to modify.
        mask : np.ndarray
            Boolean mask to OR. Can be 1D (broadcast to all rows) or 2D (per-row).
        """
        rows_array = np.asarray(rows).ravel()
        mask = np.asarray(mask, dtype=bool)

        if mask.ndim == 1:
            # 1D mask: if applied to all rows, store as a broadcast
            if len(rows_array) == self._nrows:
                self._broadcasts.append(mask.copy())
                return
            # Otherwise store per-row
            for row_idx in rows_array:
                row_idx_int = int(row_idx)
                if row_idx_int in self._modified:
                    self._modified[row_idx_int] = self._modified[row_idx_int] | mask
                else:
                    # Need to read existing flags for this row to OR with
                    existing = self._read_flags_rows([row_idx_int])[0]
                    self._modified[row_idx_int] = existing | mask
        elif mask.ndim == 2:
            if mask.shape[0] != len(rows_array):
                raise ValueError(f"mask rows {mask.shape[0]} != number of rows {len(rows_array)}")
            for i, row_idx in enumerate(rows_array):
                row_idx_int = int(row_idx)
                if row_idx_int in self._modified:
                    self._modified[row_idx_int] = self._modified[row_idx_int] | mask[i]
                else:
                    existing = self._read_flags_rows([row_idx_int])[0]
                    self._modified[row_idx_int] = existing | mask[i]
        else:
            raise ValueError(f"mask must be 1D or 2D, got {mask.ndim}D")

    def merge_from(self, other):
        """Merge another LazyFlagArray's modifications into this one (sparse-to-sparse).

        Parameters
        ----------
        other : LazyFlagArray
            The array whose modifications to merge.
        """
        if not isinstance(other, LazyFlagArray):
            raise TypeError(f"Expected LazyFlagArray, got {type(other)}")

        # Merge broadcast masks
        self._broadcasts.extend(bmask.copy() for bmask in other._broadcasts)

        # If other has a FLAGS column backing and no modifications, nothing to merge
        # unless it has a flags column (which means base data that should be OR'd in)
        if other._has_flags_column:
            # The other array has base data from FITS; we need to read it for all rows
            # But this should be rare in practice (additional_channel_mask never has FITS backing)
            logger.warning("merge_from with FITS-backed other; reading all flags")
            for row_idx in range(other._nrows):
                row_data = other[row_idx]
                if np.any(row_data):
                    if row_idx in self._modified:
                        self._modified[row_idx] |= row_data
                    else:
                        existing = self._read_flags_rows([row_idx])[0]
                        self._modified[row_idx] = existing | row_data
        else:
            # Fast path: only merge the sparse overlay
            for row_idx, row_mask in other._modified.items():
                if row_idx in self._modified:
                    self._modified[row_idx] = self._modified[row_idx] | row_mask
                else:
                    existing = self._read_flags_rows([row_idx])[0]
                    self._modified[row_idx] = existing | row_mask

    def to_dense(self):
        """Materialize the full (nrows, nchan) boolean array.

        For large arrays (> memmap_threshold bytes), uses numpy memmap
        backed by a temp file to avoid RAM spikes.

        Returns
        -------
        np.ndarray
            Dense boolean array of shape (nrows, nchan).
        """
        estimated_bytes = self._nrows * self._nchan
        use_memmap = estimated_bytes > self._memmap_threshold

        if use_memmap:
            logger.info(f"LazyFlagArray.to_dense: using memmap for {estimated_bytes / 1024**3:.1f} GB array")
            tmpfile = tempfile.NamedTemporaryFile(suffix=".flags", delete=False)
            tmpfile.close()
            result = np.memmap(tmpfile.name, dtype=bool, mode="w+", shape=(self._nrows, self._nchan))
        else:
            result = np.zeros((self._nrows, self._nchan), dtype=bool)

        # Fill from FITS FLAGS column in chunks if it exists
        if self._has_flags_column and self._filename is not None:
            for chunk_start in range(0, self._nrows, CHUNK_SIZE):
                chunk_end = min(chunk_start + CHUNK_SIZE, self._nrows)
                chunk_rows = np.arange(chunk_start, chunk_end)
                flags_data = self._read_flags_rows(chunk_rows)
                result[chunk_start:chunk_end] = flags_data

        # Apply broadcast masks
        for bmask in self._broadcasts:
            result |= bmask

        # Apply sparse overlay
        for row_idx, row_mask in self._modified.items():
            result[row_idx] |= row_mask

        return result

    def rows_as_uint8(self, rows):
        """Build a uint8 flag array for specific rows, in chunks.

        Avoids materializing the full (nrows, nchan) dense array.
        Used by the write path to build the FLAGS column efficiently.

        Parameters
        ----------
        rows : array-like
            Row indices to materialize.

        Returns
        -------
        np.ndarray
            uint8 array of shape (len(rows), nchan).
        """
        rows_array = np.asarray(rows)
        nrows_out = len(rows_array)
        result = np.zeros((nrows_out, self._nchan), dtype=np.uint8)

        # Precompute combined broadcast mask as uint8
        if self._broadcasts:
            combined_broadcast = np.zeros(self._nchan, dtype=np.uint8)
            for bmask in self._broadcasts:
                combined_broadcast |= bmask.astype(np.uint8)
            result[:] = combined_broadcast

        # Fill from FITS FLAGS column in chunks if it exists
        if self._has_flags_column and self._filename is not None:
            for chunk_start in range(0, nrows_out, CHUNK_SIZE):
                chunk_end = min(chunk_start + CHUNK_SIZE, nrows_out)
                chunk_row_indices = rows_array[chunk_start:chunk_end]
                flags_data = self._read_flags_rows(chunk_row_indices)
                result[chunk_start:chunk_end] |= flags_data.astype(np.uint8)

        # Apply sparse overlay
        # Build a reverse lookup: original row index -> position in rows_array
        if self._modified:
            row_to_pos = {}
            for i, r in enumerate(rows_array):
                row_to_pos[int(r)] = i
            for row_idx, row_mask in self._modified.items():
                if row_idx in row_to_pos:
                    result[row_to_pos[row_idx]] |= row_mask.astype(np.uint8)

        return result

    def astype(self, dtype):
        """Convert to dense array with the given dtype.

        Parameters
        ----------
        dtype : numpy dtype
            Target dtype.

        Returns
        -------
        np.ndarray
            Dense array with the requested dtype.
        """
        return self.to_dense().astype(dtype)

    def __array__(self, dtype=None, copy=None):
        """Support numpy auto-conversion (e.g., when used as MaskedArray mask)."""
        result = self.to_dense()
        if dtype is not None:
            result = result.astype(dtype)
        return result

    def copy(self):
        """Return a deep copy of this LazyFlagArray."""
        new = LazyFlagArray(
            self._nrows,
            self._nchan,
            filename=self._filename,
            hdu_index=self._hdu_index,
            has_flags_column=self._has_flags_column,
            memmap_threshold=self._memmap_threshold,
        )
        new._modified = {k: v.copy() for k, v in self._modified.items()}
        new._broadcasts = [b.copy() for b in self._broadcasts]
        return new

    def __eq__(self, other):
        """Element-wise equality (materializes both arrays)."""
        if isinstance(other, LazyFlagArray):
            return np.array_equal(np.asarray(self), np.asarray(other))
        return np.asarray(self) == other

    def __or__(self, other):
        """Element-wise OR (materializes both arrays, returns dense)."""
        return np.asarray(self) | np.asarray(other)

    def __ior__(self, other):
        """In-place OR: merge other into this LazyFlagArray.

        If other is a LazyFlagArray, uses sparse merge.
        Otherwise materializes and ORs as dense.
        """
        if isinstance(other, LazyFlagArray):
            self.merge_from(other)
            return self
        # Dense fallback
        other_arr = np.asarray(other)
        if other_arr.ndim == 2:
            for i in range(other_arr.shape[0]):
                if np.any(other_arr[i]):
                    if i in self._modified:
                        self._modified[i] = self._modified[i] | other_arr[i]
                    else:
                        existing = self._read_flags_rows([i])[0]
                        self._modified[i] = existing | other_arr[i]
        return self

    def __repr__(self):
        n_modified = len(self._modified)
        n_broadcasts = len(self._broadcasts)
        backing = "FITS-backed" if self._has_flags_column else "zero-backed"
        parts = f"LazyFlagArray(nrows={self._nrows}, nchan={self._nchan}, {backing}, {n_modified} modified rows"
        if n_broadcasts:
            parts += f", {n_broadcasts} broadcasts"
        return parts + ")"


class LazyFlagContainer:
    """A container that replaces np.empty(N, dtype=object) for flag arrays.

    Indexed by bintable number; each element is a LazyFlagArray.

    Parameters
    ----------
    size : int
        Number of bintables.
    """

    def __init__(self, size):
        self._arrays = [None] * size

    def __getitem__(self, key):
        return self._arrays[key]

    def __setitem__(self, key, value):
        self._arrays[key] = value

    def __len__(self):
        return len(self._arrays)

    def __iter__(self):
        return iter(self._arrays)

    def __ior__(self, other):
        """In-place OR: merge another container's flags element-wise."""
        if isinstance(other, LazyFlagContainer):
            for i in range(len(self._arrays)):
                if self._arrays[i] is not None and other._arrays[i] is not None:
                    self._arrays[i] |= other._arrays[i]
        else:
            # Support numpy object array (legacy)
            for i in range(len(self._arrays)):
                if self._arrays[i] is not None and other[i] is not None:
                    self._arrays[i] |= other[i]
        return self

    def copy(self):
        """Return a deep copy of this container."""
        new = LazyFlagContainer(len(self._arrays))
        for i, arr in enumerate(self._arrays):
            if arr is not None:
                new._arrays[i] = arr.copy()
        return new

    def __repr__(self):
        parts = []
        for i, arr in enumerate(self._arrays):
            parts.append(f"  [{i}]: {arr!r}")
        return "LazyFlagContainer(\n" + "\n".join(parts) + "\n)"
