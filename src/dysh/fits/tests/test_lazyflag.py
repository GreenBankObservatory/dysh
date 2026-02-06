"""Tests for LazyFlagArray and LazyFlagContainer."""

import numpy as np
import pytest

from dysh.fits.lazyflag import LazyFlagArray, LazyFlagContainer


class TestLazyFlagArray:
    """Tests for the LazyFlagArray class."""

    def test_empty_is_all_false(self):
        """An unmodified LazyFlagArray should return all False."""
        arr = LazyFlagArray(100, 256)
        result = arr[0]
        assert result.shape == (256,)
        assert not np.any(result)

    def test_empty_slice_all_false(self):
        """Slicing all rows of an unmodified array returns all False."""
        arr = LazyFlagArray(10, 32)
        result = arr[:]
        assert result.shape == (10, 32)
        assert not np.any(result)

    def test_shape_property(self):
        arr = LazyFlagArray(50, 128)
        assert arr.shape == (50, 128)
        assert len(arr) == 50
        assert arr.ndim == 2
        assert arr.dtype == np.dtype(bool)

    def test_or_rows_1d_broadcast(self):
        """A 1D channel mask applied to multiple rows."""
        arr = LazyFlagArray(10, 8)
        mask = np.array([True, False, True, False, False, False, False, True])
        rows = np.array([2, 5, 7])
        arr.or_rows(rows, mask)

        # Check flagged rows
        for r in rows:
            result = arr[r]
            np.testing.assert_array_equal(result, mask)

        # Check unflagged row
        result = arr[0]
        assert not np.any(result)

    def test_or_rows_2d(self):
        """Per-row 2D masks."""
        arr = LazyFlagArray(10, 4)
        rows = np.array([1, 3])
        mask = np.array(
            [
                [True, False, False, True],
                [False, True, True, False],
            ]
        )
        arr.or_rows(rows, mask)

        np.testing.assert_array_equal(arr[1], [True, False, False, True])
        np.testing.assert_array_equal(arr[3], [False, True, True, False])
        assert not np.any(arr[0])

    def test_or_rows_accumulates(self):
        """Multiple or_rows calls should accumulate (OR) flags."""
        arr = LazyFlagArray(5, 4)
        arr.or_rows([0], np.array([True, False, False, False]))
        arr.or_rows([0], np.array([False, False, True, False]))
        np.testing.assert_array_equal(arr[0], [True, False, True, False])

    def test_getitem_single_row(self):
        arr = LazyFlagArray(5, 3)
        arr.or_rows([2], np.array([True, True, False]))
        result = arr[2]
        assert result.shape == (3,)
        np.testing.assert_array_equal(result, [True, True, False])

    def test_getitem_array(self):
        arr = LazyFlagArray(10, 4)
        arr.or_rows([1], np.array([True, False, False, False]))
        arr.or_rows([3], np.array([False, True, False, False]))

        result = arr[np.array([1, 3])]
        assert result.shape == (2, 4)
        np.testing.assert_array_equal(result[0], [True, False, False, False])
        np.testing.assert_array_equal(result[1], [False, True, False, False])

    def test_getitem_slice(self):
        arr = LazyFlagArray(10, 4)
        arr.or_rows([2], np.array([True, False, False, True]))
        result = arr[1:4]
        assert result.shape == (3, 4)
        np.testing.assert_array_equal(result[1], [True, False, False, True])  # row 2
        assert not np.any(result[0])  # row 1
        assert not np.any(result[2])  # row 3

    def test_getitem_tuple_row_col(self):
        """Test (rows, columns) tuple indexing."""
        arr = LazyFlagArray(10, 8)
        arr.or_rows([3], np.array([True, False, True, False, True, False, True, False]))
        result = arr[3, :4]
        np.testing.assert_array_equal(result, [True, False, True, False])

    def test_getitem_rows_and_col_slice(self):
        """Test [rows, col_slice] pattern used by _check_no_data_to_calibrate."""
        arr = LazyFlagArray(10, 8)
        arr.or_rows([2, 4], np.array([True] * 8))
        result = arr[np.array([2, 4]), slice(2, 6)]
        assert result.shape == (2, 4)
        assert np.all(result)

    def test_to_dense(self):
        arr = LazyFlagArray(5, 3)
        arr.or_rows([1], np.array([True, False, True]))
        dense = arr.to_dense()
        assert dense.shape == (5, 3)
        assert dense.dtype == bool
        np.testing.assert_array_equal(dense[1], [True, False, True])
        assert not np.any(dense[0])

    def test_astype(self):
        arr = LazyFlagArray(3, 2)
        arr.or_rows([0], np.array([True, False]))
        result = arr.astype(np.uint8)
        assert result.dtype == np.uint8
        np.testing.assert_array_equal(result[0], [1, 0])

    def test_array_protocol(self):
        """Test __array__ for numpy auto-conversion."""
        arr = LazyFlagArray(3, 2)
        arr.or_rows([1], np.array([True, True]))
        result = np.asarray(arr)
        assert result.shape == (3, 2)
        assert result.dtype == bool

    def test_copy(self):
        arr = LazyFlagArray(5, 3)
        arr.or_rows([2], np.array([True, False, True]))
        arr2 = arr.copy()
        # Modify original
        arr.or_rows([2], np.array([False, True, False]))
        # Copy should be unchanged
        np.testing.assert_array_equal(arr2[2], [True, False, True])
        np.testing.assert_array_equal(arr[2], [True, True, True])

    def test_eq(self):
        arr1 = LazyFlagArray(3, 2)
        arr2 = LazyFlagArray(3, 2)
        assert arr1 == arr2

        arr1.or_rows([0], np.array([True, False]))
        assert not (arr1 == arr2)

    def test_or_operator(self):
        arr1 = LazyFlagArray(3, 2)
        arr1.or_rows([0], np.array([True, False]))
        arr2 = LazyFlagArray(3, 2)
        arr2.or_rows([0], np.array([False, True]))
        result = arr1 | arr2
        assert isinstance(result, np.ndarray)
        np.testing.assert_array_equal(result[0], [True, True])

    def test_ior_lazy(self):
        """Test |= with another LazyFlagArray."""
        arr1 = LazyFlagArray(5, 3)
        arr1.or_rows([1], np.array([True, False, False]))
        arr2 = LazyFlagArray(5, 3)
        arr2.or_rows([1], np.array([False, True, False]))
        arr2.or_rows([3], np.array([False, False, True]))
        arr1 |= arr2
        np.testing.assert_array_equal(arr1[1], [True, True, False])
        np.testing.assert_array_equal(arr1[3], [False, False, True])

    def test_merge_from(self):
        arr1 = LazyFlagArray(5, 3)
        arr1.or_rows([0], np.array([True, False, False]))
        arr2 = LazyFlagArray(5, 3)
        arr2.or_rows([0], np.array([False, True, False]))
        arr2.or_rows([2], np.array([False, False, True]))

        arr1.merge_from(arr2)
        np.testing.assert_array_equal(arr1[0], [True, True, False])
        np.testing.assert_array_equal(arr1[2], [False, False, True])

    def test_or_rows_2d_mismatch_raises(self):
        arr = LazyFlagArray(5, 3)
        with pytest.raises(ValueError, match="mask rows"):
            arr.or_rows([0, 1], np.array([[True, False, False]]))

    def test_memory_efficiency(self):
        """Verify no large allocation when unflagged."""
        import sys

        arr = LazyFlagArray(100_000, 16384)
        # The LazyFlagArray itself should be tiny
        # (no nrows*nchan allocation)
        size = sys.getsizeof(arr) + sys.getsizeof(arr._modified)
        assert size < 1024, f"LazyFlagArray overhead too large: {size} bytes"

    def test_repr(self):
        arr = LazyFlagArray(10, 4)
        r = repr(arr)
        assert "nrows=10" in r
        assert "nchan=4" in r
        assert "zero-backed" in r


class TestLazyFlagContainer:
    """Tests for the LazyFlagContainer class."""

    def test_basic_operations(self):
        container = LazyFlagContainer(3)
        assert len(container) == 3
        assert container[0] is None

    def test_setitem_getitem(self):
        container = LazyFlagContainer(2)
        arr = LazyFlagArray(5, 3)
        container[0] = arr
        assert container[0] is arr

    def test_iter(self):
        container = LazyFlagContainer(3)
        for i in range(3):
            container[i] = LazyFlagArray(5, 3)
        items = list(container)
        assert len(items) == 3

    def test_copy(self):
        container = LazyFlagContainer(2)
        container[0] = LazyFlagArray(5, 3)
        container[0].or_rows([1], np.array([True, False, True]))
        container[1] = LazyFlagArray(5, 3)

        copy = container.copy()
        # Modify original
        container[0].or_rows([1], np.array([False, True, False]))
        # Copy should be independent
        np.testing.assert_array_equal(copy[0][1], [True, False, True])

    def test_ior_lazy_containers(self):
        c1 = LazyFlagContainer(2)
        c1[0] = LazyFlagArray(3, 2)
        c1[0].or_rows([0], np.array([True, False]))
        c1[1] = LazyFlagArray(3, 2)

        c2 = LazyFlagContainer(2)
        c2[0] = LazyFlagArray(3, 2)
        c2[0].or_rows([0], np.array([False, True]))
        c2[1] = LazyFlagArray(3, 2)
        c2[1].or_rows([1], np.array([True, True]))

        c1 |= c2
        np.testing.assert_array_equal(c1[0][0], [True, True])
        np.testing.assert_array_equal(c1[1][1], [True, True])

    def test_container_is_not_none(self):
        """LazyFlagContainer should be truthy (not None)."""
        c = LazyFlagContainer(1)
        assert c is not None
