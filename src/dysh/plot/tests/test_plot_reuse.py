"""Test plot reuse and cleanup"""

import gc

import matplotlib.pyplot as plt

from dysh.spectra import Spectrum


class TestPlotReuse:
    """Tests for proper cleanup when reusing plot objects."""

    def setup_method(self):
        """Disable interactive plotting for tests."""
        plt.ioff()

    def teardown_method(self):
        """Close all figures after each test."""
        plt.close("all")

    def test_figure_cleanup_on_replot(self):
        """
        Test that figures don't accumulate in pyplot registry when calling plot() multiple times.

        This test reproduces Issue #813 where repeatedly calling plot() on the same Spectrum
        object would accumulate figures in pyplot's global registry, eventually causing
        the interactive plotter to hang.
        """
        # Disable GC to force accumulation (makes test more reliable)
        gc.disable()
        try:
            s = Spectrum.fake_spectrum()

            # Track figure and selector IDs to verify they change
            figure_ids = []
            selector_ids = []

            # Create plot 3 times (simulates user repeatedly plotting same spectrum)
            for i in range(3):
                plotter = s.plot(select=True)  # Enable selector to test event handler cleanup

                # Record IDs at this moment
                figure_ids.append(id(plotter._figure))
                if plotter._selector:
                    selector_ids.append(id(plotter._selector))

                # Critical assertion: pyplot should only track 1 figure, not accumulate
                fig_nums = plt.get_fignums()
                assert len(fig_nums) == 1, (
                    f"Iteration {i + 1}: Expected 1 figure in pyplot registry, "
                    f"got {len(fig_nums)}. This indicates figures are accumulating "
                    f"without cleanup (Issue #813)."
                )

            # Verify that new Figure objects were created each time (not reused broken objects)
            assert len(set(figure_ids)) == 3, (
                "Expected 3 different Figure objects across 3 plot() calls, "
                f"but got {len(set(figure_ids))} unique IDs. "
                "This indicates Figure objects are being improperly reused."
            )

            # Verify that new Selector objects were created each time
            if selector_ids:
                assert len(set(selector_ids)) == 3, (
                    "Expected 3 different Selector objects across 3 plot() calls, "
                    f"but got {len(set(selector_ids))} unique IDs. "
                    "This indicates Selector objects are being improperly reused "
                    "without proper cleanup of event handlers."
                )

        finally:
            gc.enable()

    def test_selector_event_handler_cleanup(self):
        """
        Test that InteractiveSpanSelector event handlers are properly disconnected.

        This verifies that the disconnect() method properly cleans up all 4 event handlers
        to prevent them from accumulating or referencing closed canvases.
        """
        s = Spectrum.fake_spectrum()

        # Create initial plot with selector
        plotter = s.plot(select=True)
        selector1 = plotter._selector

        # Manually disconnect the selector
        selector1.disconnect()

        # Verify connection IDs were cleared
        assert selector1.cid_press is None, "cid_press should be None after disconnect()"
        assert selector1.cid_release is None, "cid_release should be None after disconnect()"
        assert selector1.cid_motion is None, "cid_motion should be None after disconnect()"
        assert selector1.cid_key is None, "cid_key should be None after disconnect()"

    def test_no_freeze_after_multiple_plots(self):
        """
        Test that creating multiple plots doesn't cause pyplot state corruption.

        This is a regression test for the freeze that occurred in Strategy 1 testing,
        where creating multiple plots would cause TkAgg to freeze when closing a figure.
        """
        gc.disable()
        try:
            s = Spectrum.fake_spectrum()
            plotters = []

            # Create 5 plots (strategy 1 used 10, but 5 is sufficient for test)
            for _i in range(5):
                plotter = s.plot(select=True)
                plotters.append(plotter)

            # After all plots, only 1 should exist in pyplot registry
            fig_nums = plt.get_fignums()
            assert len(fig_nums) == 1, (
                f"After creating 5 plots, expected 1 figure in pyplot registry, "
                f"got {len(fig_nums)}. Figures are accumulating without cleanup."
            )

            # Simulate closing the figure (in real scenario, user clicks X)
            # This should NOT cause a freeze or error
            plt.close(plotters[-1]._figure)

            # After closing, pyplot should have 0 figures
            assert len(plt.get_fignums()) == 0, "Figure should be removed from pyplot registry after close"

        finally:
            gc.enable()

    def test_plotter_object_reuse(self):
        """
        Test that the SpectrumPlot object itself is reused across plot() calls.

        This verifies that s.plot() returns the same SpectrumPlot object each time,
        but properly cleans up and recreates internal Figure/Selector objects.
        """
        s = Spectrum.fake_spectrum()

        plotter1 = s.plot()
        plotter2 = s.plot()
        plotter3 = s.plot()

        # Same plotter object should be reused
        assert plotter1 is plotter2 is plotter3, "SpectrumPlot object should be reused across plot() calls"

        # Note: at this point, all three variables point to same plotter with same ._figure,
        # so we can't compare them directly. This test mainly verifies no crash occurs.

        # Final state: only 1 figure should exist
        assert len(plt.get_fignums()) == 1, "Only final figure should remain in pyplot registry"
