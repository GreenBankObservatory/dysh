import astropy.units as u
import numpy as np
import pytest
from astroquery.splatalogue import Splatalogue

from dysh.line import SpectralLineSearch


@pytest.mark.requires_internet
class TestSearch:
    """Test dysh.line search functions"""

    def test_local_file(self):
        z = SpectralLineSearch.query_lines(
            cat="gbtlines",
            min_frequency=89.18 * u.GHz,
            max_frequency=89.2 * u.GHz,
            chemical_name="HCO+",
        )
        assert len(z) == 2

    def test_remote_query(self):
        Splatalogue.TIMEOUT = 90  # increase default from 60
        z = SpectralLineSearch.query_lines(
            min_frequency=115 * u.GHz,
            max_frequency=116 * u.GHz,
            chemical_name=" CO ",
            only_NRAO_recommended=True,
        )
        assert len(z) == 1

    def test_keywords(self):
        z = SpectralLineSearch.query_lines(
            min_frequency=1 * u.GHz,
            max_frequency=10 * u.GHz,
            cat="splat",  # also test minimum match for catalogs
            chemical_name="Methyl Formate",
            only_NRAO_recommended=True,
            intensity_lower_limit=20,
            intensity_type="sij",
        )
        assert len(z) == 22
        assert min(z["sijmu2"]) >= 20

    def test_recomb(self):
        z = SpectralLineSearch.recomb(line="Calpha", min_frequency=1 * u.GHz, max_frequency=10 * u.GHz)
        assert len(z) == 100
        # check that requesting only certain column names works
        columns = ["name", "rest_frequency"]
        z = SpectralLineSearch.recomb(
            min_frequency=2 * u.GHz,
            max_frequency=8.4 * u.GHz,
            line="Hbeta",
            only_NRAO_recommended=False,
            columns=columns,
        )
        assert len(z) == 71
        assert z.colnames == columns

        # all recombination lines from a local GBT specific catalog
        z = SpectralLineSearch.recomball(min_frequency=500 * u.MHz, max_frequency=1 * u.GHz, cat="gbtrecomb")
        assert len(z) == 867

    def test_search_with_redshift(self):
        redshift = 0
        z = SpectralLineSearch.recomb(
            line="Halpha", min_frequency=1.022 * u.GHz, max_frequency=8.322 * u.GHz, redshift=redshift
        )
        diff = (z["obs_frequency"] * (1.0 + redshift) - z["rest_frequency"]).data
        redshift = 1.5
        z = SpectralLineSearch.recomb(
            line="Calpha", min_frequency=1 * u.GHz, max_frequency=10 * u.GHz, redshift=redshift
        )
        diff = (z["obs_frequency"] * (1.0 + redshift) - z["rest_frequency"]).data
        assert np.all(np.isclose(diff, 0.0, rtol=1e-8))
        redshift = 5.0
        z = SpectralLineSearch.query_lines(
            cat="gbtlines", min_frequency=1 * u.GHz, max_frequency=10 * u.GHz, redshift=redshift
        )
        diff = (z["obs_frequency"] * (1.0 + redshift) - z["rest_frequency"]).data
        assert np.all(np.isclose(diff, 0.0, rtol=1e-8))
        redshift = 1
        z = SpectralLineSearch.query_lines(
            chemical_name="Carbon Monoxide", min_frequency=100 * u.GHz, max_frequency=120 * u.GHz, redshift=redshift
        )
