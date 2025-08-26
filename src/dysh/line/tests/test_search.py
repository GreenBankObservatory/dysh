import astropy.units as u

from dysh.line.search import SpectralLineSearch


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
            chemical_name="Methyl Formate",
            only_NRAO_recommended=True,
            intensity_lower_limit=20,
            intensity_type="sij",
        )
        assert len(z) == 22
        assert min(z["sijmu2"]) >= 20
