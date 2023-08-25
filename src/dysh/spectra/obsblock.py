class Obsblock:
    """Class that holds a series of spectra on which bulk operations can be performed"""

    def __init__(self, speclist, index):
        self._speclist = speclist
        self._index = index  # pandas dataframe

    def __getitem__(self, i):
        return self._speclist[i]

    def __len__(self):
        return len(self._speclist)

    def __op__(self, opname):
        pass

    def baseline(self, order, exclude=None, **kwargs):
        """compute and optionally remove a baseline"""
        pass
