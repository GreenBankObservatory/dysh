"""
Functions for handling notebook tests.
"""

import random
import string
import nbformat

from nbclient import NotebookClient
from nbformat.v4 import new_code_cell


class NotebookTester(NotebookClient):


    def __init__(self, nb, **kw):
        super().__init__(nb, **kw)


    def get_ref(self, name: str):
        """
        """

        cell = new_code_cell(name)
        self.execute_cell(cell=cell, cell_index=0)
        save_varname = random_varname()
        code = f"""
                import json
                from IPython import get_ipython
                from IPython.display import JSON
        
                {save_varname} = get_ipython().last_execution_result.result
        
                json.dumps({save_varname})
                JSON({{"value" : {save_varname}}})
                """
        output = self.inject(code)
        return output["outputs"][0]["data"]['application/json']['value']


    def inject(self, code: str, pop: bool = False):
        """
        """

        inject_idx = len(self.nb.cells)

        code_cell = new_code_cell(code)
        self.nb.cells.insert(inject_idx, code_cell)

        out = self.execute_cell(cell=code_cell, cell_index=inject_idx)

        if pop:
            self.nb.cells.pop(inject_idx)

        return out


    def get_cell_idx_from_tag(self, tag: str):
        """
        """

        indices = []

        for i,cell in enumerate(self.nb.cells):
            try:
                if tag in cell["metadata"]["tags"]:
                    indices.append(i)
            except KeyError:
                continue

        return indices


def random_varname(length=10):
    """
    Creates a random variable name as string of a given length.
    This is used in testbook to generate temporary variables within the notebook.

    Parameters
    ----------
        length (int)

    Returns:
    --------
        random variable name as string of given length
    """
    return ''.join(random.choice(string.ascii_lowercase) for _ in range(length))
