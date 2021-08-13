""" BioSimulators-compliant command-line interface to the `CBMPy <http://cbmpy.sourceforge.net/>`_ simulation program.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from . import get_simulator_version
from ._version import __version__
from .core import exec_sedml_docs_in_combine_archive
from biosimulators_utils.simulator.cli import build_cli

App = build_cli('biosimulators-cbmpy', __version__,
                'CBMPy', get_simulator_version(), 'http://cbmpy.sourceforge.net/',
                exec_sedml_docs_in_combine_archive)


def main():
    with App() as app:
        app.run()
