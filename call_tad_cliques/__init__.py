# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from __future__ import absolute_import

import call_tad_cliques.io
import call_tad_cliques.cis
import call_tad_cliques.trans
import call_tad_cliques.run_nchg
import call_tad_cliques.call_cliques
import call_tad_cliques.preprocessing

from importlib.metadata import version
import re

__version__ = version(re.match(r"^[^.]+", __name__)[0])