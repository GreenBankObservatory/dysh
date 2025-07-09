#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 15:39:39 2025

@author: mpound
"""

import sys

from astropy.table import Table

for file in sys.argv[1:]:
    print(f"========== {file} ============")
    t = Table.read(file, format="ascii.ecsv")
    t.pprint_all()
