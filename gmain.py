#!/usr/bin/env python

import sys
from loom.__main__ import run_with_sys_argv

argv = sys.argv[1:]
argv += ['--gui-mode', 'True']

run_with_sys_argv(argv)
