#!/usr/bin/python
usage = """wrapper.py [--options] exe seg params"""
description = """a script that cleans up after Omicron processes and moves trigger files into the desired directory structure"""
author = "R. Essick (reed.essick@ligo.org)"

import os
import glob
import subprocess
import time

from optparse import OptionParser

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-C", "--clean-ffl", default=False, action="store_true", help="delete ffl files older than -ffl-shelflife")
parser.add_option("", "--ffl-shelflife", default=600, type="float", help="keep ffl files for at least this long")
parser.add_option("-t", "--trgdir", default="./", type="string", help="the directory from which ffl files will be cleansed")

opts, args = parser.parse_args()

if len(args) != 5:
	raise ValueError("please supply exactly 5 arguments: exe seg params")

#=================================================

print "WARNING: if we move to this solution, we'll need to change the condor sub files, cmd statments in crawler.py!"

subprocess.Popen(args).wait()

#=================================================

### move output files to the correct destination
raise StandardError("WRITE ME")

#=================================================

### clean up ffl files
to = time.time()
if opts.clean_ffl:
	for ffl in sorted(glob.glob("%s/ffconvert.*.ffl"%opts.trgdir)):
		if (to - int(os.path.getmtime( ffl )) > opts.ffl_shelflife):
			os.remove( ffl )


