#!/usr/bin/python
usage = """wrapper.py [--options] exe seg params"""
description = """a script that cleans up after Omicron processes and moves trigger files into the desired directory structure"""
author = "R. Essick (reed.essick@ligo.org)"

from optparse import OptionParser

#=================================================

parser = OptionParser(usage=usage, description=description)

opts, args = parser.parse_args()

if len(args) != 3:
	raise ValueError("please supply exactly 3 arguments: exe seg params")
cmd = " ".join(args)

#=================================================


"""
run the command
wait
when finished
	remove ffl file
	move trg.xml into requested directory structure
"""

raise StandardError("WRITE ME")

