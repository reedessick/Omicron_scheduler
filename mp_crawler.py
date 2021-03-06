#!/usr/bin/python
usage = "crawler.py [--options] config.ini"
description = "a script that lives persistently and forks off Omicron jobs periodically"
author = "R. Essick (reed.essick@ligo.org)"

#import lal
from lal import gpstime
from pylal import Fr

import os
import sys ### only needed if we touch sys.stdout?
import multiprocessing as mp
import subprocess as sp
import numpy as np
import time
import logging
from ConfigParser import SafeConfigParser

from optparse import OptionParser

#=================================================
def __safe_fork(cmd, stdin=None, stdout=None, stderr=None):
	"""
	helper for a double fork. This is called by multiprocessing to orphan the actual processes
	"""
	if stdin:
		stdin_obj = open(stdin, "r")
	if stdout:
		stdout_obj = open(stdout, "w")
	if stderr:
		stderr_obj = open(stderr, "w")

	if stdin and stdout and stderr:
		p = sp.Popen(cmd, stdin=stdin_obj, stdout=stdout_obj, stderr=stderr_obj) ### launch subprocess
		### we don't wait, communicate or call
		### by returning, this process will die and orphan the subprocess
	elif stdin and stdout:
		p = sp.Popen(cmd, stdin=stdin_obj, stdout=stdout_obj)
	elif stdin and stderr:
		p = sp.Popen(cmd, stdin=stdin_obj, stderr=stderr_obj)
	elif stdout and stderr:
		p = sp.Popen(cmd, stdout=stdout_obj, stderr=stderr_obj)
	elif stdin:
		p = sp.Popen(cmd, stdin=stdin_obj)
	elif stdout:
		p = sp.Popen(cmd, stdout=stdout_obj)
	elif stderr:
		p = sp.Popen(cmd, stderr=stderr_obj)
	else:
		p = sp.Popen(cmd)

	if stdin:
		stdin_obj.close()
	if stdout:
		stdout_obj.close()
	if stderr:
		stderr_obj.close()

###
def safe_fork(cmd, stdin=None, stdout=None, stderr=None):
	"""
	a wrapper for a double fork
	"""
	p = mp.Process(target=__safe_fork, args=(cmd, stdin, stdout, stderr)) ### define job
	p.start() ### launch
        p.join() ### call, will orphan subprocess when it finishes

###
def report(statement, verbose):
	"""
	wrapper for reporting output
	"""
	if verbose:
		print statement
	logger.info(statement)

###	
def find_frames(ldr_server, ldr_url_type, ldr_type, ifo, start, stride, verbose=False):
	"""
	wrapper for ligo_data_find
	"""
	end = start+stride

        cmd = "ligo_data_find --server=%s --url-type=%s --type=%s --observatory=%s --gps-start-time=%d --gps-end-time=%d"%(ldr_server, ldr_url_type, ldr_type, ifo, start, end)
        report(cmd, verbose)

	p = sp.Popen(cmd.split(), stdout=sp.PIPE, stderr=sp.STDOUT)
	
	frames = p.communicate()[0].replace("No files found!", "").replace("\n", " ") ### handle an empty list appropriately

	frames = frames.replace("file://localhost","")

	return [l for l in frames.split() if l.endswith(".gwf")]

###
def coverage(frames, start, stride):
	"""
	determines the how much of [start, start+stride] is covered by these frames

	assumes non-overlapping frames!
	"""
	### generate segments from frame names
	segs = [[float(l) for l in frame.strip(".gwf").split("-")[-2:]] for frame in sorted(frames)]

	### check whether segments overlap with desired time range
	covered = 1.0*stride

	end = start + stride
	for s, d in segs:
		e = s+d

		if (s < end) and (start < e): ### at least some overlap
			covered -= min(e, end) - max(s, start) ### subtract the overlap

		if covered <= 0:
			break

	return 1 - covered/stride ### return fraction of coverage

###
def str_framecache(frames, ifo, type):
	"""
	build a string for the framecache
	"""
	S = ""
	for frame in frames:
		s, d = frame.strip(".gwf").split("-")[-2:]
		S += "%s %s %s %s %s\n"%(ifo, type, s, d, frame)
	return S

###
def extract_scisegs(frames, channel, bitmask, start, stride):
	"""
	extract scisegs from channel in frames using bitmask
	"""
	if not frames: ### empty list, so no segments
		return []

	### extract vectors and build segments
	segset = []
	for frame in frames:
		### extract the vector from the frame
		vect, s, ds, dt, xunit, yunit = Fr.frgetvect1d(frame, channel)		
		n = len(vect)

		### build time vector        add starting time
		t = np.arange(0, dt*n, dt) + s+ds

		### determine whether state acceptable
		### add "False" buffers to set up the computation of start and end time
		state = np.concatenate( ([False], vect == bitmask, [False]))

		### determine beginning of segments
		###      i=False      i+1 = True  strip the trailing buffer
		b = ( (1-state[:-1])*(state[1:]) )[:-1].astype(bool)
		b = t[b] ### select out times

		### determine end of segments
                ###     i=True     i+1=False      strip the leading buffer
		e = ( (state[:-1])*(1-state[1:]) )[1:].astype(bool) 
		e = t[e] + dt ### select out times
		              ### extra dt moves these markers to the end of segments

		### stitch together start and end times, append to global list
		segset += list( np.transpose( np.array( [b, e] ) ) )

	if not segset: ### empty list
		return []

	### clean up segs!
	segs = []
	seg1 = segset[0]
	for seg2 in segset[1:]:
		if seg1[1] == seg2[0]:
			seg1[1] = seg2[1] ### join the segments
		else:
			segs.append( list(seg1) )
			seg1 = seg2
	segs.append( list(seg1) )

	### return final list of lists!
	return segs

###
def str_omicron_config(framecache, channels, samplefrequency=4096, chunkduration=32, blockduration=32, overlapduration=4, windows=[2,4], fftplan="ESTIMATE", frequencyrange=[32,2048], qrange=[3,141], mismatch=0.2, snrthreshold=5.5, nmax=1e6, clustering="time", outputdir="./", format=["xml"], verbosity=0, writepsd=0, writetimeseries=0, writewhiteneddata=0, plotstyle="GWOLLUM"):
	"""
	builds the string that represents the omicron parameter file
	WARNING: may be sub-optimal if required extremely repetitively because of the way we concatenate strings
		(strings are immutable in python, so we create many objects)
	"""
	s = ""

	### data 
	s += "DATA\tLCF\t%s\n"%framecache
	for item in channels.items():
		s += "DATA\tCHANNELS\t%s\nDATA\tNATIVEFREQUENCY\t%d\n"%item

	s += "DATA\tSAMPLEFREQUENCY\t%d\n"%samplefrequency

	### parameters
	s += "PARAMETER\tCHUNKDURATION\t%d\n"%chunkduration
	s += "PARAMETER\tBLOCKDURATION\t%d\n"%blockduration
	s += "PARAMETER\tOVERLAPDURATION\t%d\n"%overlapduration
	s += "PARAMETER\tWINDOWS\t%s\n"%" ".join([str(w) for w in windows])
	s += "PARAMETER\tFFTPLAN\t%s\n"%fftplan
	s += "PARAMETER\tFREQUENCYRANGE\t%.4f\t%.4f\n"%tuple(frequencyrange)
	s += "PARAMETER\tQRANGE\t%.4f\t%.4f\n"%tuple(qrange)
	s += "PARAMETER\tMISMATCHMAX\t%.4f\n"%mismatch

	### triggers
	s += "TRIGGER\tSNRTHRESHOLD\t%.4f\n"%snrthreshold
	s += "TRIGGER\tNMAX\t%d\n"%int(nmax)
	s += "TRIGGER\tCLUSTERING\t%s\n"%clustering

	### output
	s += "OUTPUT\tDIRECTORY\t%s\n"%outputdir
	s += "OUTPUT\tFORMAT\t%s\n"%(",".join(format))
	s += "OUTPUT\tVERBOSITY\t%d\n"%verbosity
	s += "OUTPUT\tWRITEPSD\t%d\n"%writepsd
	s += "OUTPUT\tWRITETIMESERIES\t%d\n"%writetimeseries
	s += "OUTPUT\tWRITEWHITENEDDATA\t%d\n"%writewhiteneddata
	s += "OUTPUT\tPLOTSTYLE\t%s\n"%plotstyle

	return s

###
def str_omicron_sub(universe, executable, arguments, log, output, error, getenv=True, notification="never"):
	"""
	builds a string that represents a condor sub file for omicron
	"""
	s = ""

	s += "universe     = %s\n"%universe
	s += "executable   = %s\n"%executable	
	s += "arguments    = %s\n"%" ".join(arguments)
#	s += "arguments    = %s\n"%" ".join(["$(%s)"%a for a in arguments])
	s += "log          = %s\n"%log
	s += "output       = %s\n"%output
	s += "error        = %s\n"%error
	s += "notification = %s\n"%notification
	s += "queue 1\n"

	return s

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--no-robot-cert", default=False, action="store_true", help="do not use robot cert specified in config file")

parser.add_option("-s", "--gps-start", default=None, type="int")
parser.add_option("-e", "--gps-end", default=np.infty, type="float")

opts, args = parser.parse_args()

if len(args) != 1:
	raise StandardError("Please supply only a single argument")
configfile = args[0]

#=================================================

config = SafeConfigParser()
config.read(configfile)

#=================================================
### setup logger to record processes
logfilename = config.get("general","logfile")

### ensure that path to log will exist
logpath = "/".join(logfilename.split("/")[:-1])
if logfilename[0] == "/":
	logpath = "/%s"%logpath
if not os.path.exists(logpath):
	os.makedirs(logpath)

global logger
logger = logging.getLogger('crawler_log')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(message)s')

### redirect stderr into logger
hdlr = logging.FileHandler(logfilename)
hdlr.setFormatter(formatter)
hdlr.setLevel(logging.INFO)
logger.addHandler(hdlr)

#=================================================
### source environment scripts
report("sourcing environment", opts.verbose)

for reason, script in config.items("environment"):
	cmd = "source %s"%(script)
	report(cmd, opts.verbose)

	p = sp.Popen(cmd.split(), executable="/bin/bash", stdout=sp.PIPE)

	report(p.communicate()[0].strip("\n"), opts.verbose)
	if 0 != p.returncode:
		raise StandardError("failed to source %s"%script)

### robot cert
if not opts.no_robot_cert:
	if os.environ.has_key("X509_USER_PROXY"):
		del os.environ['X509_USER_PROXY']

	### set cert and key
	os.environ['X509_USER_CERT'] = config.get('robot cert', 'robot_certificate')
	os.environ['X509_USER_KEY'] = config.get('robot cert', 'robot_key')


#=================================================
report("pulling out parameters from : %s"%configfile, opts.verbose)

### pull out basic parameters
ifo = config.get("general", "ifo")
outputdir = config.get("general", "outputdir")
stride = config.getint("general", "stride")
delay = config.getint("general", "delay")
padding = config.getint("general", "padding")

max_wait = config.getint("general", "max_wait")

#========================

### pull out ldr params
ldr_server = config.get("ligo_data_find", "server")
ldr_url_type = config.get("ligo_data_find", "url-type")
ldr_type = config.get("ligo_data_find", "type")

#========================

### pull out sciseg params
sciseg_channel = config.get("scisegs","channel")
sciseg_bitmask = config.getint("scisegs","bitmask")


#========================

### pull out omicron run parameters
block = config.getboolean("omicron", "block") ### whether to block
condor = config.getboolean("omicron", "condor") ### whether to use condor
scisegs = config.getboolean("omicron","scisegs") ### whether to use scisegs
executable = config.get("omicron", "executable")

### output formatting
format = eval(config.get("omicron","format"))
verbosity = config.getint("omicron","verbosity")
writepsd = config.getint("omicron","writepsd")
writetimeseries = config.getint("omicron","writetimeseries")
writewhiteneddata = config.getint("omicron","writewhiteneddata")
plotstyle = config.get("omicron","plotstyle")

### set up params_string for each channel_set
params_strings = []
for section_name in sorted(config.options("channel sets")):
	report("setting up template omicron params file for : %s"%section_name, opts.verbose)

	params_strings.append( ( section_name, 
	                         str_omicron_config( "%s", # will be filled in later
	                                             dict( ("%s1:%s"%(ifo,key.upper()), float(value)) for key, value in config.items("%s channels"%section_name) ),
	                                             samplefrequency = config.getint(section_name,"samplefrequency"),
	                                             chunkduration = config.getint(section_name,"chunkduration"),
	                                             blockduration = config.getint(section_name,"blockduration"),
	                                             overlapduration = config.getint(section_name,"overlapduration"),
	                                             windows = eval(config.get(section_name,"windows")),
	                                             fftplan = config.get(section_name,"fftplan"),
	                                             frequencyrange = eval(config.get(section_name,"frequencyrange")),
	                                             qrange = eval(config.get(section_name,"qrange")),
	                                             mismatch = config.getfloat(section_name,"mismatch"),
	                                             snrthreshold = config.getfloat(section_name,"snrthreshold"),
	                                             nmax = config.getint(section_name,"nmax"),
	                                             clustering = config.get(section_name,"clustering"),
	                                             outputdir="%s",  # will be filled in later
	                                             format=format, 
	                                             verbosity=verbosity, 
	                                             writepsd=writepsd,
	                                             writetimeseries=writetimeseries,
	                                             writewhiteneddata=writewhiteneddata,
	                                             plotstyle=plotstyle
	                                           )
	                       )
	                     )

### set up condor files if needed
if condor:
	report("setting up condor sub files", opts.verbose)

	### write sub template
	report("building sub template", opts.verbose)
	sub_string = str_omicron_sub("vanilla", executable, ["%s", "%s"], "%s", "%s", "%s", getenv=True, notification="never")

#=================================================

### setting up initial time
report("", opts.verbose)
if opts.gps_start == None:
	t = ( int(gpstime.gps_time_now()) / stride)*stride
else:
	t = (opts.gps_start/stride)*stride ### round to integer number of strides

#=================================================
# LOOP until we "finish"
#=================================================
while t < opts.gps_end:

	report("=========================================================================", opts.verbose)
	report("processing stride: [%d-%d, %d+%d]"%(t, padding, t+stride, padding), opts.verbose)

	### wait to analyze this stride
	nowgps = float( gpstime.gps_time_now() )
	wait = (t+stride+padding) + delay - nowgps 
	if wait > 0:
		report("sleeping for %d sec"%wait, opts.verbose)
		time.sleep(wait)

	### build directories
	t5 = t/100000
	segdir = "%s/segments/%s-%d/"%(outputdir, ifo, t5)
	logdir = "%s/logs/%s-%d/"%(outputdir, ifo, t5)
	framedir = "%s/frames/%s-%d/"%(outputdir, ifo, t5)
	trgdir = "%s/triggers/%s-%d/"%(outputdir, ifo, t5)
	for directory in [outputdir, segdir, logdir, framedir, trgdir]:
        	if not os.path.exists(directory):
			report("building directory : %s"%directory, opts.verbose)
                	os.makedirs(directory)

	if condor:
		condordir = "%s/condor/%s-%d/"%(outputdir, ifo, t5)
		if not os.path.exists(condordir):
			report ("building directory : %s"%condordir, opts.verbose)
			os.makedirs(condordir)
	
	### find frames within time window
	report("finding frames within stride", opts.verbose)
	frames = find_frames(ldr_server, ldr_url_type, ldr_type, ifo, t-padding, stride+2*padding, verbose=opts.verbose) 
	covered = coverage( frames, t-padding, stride+2*padding) ### find out the coverage

	### keep looking every second until we either find frames or time out
	if covered < 1.0:
		report("coverage = %.5f < 1.0, we'll check every second for more frames and wait at most %d seconds before proceeding."%(covered, max_wait), opts.verbose)

	while (covered < 1.0) and ( (float(gpstime.gps_time_now()) - ( (t+stride+padding) + delay ) ) < max_wait ):

		###
		time.sleep( 1 ) # don't break the file system
		###

		frames = find_frames(ldr_server, ldr_url_type, ldr_type, ifo, t-padding, stride+2*padding, verbose=False) ### don't report this every time in the loop
		covered = coverage( frames, t-padding, stride+2*padding) ### find out the coverage

		if covered >= 1.0:
			report("covered >= 1.0")

	if covered < 1.0:
		report("coverage = %.5f < 1.0, but we've timed out after waiting at least %d seconds."%(covered, max_wait), opts.verbose)

	### write framecache
	framecache = "%s/%s_%d-%d.lcf"%(framedir, ifo, t, stride)
	report("writing framecache : %s"%framecache, opts.verbose)

	framecache_obj = open(framecache, "w")
	framecache_obj.write( str_framecache(frames, ifo, ldr_type) )
	framecache_obj.close()

	### if we have data, process it!
	if not frames:
		report("no frames found! skipping...", opts.verbose)

	else:
		### find scisegs
		segfile = "%s/%s_%d-%d.seg"%(segdir, ifo, t, stride)
		if scisegs: ### extract from ODC vector in frames
			report("extracting scisegs to : %s"%(segfile), opts.verbose)
			segs = extract_scisegs(frames, "%s1:%s"%(ifo, sciseg_channel), sciseg_bitmask, t-padding, t+stride+padding)

		else: ### use entire segment for analysis
			report("using analysis segment as scisegs", opts.verbose)
			segs = [(t-padding, t+stride+padding)]

		report("writing scisegs : %s"%segfile, opts.verbose)
		file_obj = open(segfile, "w")
		for a, b in segs:
			file_obj.write("%d %d"%(a, b))
		file_obj.close()

		### launch jobs

		if block: ### run from here and block
			procs = []
			for channel_set, params_string in params_strings: ### iterate over separate jobs
				### write omicron params file
				params = "%s/%s_%s-%d-%d.params"%(logdir, ifo, channel_set, t, stride)
				report("writing params : %s"%params, opts.verbose)

				file_obj = open(params, 'w')
				file_obj.write(params_string%(framecache, trgdir))
				file_obj.close()

				out = "%s/%s_%s-%d-%d.out"%(logdir, ifo, channel_set, t, stride)
				err = "%s/%s_%s-%d-%d.err"%(logdir, ifo, channel_set, t, stride)
				cmd = "%s %s %s"%(executable, segfile, params)
				report(cmd, opts.verbose)
				report("out: %s"%out, opts.verbose)
				report("err: %s"%err, opts.verbose)

				out_obj = open(out, "w")
				err_obj = open(err, "w")
				procs.append( (cmd, sp.Popen(cmd.split(), stdout=out_obj, stderr=err_obj)) )
				out_obj.close()
				err_obj.close()

			while procs:
				cmd, p = procs.pop(0)
				p.wait() ### block!
			
		elif condor: ### run under condor
			procs = []
			for channel_set, params_string in params_strings: ### iterate over separate jobs
				### write omicron params file
				params = "%s/%s_%s-%d-%d.params"%(logdir, ifo, channel_set, t, stride)
				report("writing params : %s"%params, opts.verbose)

				file_obj = open(params, 'w')
				file_obj.write(params_string%(framecache, trgdir))
				file_obj.close()

				log = "%s/%s_%s-%d-%d.log"%(logdir, ifo, channel_set, t, stride)
				out = "%s/%s_%s-%d-%d.out"%(logdir, ifo, channel_set, t, stride)
				err = "%s/%s_%s-%d-%d.err"%(logdir, ifo, channel_set, t, stride)
			
				sub = "%s/%s_%s-%d-%d.sub"%(condordir, ifo, channel_set, t, stride)

				report("sub: %s"%sub, opts.verbose)
				report("log: %s"%log, opts.verbose)
				report("out: %s"%out, opts.verbose)
				report("err: %s"%err, opts.verbose)

				### write sub
				file_obj = open(sub, "w")
				file_obj.write(sub_string%(segfile, params, log, out, err))
				file_obj.close()

				### submit through condor
				cmd = "condor_submit %s"%sub
				report(cmd, opts.verbose)

				procs.append( (cmd, sp.Popen(cmd.split(), stdout=sys.stdout, stderr=sys.stderr)) )

			while procs:
				cmd, p = procs.pop(0)
				p.wait()

		else: ### run from here but do not block		
			for channel_set, params_string in params_strings: ### iterate over separate jobs
				### write omicron params file
				params = "%s/%s_%s-%d_%d-%d.params"%(logdir, ifo, channel_set, t, stride)
				report("writing params : %s"%params, opts.verbose)

				file_obj = open(params, 'w')
				file_obj.write(params_string%(framecache, trgdir))
				file_obj.close()

				out = "%s/%s_%s-%d-%d.out"%(logdir, ifo, channel_set, t, stride)
				err = "%s/%s_%s-%d-%d.err"%(logdir, ifo, channel_set, t, stride)
				cmd = "%s %s %s"%(executable, segfile, params)
				report(cmd, opts.verbose)
				report("out: %s"%out, opts.verbose)
				report("err: %s"%err, opts.verbose)

				safe_fork(cmd.split(), stdout=out, stderr=err)

	report("Done with stride: [%d-%d, %d+%d]"%(t, padding, t+stride, padding), opts.verbose)

	### increment!
	t += stride
