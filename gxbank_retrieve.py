#!/GWD/appbase/tools/bin/python2.7
from optparse import OptionParser
## deprecated, but argparse requires python >= 2.7 does not work on GSK servers
import logging
#import logging.handlers
import sys
import re
import gzip
import os
import pwd
#from datetime import datetime
import subprocess
#from itertools import izip_longest
import bankfunctions

#define paths
import bankconstants

#extract arguments
parser = OptionParser(description = 'usage: %prog OPTIONS')
parser.add_option('-s', '--sample-list', help = 'Path to the tab delimited list of samples, see FILL IN LATER for guidance'',
                  action = 'store', type = 'string', dest = 'samples', default = '')
parser.add_option('-b', '--bed-file', help = 'Path to the optional BED file specifying hg19 regions to extract, see FILL IN LATER for guidance',
                  action = 'store', type = 'string', dest = 'bed', default = '')
parser.add_option('-o', '--vcf-path', help = 'Path of directory where to write the gzip\'d VCF file [default: working directory]',
                  action = 'store', type = 'string', dest = 'vcfpath', default = os.getcwd())
parser.add_option('-v', '--vcf', help = 'Name of the output VCF',
                  action = 'store', type = 'string', dest = 'vcf', default = '')
parser.add_option('-d', '--data-bank', help = 'Path to data bank [default: Production]',
                  action = 'store', type = 'string', dest = 'bank', default = bankconstants.prodbank)
parser.add_option('-f', '--force', help = 'Force overwrite if output VCF already exists',
                  action = 'store_true', dest = 'force', default = False)
"""
#Will add this option (and a consensus by coordinate option) later
parser.add_option('-c', '--skip-consensus', help = 'Do not report consensus genotypes by subject [default: consensus]',
                  action = 'store_false', dest = 'consensus', default = True)
"""
(options, args) = parser.parse_args()
#Get absolute paths for inclusion in log / error messages
out_path = os.path.abspath(options.vcfpath)
bank_path = os.path.abspath(options.bank)
sample_path = os.path.abspath(options.samples)
bed_path = os.path.abspath(options.bed)

#validate out vcf file name: allow alphanum, dashes, and underscores
if options.vcf == '' :
	raise Exception('ERROR: vcf file name cannot be blank')
elif re.search(r"^([\w\-])+$", options.vcf) == None :
	raise Exception('ERROR: vcf file name may only contain alphanumeric characters, dash or underscore')
#validate outpath - must be valid path to a writable directory
if not os.access(out_path, os.W_OK) :
	raise Exception('ERROR: cannot write to output directory [' + out_path + ']')
	
#check if vcf already exists and not forcing
if os.path.isfile(out_path + '.vcf.gz') and not options.force :
	raise Exception('ERROR: vcf [' + out_path + '.vcf.gz] already exists')

#check sample list is readable and validate
errsampleheader = "ERROR: Sample list [{0}] missing {1} column"
errsamplebatch = "ERROR: Sample list[{0}] invalid, batch [{1}] in platform [{2}] is not in the bank [{3}]"
sample_header = {}
expected_fields = ['USUBJID','LabID','Platform','Batch']
if not os.access(sample_path, os.R_OK) :
	raise Exception ('ERROR: Could not find or read sample list [' + sample_path + ']')
	with open(sample_path,'r') as samples :
		#validate header
		header = bankfunctions.read_through_headers(samples,'#')
		fields = header.strip().split('\t')
		for i,field in enumerate(fields) :
			sample_header[field] = i
		for field in expected_fields :
			if field not in sample_header.keys() :
				raise Exception(errsampleheader.format(sample_path,field))
		#read sample list
		sample_by_platform_batch = {}
		for record in samples :
			if record[0] != '#' :
				fields = record.strip().split('\t')
				platform = fields[sample_header['Platform']]
				batch = fields[sample_header['Batch']]
				labid = fields[sample_header['LabID']]
				if platform not in sample_by_platform_batch.keys() :
					sample_by_platform_batch[platform] = {}
				if batch not in sample_by_platform_batch[platform].keys :
					sample_by_platform_batch[platform][batch] = set()
				sample_by_platform_batch[platform][batch].add(labid)
		#validate platform & batch
		for platform in batch_by_platform.keys() :
			for batch in batch_by_platform[platform].keys() :
				#just checking for presence of vcf, will assume tbi is also present
				batch_path = '/'.join[bank_path,'Data',platform,batch] + 'vcf.gz'
				if not os.access(batch_path, os.R_OK) :
					Raise Exception(errsamplebatch.format(sample_path,platform,batch,bank_path))


#if bed file specified, check readable & valid
errbedchr = "ERROR: BED file [{0}] must use hg19 chromosomes like chr1 [{1}] on line {2} is invalid"
errbedpos = "ERROR: BED file [{0}] position must only contain digits [{1}] on line {2} is invalid"
errbedrange = "ERROR: BED file [{0}] must have end position >= start. End [{1}] on line {2} is < start [{3}]"
if options.bed != '' :
	if not os.access(bed_path, os.R_OK) :
		raise Exception ('ERROR: Could not find or read bed file [' + bed_path + ']')
	#looked at bedops as validator but requires end > start which is not necessarily required for vcftools
	with open(bed_path,'r') as bed :
		for linenum, line in enumerate(bed) :
			if line[0] != '#' :
				fields = line.strip().split('\t')
				if (fields[0][0:3] != 'chr') :
					raise Exception (errbedchr.format(bed_path,fields[0],linenum + 1))
				if re.match(r"^\d+$",fields[1]) == None :
					raise Exception (errbedpos.format(bed_path,fields[1],linenum + 1))
				if re.match(r"^\d+$",fields[2]) == None :
					raise Exception (errbedpos.format(bed_path,fields[2],linenum + 1))
				if int(fields[2]) < int(fields[1]) :
					raise Exception (errbedrange.format(bed_path,fields[2],linenum + 1,fields[1]))

#create rolling log
RetrieveLogger = logging.getLogger()
RetrieveLogger.setLevel(logging.DEBUG)
#Doesn't seem to be a way to set backupCount to infinite, tried 0 and -1, both effectively mean no rolling, just single file that gets overwritten when full
#Setting to unreasonably high # to account for this
RollingLog = logging.handlers.RotatingFileHandler(options.bank + '/Logs/Retrieve.log', maxBytes=50000000, backupCount=5000)
#will this die if can't write to log?
LogRecordFormat = logging.Formatter(fmt='%(process)d\t%(levelname)s\t%(asctime)s\t%(message)s')
RollingLog.setFormatter(LogRecordFormat)
RetrieveLogger.addHandler(RollingLog)

#start log with user and arguments - will validate file contents as process
if options.force :
	logging.info('[%s;%s] is retrieving samples from [%s] and writing to [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], sample_path, out_path + '.vcf.gz')
else :
	logging.info('[%s;%s] is retrieving samples from [%s] and overwriting [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], sample_path, out_path + '.vcf.gz')
if options.bed != '' :
	logging.info('Restricting data to regions specified in [%s]', bed_path)
logging.info('Using data in bank [%s]', bank_path)
	
for platform in sample_by_platform_batch.keys() :
	for batch in sample_by_platform_batch[platform].keys() :
		with open(out_path + '/list.samples.list','w') as samples :
			for sample in sample_by_platform_batch[platform][batch] :
				samples.write(sample + '\n')

