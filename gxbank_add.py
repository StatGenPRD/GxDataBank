#!/GWD/appbase/tools/bin/python2.7
from optparse import OptionParser
## deprecated, but argparse requires python >= 2.7 does not work on GSK servers
import logging
import logging.handlers
import sys
import re
import gzip
import os
import pwd
from datetime import datetime
import subprocess
from itertools import izip_longest

#define paths
import bankconstants

#extract arguments
parser = OptionParser(description = 'usage: %prog OPTIONS')
parser.add_option('-b', '--batch-name', help = 'Name of batch of data as supplied by the lab and recorded in Reflection',
                  action = 'store', type = 'string', dest = 'batch', default = '')
parser.add_option('-p', '--platform-name', help = 'Name that uniquely identifies the platform for the data',
                  action = 'store', type = 'string', dest = 'platform', default = '')
parser.add_option('-v', '--vcf', help = 'Path and name of the gzip\'d VCF file to add',
                  action = 'store', type = 'string', dest = 'vcf', default = '')
parser.add_option('-d', '--data-bank', help = 'Path to data bank [default: Production]',
                  action = 'store', type = 'string', dest = 'bank', default = bankconstants.prodbank)
parser.add_option('-f', '--force', help = 'Force overwrite if VCF already exists in bank',
                  action = 'store_true', dest = 'force', default = False)
(options, args) = parser.parse_args()


#validate batch name: allow alphanum, dashes, and underscores
if options.batch == '' :
	raise Exception('ERROR: batch-name cannot be blank')
elif re.search(r"^([\w\-])+$", options.batch) == None :
	raise Exception('ERROR: batch-name may only contain alphanumeric characters, dash or underscore')

#validate platform name - must match an existing platform definition
platform_path = options.bank + '/PlatformDefinitions/' + options.platform + '.vcf.gz'
if not os.access(platform_path, os.R_OK) :
	raise Exception ('ERROR: Could not find or read platform definition [' + platform_path + ']')

#validate vcf path - must be valid path to a file that can be read
if not os.access(options.vcf, os.R_OK) :
	raise Exception('ERROR: Could not find or read vcf [' + options.vcf + ']')

#validate data bank path
#must be valid path to a writable directory
if not os.access(options.bank + '/Data', os.W_OK) :
	raise Exception('ERROR: user cannot write to bank [' + options.bank + '/Data]')
	
#check if batch name already exists in bank and not forcing
if os.path.isfile(options.bank + '/Data/' + options.batch + '.bcf') and not options.force :
	raise Exception('ERROR: batch-name [' + options.batch + '] already loaded to this bank [' + options.bank + '/Data]')
	
#must pass vcf-validator
validator_result = subprocess.check_output(["perl", "-I", bankconstants.vcftools_base + '/perl', bankconstants.vcf_validator, options.vcf],stderr=subprocess.STDOUT)
validator_result = validator_result.strip().split('\n')
for error in validator_result :
	print error
	if re.search(r"Not required", error) == None and error != '':
		raise Exception('ERROR: vcf not valid [' + error + ']')

#must match platform definition - no extra/missing records
#CHROM,POS,ID,REF,ALT must match
#Same sort order - may be able to relax this requirement later
err1 = "ERROR: VCF file [{0}] does not match PLATFORM file [{1}] on {2} at variant record {3}. VCF={4}; PLATFORM={5}"
err2 = "ERROR: VCF file [{0}] does not match PLATFORM file [{1}]. Different line counts."
fields = ["CHROM","POS","ID","REF","ALT"]
with gzip.open(options.vcf, 'rb') as vcf, gzip.open(platform_path, 'rb') as platform :
	#remove headers from each file
	while vcf.readline()[0:2] == '##' :
		continue
	while platform.readline()[0:2] == '##' :
		continue
	
	paired_lines = izip_longest(vcf,platform,fillvalue='MISSING')
	for line_num, (record1,record2) in enumerate(paired_lines) :
		if "MISSING" in (record1,record2) :
			raise Exception(err2.format(options.vcf,platform_path))
		vals1 = record1.strip().split('\t')
		vals2 = record2.strip().split('\t')
		for i, field in enumerate(fields) :
			if vals1[i] != vals2[i] :
				raise Exception(err1.format(options.vcf,platform_path,field,line_num + 1,vals1[i],vals2[i]))

#lock all files (vcf, batch bcf, platform def)
#should log be locked? or want simultaneous adds of independent batches to have log entries intermixed?
#does this handle instances of process interruption (e.g. while vcf being converted to bcf)?

#create rolling log
AddLogger = logging.getLogger()
AddLogger.setLevel(logging.DEBUG)
#Doesn't seem to be a way to set backupCount to infinite, tried 0 and -1, both effectively mean no rolling, just single file that gets overwritten when full
#Setting to unreasonably high # to account for this
RollingLog = logging.handlers.RotatingFileHandler(options.bank + '/Logs/Add.log', maxBytes=50000000, backupCount=5000)
#will this die if can't write to log?
LogRecordFormat = logging.Formatter(fmt='%(process)d\t%(levelname)s\t%(asctime)s\t%(message)s')
RollingLog.setFormatter(LogRecordFormat)
AddLogger.addHandler(RollingLog)

#start log with user and arguments
if options.force :
	AddLogger.info('[%s;%s] is overwriting batch [%s] in bank [%s] from platform [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], options.batch, options.bank, options.platform)
else :
	AddLogger.info('[%s;%s] is adding batch [%s] to bank [%s] from platform [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], options.batch, options.bank, options.platform)

#Attempt to copy
try:
	subprocess.check_output(["cp", options.vcf, options.bank + '/Data/' + options.batch + '.vcf.gz'],stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
	AddLogger.error('Failed to copy [%s] into bank [%s] as [%s]', options.vcf, options.bank + '/Data/', options.batch + '.vcf.gz')
	AddLogger.error('Subprocess copy error: [%s]',e.output.strip())
	raise Exception('Failed to copy [%s] into bank [%s] as [%s] error: [%s]', options.vcf, options.bank + '/Data/', options.batch + '.vcf.gz',e.output)

#Attempt to tabix
try:
	subprocess.check_output(["tabix", "-p", "vcf", options.bank + '/Data/' + options.batch + '.vcf.gz'],stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
	AddLogger.error('Failed to tabix [%s]', options.bank + '/Data/' + options.batch + '.vcf.gz')
	AddLogger.error('Subprocess tabix error: [%s]',e.output.strip())
	raise Exception('Failed to tabix [%s] error: [%s]', options.bank + '/Data/' + options.batch + '.vcf.gz',e.output)
	
os.chown(options.bank + '/Data/' + options.batch + '.vcf.gz', -1, 2593)
os.chmod(options.bank + '/Data/' + options.batch + '.vcf.gz', 0660)

#Compare input and output
try:
	subprocess.check_output(["diff", options.vcf, options.bank + '/Data/' + options.batch + '.vcf.gz'],stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
	AddLogger.error('File in bank [%s] does not match inpput [%s]', options.bank + '/Data/' + options.batch + '.vcf.gz', options.vcf)
	AddLogger.error('Subprocess diff error: [%s]',e.output.strip())
	raise Exception('File in bank [%s] does not match input [%s] error: [%s]', options.bank + '/Data/' + options.batch + '.vcf.gz', options.vcf, e.output)
	
#Count sample and variant records
try:
	count = subprocess.check_output([bankconstants.vcftools, "--gzvcf", options.bank + '/Data/' + options.batch + '.vcf.gz'],stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
	AddLogger.error('Error counting records from [%s]', options.bank + '/Data/' + options.batch + '.vcf.gz')
	AddLogger.error('Subprocess vcftools error: [%s]',e.output.strip())
	raise Exception('Error counting records from [%s] error: [%s]', options.bank + '/Data/' + options.batch + '.vcf.gz', e.output)
count = count.strip().split('\n')
for log in count :
	if re.search(r"kept (\d+) out of \d+ Individuals",log) :
		sample = re.search(r"kept (\d+) out of \d+ Individuals",log)
	elif re.search(r"kept (\d+) out of a possible \d+ Sites",log) :
		var = re.search(r"kept (\d+) out of a possible \d+ Sites",log)
if sample and var :
	AddLogger.info('batch [%s] added with [%s] samples and [%s] variant records', options.batch , sample.group(1), var.group(1))
else :
	AddLogger.error('Error counting records from [%s]', options.bank + '/Data/' + options.batch + '.vcf.gz')
	AddLogger.error('vcftools log where expecting counts: [%s]',re.sub(r"\t"," ",";".join(count))
	raise Exception('Error counting records from [%s] vcftools log: [%s]', options.bank + '/Data/' + options.batch + '.vcf.gz', count)


