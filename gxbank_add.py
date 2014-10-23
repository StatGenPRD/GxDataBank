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
import bankfunctions

#define paths
import bankconstants

#extract arguments
parser = OptionParser(description = 'usage: %prog OPTIONS')
parser.add_option('-b', '--batch-name', help = 'Name of batch of data as supplied by the lab and recorded in Reflection',
                  action = 'store', type = 'string', dest = 'batch', default = '')
parser.add_option('-v', '--vcf', help = 'Path and name of the bgzip\'d VCF file to add',
                  action = 'store', type = 'string', dest = 'vcf', default = '')
parser.add_option('-s', '--summary-data', help = 'Path and name of the gzip\'d AxiomGT1.summary.txt file containing the intensity data for the clusterset used to generate the VCF file',
                  action = 'store', type = 'string', dest = 'summary', default = '')
parser.add_option('-p', '--posterior-data', help = 'Path and name of the gzip\'d AxiomGT1.posteriors.txt file containing the posteriors for the clusterset used to generate the VCF file',
                  action = 'store', type = 'string', dest = 'posteriors', default = '')
parser.add_option('-d', '--data-bank', help = 'Path to data bank [default: Production]',
                  action = 'store', type = 'string', dest = 'bank', default = bankconstants.prodbank)
parser.add_option('-f', '--force', help = 'Force overwrite if VCF already exists in bank',
                  action = 'store_true', dest = 'force', default = False)
(options, args) = parser.parse_args()

vcf_path = os.path.realpath(options.vcf)
tbi_path = vcf_path + '.tbi'
summary_path = os.path.realpath(options.summary)
posteriors_path = os.path.realpath(options.posteriors)

bank_path = os.path.realpath(options.bank)
bank_vcf_path = os.path.join(bank_path, 'Data', 'VCF', options.batch + '.vcf.gz')
bank_tbi_path = bank_vcf_path + '.tbi'
bank_summary_path = os.path.join(bank_path, 'Data', 'Other', options.batch + '.summary.txt.gz')
bank_posteriors_path = os.path.join(bank_path, 'Data', 'Other', options.batch + '.posteriors.txt.gz')

source_files = [vcf_path, summary_path, posteriors_path]
dest_files = [bank_vcf_path, bank_summary_path, bank_posteriors_path, bank_tbi_path]

#validate batch name: allow alphanum, dashes, and underscores
if options.batch == '' :
	raise Exception('ERROR: batch-name cannot be blank')
elif re.search(r"^([\w\-])+$", options.batch) == None :
	raise Exception('ERROR: batch-name may only contain alphanumeric characters, dash or underscore')

#validate source files - must be valid paths to files that can be read
for source_file in source_files:
        if not os.access(source_file, os.R_OK) :
                raise Exception('ERROR: Could not find or read file [' + source_file + ']')
	
#validate data bank paths
#must be valid path to a writable directory
if not os.access(os.path.join(bank_path, 'Data', 'VCF'), os.W_OK) :
	raise Exception('ERROR: user cannot write to bank [' + bank_path + '/Data/VCF]')
if not os.access(os.path.join(bank_path, 'Data', 'Other'), os.W_OK) :
	raise Exception('ERROR: user cannot write to bank [' + bank_path + '/Data/Other]')

#check if batch name already exists in bank and not forcing
for filename in [bank_vcf_path,bank_tbi_path,bank_summary_path] :
	if os.path.isfile(filename) and not options.force :
		raise Exception('ERROR: batch-name [' + options.batch + '] already loaded to this bank [' + filename + ']')

#Attempt to tabix - doing this outside bank as it can indicate malformed input (e.g. unsorted or non-bgzip'd vcf)
try:
	subprocess.check_output(["tabix", "-fp", "vcf", vcf_path],stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
	raise Exception('Failed to tabix [%s] error: [%s]' % (vcf_path, e.output))
source_files.append(vcf_path + '.tbi')	

#must pass vcf-validator
validator_result = subprocess.check_output(["perl", "-I", bankconstants.vcftools_base + '/perl', bankconstants.vcf_validator, vcf_path],stderr=subprocess.STDOUT)
validator_result = validator_result.strip().split('\n')
for error in validator_result :
	#Acceptable warning messages for missing contig definitions like:
	#The header tag 'contig' not present for CHROM=chr1. (Not required but highly recommended.)
	#and for un-defined variants like:
	#chr0:0 .. REF allele listed in the ALT field??
	#Our convention is to define these variants on chr0
	if re.search(r"Not required", error) == None and error != 'chr0:0 .. REF allele listed in the ALT field??' and error != '':
		raise Exception('ERROR: vcf [' + vcf_path + '] not valid [' + error + ']')

#must match platform definition - no extra/missing records
#CHROM,POS,ID,REF,ALT must match
#Same sort order - may be able to relax this requirement later
#samples must also match between the summary and vcf
err1 = "ERROR: VCF file [{0}] does not match PLATFORM file [{1}] on {2} at variant record {3}. VCF={4}; PLATFORM={5}"
err2 = "ERROR: VCF file [{0}] does not match PLATFORM file [{1}]. Different line counts."
fields = ["CHROM","POS","ID","REF","ALT"]
with gzip.open(vcf_path, 'rb') as vcf, gzip.open(platform_path, 'rb') as platform, gzip.open(summary_path, 'rb') as summary, gzip.open(posteriors_path, 'rb') as posteriors:
	#remove header comments from each file and capture header row
	vcfhead = bankfunctions.read_through_headers(vcf,'##')
	platformhead = bankfunctions.read_through_headers(platform,'##')
	summaryhead = bankfunctions.read_through_headers(summary,'#')
	posteriorshead = bankfunctions.read_through_headers(posteriors,'#')
	
	#check samples match between vcf and summary
	vcf_samples = vcfhead.strip().split('\t')
	summary_samples = summaryhead.strip().split('\t')
	posteriors_samples = posteriorshead.strip().split('\t')
	if vcf_samples[9:] != summary_samples[1:] :
		raise Exception('ERROR: vcf + [' + vcf_path + '] samples do not match samples in summary [' + summary_path + ']')
        if vcf_samples[9:] != posterior_samples[1:] :
		raise Exception('ERROR: vcf + [' + vcf_path + '] samples do not match samples in posteriors file [' + posteriors_path + ']')
	sample_count = len(summary_samples) - 1
	
	#check variants match between vcf and platform
	var_count = 0
	paired_lines = izip_longest(vcf,platform,fillvalue='MISSING')
	for line_num, (record1,record2) in enumerate(paired_lines) :
		if "MISSING" in (record1,record2) :
			raise Exception(err2.format(vcf_path,platform_path))
		var_count += 1
		vals1 = record1.strip().split('\t')
		vals2 = record2.strip().split('\t')
		for i, field in enumerate(fields) :
			if vals1[i] != vals2[i] :
				raise Exception(err1.format(vcf_path,platform_path,field,line_num + 1,vals1[i],vals2[i]))
				
#lock all files (vcf, batch bcf, platform def)
#does this handle instances of process interruption (e.g. while vcf being converted to bcf)?

#create rolling log
AddLogger = logging.getLogger()
AddLogger.setLevel(logging.DEBUG)
RollingLog = logging.handlers.RotatingFileHandler(options.bank + '/Logs/Add.log', maxBytes=50000000, backupCount=5000)
LogRecordFormat = logging.Formatter(fmt='%(process)d\t%(levelname)s\t%(asctime)s\t%(message)s')
RollingLog.setFormatter(LogRecordFormat)
AddLogger.addHandler(RollingLog)

#start log with user and arguments
if options.force :
	AddLogger.info('[%s;%s] is overwriting batch [%s] in bank [%s] from platform [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], options.batch, bank_path, options.platform)
else :
	AddLogger.info('[%s;%s] is adding batch [%s] to bank [%s] from platform [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], options.batch, bank_path, options.platform)

#Attempt to copy vcf, tbi, and summary into bank
for source,dest in zip(source_files, dest_files) :
	try:
		subprocess.check_output(["cp", '-f', source, dest],stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError as e:
		AddLogger.error('Failed to copy [%s] into bank as [%s]', source, dest)
		AddLogger.error('Subprocess copy error: [%s]',e.output.strip())
		raise Exception('Failed to copy [%s] into bank as [%s] error: [%s]' % (source, dest, e.output))

#Set group and permissions of files copied to bank
for filename in dest_files :
	os.chown(filename, -1, 2593)
	os.chmod(filename, 0440)

#Compare input and output
for source,dest in zip(source_files, dest_files):
	try:
		subprocess.check_output(["diff", source, dest],stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError as e:
		AddLogger.error('File in bank [%s] does not match inpput [%s]', dest, source)
		AddLogger.error('Subprocess diff error: [%s]',e.output.strip())
		raise Exception('File in bank [%s] does not match input [%s] error: [%s]' % (dest, source, e.output))

AddLogger.info('batch [%s] added with [%s] samples and [%s] variant records', options.batch , str(sample_count), str(var_count))



	
