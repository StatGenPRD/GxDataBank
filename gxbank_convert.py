#!/GWD/appbase/tools/bin/python2.7
from optparse import OptionParser
## deprecated, but argparse requires python >= 2.7 does not work on GSK servers
import logging
import sys
import re
import gzip
import os
import pwd
import subprocess
import time
from itertools import izip_longest
import bankfunctions

#define paths
import bankconstants

start_time = time.time()

#---------------------------------------------------------------------------------------------
# ARGUMENTS
#---------------------------------------------------------------------------------------------
parser = OptionParser(description = 'usage: %prog OPTIONS')
#Add flag -u to leave unsorted VCF for test cases and then add system call to delete the unsorted after sort finishes
parser.add_option('-a', '--axiom-files', help = 'Path of the directory containing the gzip\'d AxiomGT1 files to be converted and the Ps.performance.txt.gz from SNPolisher',
                  action = 'store', type = 'string', dest = 'axiom', default = '')
parser.add_option('-p', '--platform-name', help = 'Name and version of the platform e.g. GSKBB1_v1',
                  action = 'store', type = 'string', dest = 'platform', default = '')
parser.add_option('-v', '--vcf', help = 'Name of the bgzip\'d VCF file to create (just the root, .vcf.gz will be added). This will also be recorded as the batchName in the VCF file header',
                  action = 'store', type = 'string', dest = 'vcf', default = '')
parser.add_option('-o', '--vcf-path', help = 'Path of directory where to write the gzip\'d VCF file, intermediary files, and the log [default: working directory]',
                  action = 'store', type = 'string', dest = 'vcfpath', default = os.getcwd())
parser.add_option('-d', '--data-bank', help = 'Path to data bank containing platform definition [default: Production]',
                  action = 'store', type = 'string', dest = 'bank', default = bankconstants.prodbank)
parser.add_option('-f', '--force', help = 'Force overwrite if VCF file already exists',
                  action = 'store_true', dest = 'force', default = False)
(options, args) = parser.parse_args()

#---------------------------------------------------------------------------------------------
# VALIDATE INPUT
#---------------------------------------------------------------------------------------------

#Get absolute paths for inclusion in log / error messages
outdir = os.path.realpath(options.vcfpath)
bank = os.path.realpath(options.bank)
axiomdir = os.path.realpath(options.axiom)

#validate vcf name: allow alphanum, dashes, and underscores
if options.vcf == '' :
	raise Exception('ERROR: vcf name cannot be blank')
elif re.search(r"^([\w\-])+$", options.vcf) == None :
	raise Exception('ERROR: vcf name may only contain alphanumeric characters, dash or underscore')

#validate platform name - must match an existing platform definition VCF
platformVCF_path = os.path.join(options.bank, 'GxDataBankPlatforms', options.platform, 'unsorted.vcf.gz')
if not os.access(platformVCF_path, os.R_OK) :
	raise Exception ('ERROR: Could not find or read platform definition [' + platformVCF_path + ']')

#validate axiom files - files must exist & be readable
callFile = os.path.join(options.axiom, 'AxiomGT1.calls.txt.gz')
confFile = os.path.join(options.axiom, 'AxiomGT1.confidences.txt.gz')
perfFile = os.path.join(options.axiom, 'Ps.performance.txt.gz')
vcfFile = os.path.join(options.vcfpath, options.vcf + '.vcf.gz')
vcfFileUnsorted = os.path.join(options.vcfpath, options.vcf + '_unsorted.vcf.gz')

#Decided not to include intensity data in VCF as inflates size (~1 GB vs ~9 GB), requiring longer to create, sort and validate.
#Instead, will store the summary file as-is beside the VCF in the bank
files = (callFile, confFile, perfFile)
for file in files :
	if not os.access(file, os.R_OK) :
		raise Exception('ERROR: Could not find or read AxiomGT1 file [' + file + ']')

#validate vcf path - must be valid path to a writable directory
if not os.access(options.vcfpath, os.W_OK) :
	raise Exception('ERROR: user cannot write to vcf-path [' + options.outdir + ']')
	
#check if vcf already exists and not forcing
if os.path.isfile(vcfFile) and not options.force :
	raise Exception('ERROR: vcf [' + vcfFile + '] already exists in vcf-path [' + options.outdir + ']')

#---------------------------------------------------------------------------------------------
# LOGGING
#---------------------------------------------------------------------------------------------

#create log - stand alone, not part of bank log
logging.basicConfig(filename=os.path.join(options.vcfpath, options.vcf + '.log'), filemode='w', format='%(levelname)s\t%(asctime)s\t%(message)s', level=logging.DEBUG)

#start log with user and arguments - will validate file contents as process
if options.force :
	logging.info('[%s;%s] is overwriting vcf [%s.vcf.gz] in directory [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], options.vcf, outdir)
else :
	logging.info('[%s;%s] is creating vcf [%s.vcf.gz] in directory [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], options.vcf, outdir)
logging.info('From AxiomGT1 files here [%s]', axiomdir)
logging.info('Using definition of platform [%s] in bank [%s]', options.platform, bank)


##
# "UNSORT" Ps.performance for GSKBB2 so that "special" APOE probeset is at bottom as it is in the AxiomGT files
# Affy's guidance for this special probeset is to call using different clustering algorithm parameters requiring separate processing from the rest of the array
# Their guidance is to manually append the AxiomGT results from this special probeset onto the end of the AxiomGT files from the rest of the array
# When this combined set of AxiomGT files (the whole array with special probeset appended) is run through SNPolisher's Ps_Classification function, it reports results in sorted order
# Which results in the special probset no longer being at the bottom of the Ps.performance file and thus this file is inconsistent with the AxiomGT files.
##
perfFileUnsorted = os.path.join(options.vcfpath, 'Ps.performance_unsorted.txt')

if options.platform[0:6] == 'GSKBB2' :
        logging.info('Creating a modified version of Ps.performance with AX-95861335 at the bottom since this is GSKBB2')
        cmd = "zgrep -v '^AX-95861335' {0} > {1}; zgrep '^AX-95861335' {0} >> {1}; gzip -f {1}".format(perfFile,perfFileUnsorted)
        os.system(cmd)
        perfFile = perfFileUnsorted + '.gz'


vcffields = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]

err1 = "ERROR: Probeset IDs in AxiomGT1 files in [{0}] do not match at record {1}. {2}={3}; {4}={5}; {6}={7}; {8}={9}"
err2 = "ERROR: AxiomGT1 files in [{0}] have different line counts."
err3 = "ERROR: Probeset [{0}] not found in platform definition [{1}]."
err4 = "ERROR: AxiomGT1 files in [{0}] have different sample counts. {1}={2}; {3}={4}"
err5 = "ERROR: AxiomGT1 files in [{0}] have different sample orders, first mismatch found at column {1}: {2}={3}; {4}={5}"

##
# ANNOTATIONS for INFO field
##
annotation_flag = lambda key,value: key if value == '1' else None
annotation_num = lambda key,value: "{0}={1}".format(key,value) if value != 'NA' and value != '' else None
annotation_string = lambda key,value: "{0}={1}".format(key,value)

##
# FILTERS for FILTER field
##
filter_num = lambda key,value,criteria: key if value != 'NA' and value != '' and float(value) < criteria else None

##
# FILTER Thresholds
##
thresholds_homRO = {'0':0.6, '1':0.6, '2':0.3, '3':-0.9}

with gzip.open(callFile, 'rb') as calls, \
	gzip.open(confFile, 'rb') as confs, \
     gzip.open(perfFile, 'rb') as perfs, \
	  gzip.open(platformVCF_path,"rb") as plat_vcf, \
	  gzip.open(vcfFileUnsorted,"wb") as vcf :
	  
	#Remove headers
	bankfunctions.read_through_headers(perfs,'#')
	bankfunctions.read_through_headers(plat_vcf,'##')
	#Confirm order of samples in each Axiom file matches
	call = bankfunctions.read_through_headers(calls,'#').strip().split('\t')
	conf = bankfunctions.read_through_headers(confs,'#').strip().split('\t')
	#if len(call) != len(conf) or len(conf) != len(sum) :
	if len(call) != len(conf) :
		raise Exception(err4.format(axiomdir,'AxiomGT1.calls.txt.gz',len(call) - 1,'AxiomGT1.confidences.txt.gz',len(conf) - 1))
	for i in range(1,len(call)) :
		#if call[i] != conf[i] or conf[i] != sum[i] :
		if call[i] != conf[i] :
			raise Exception(err5.format(axiomdir,i,'AxiomGT1.calls.txt.gz',call[i],'AxiomGT1.confidences.txt.gz',conf[i]))
	
   #log number of samples
	logging.info('[%d] samples found in these Axiom files', len(call) - 1)

	#Start new vcf with headers
	headers = """##fileformat=VCFv4.1
##platformSource=file://us2us00013.corpnet2.com{platformVCF_path}
##axiomSource=file://us2us00013.corpnet2.com{axiomdir}
##batchName={vcf}
##fileAuthor="{author}"
##FILTER=<ID=fld,Description="FLD less than 3.6">
##FILTER=<ID=hs,Description="HetSO less than -0.1">
##FILTER=<ID=cr,Description="Call rate less than 95 percent">
##FILTER=<ID=hr,Description="HomRO less than 0.6|0.3|-0.9 for 1|2|3 clusters respectively">
##INFO=<ID=CR,Number=1,Type=Float,Description="Call rate">
##INFO=<ID=FLD,Number=1,Type=Float,Description="FLD">
##INFO=<ID=HFLD,Number=1,Type=Float,Description="HomFLD">
##INFO=<ID=HS,Number=1,Type=Float,Description="HetSO">
##INFO=<ID=HR,Number=1,Type=Float,Description="HomRO">
##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor allele count">
##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of clusters">
##INFO=<ID=HZ,Number=0,Type=Flag,Description="Hemizygous">
##INFO=<ID=SP,Number=1,Type=String,Description="SNPolisher class">
##INFO=<ID=SPB,Number=0,Type=Flag,Description="SNPolisher best probeset">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GC,Number=1,Type=Float,Description="Genotype Confidence">
"""

	vcf.write(headers.format(
		platformVCF_path  =platformVCF_path,
		axiomdir          =axiomdir,
		vcf               =options.vcf,
		author            =pwd.getpwuid(os.getuid())[4]
	))

        pattern = re.compile("\.cel", re.IGNORECASE)

        samples = [pattern.sub("", header) for header in call[1:]]
	vcf.write('\t'.join(vcffields) + '\t' + '\t'.join(samples) + '\n')

	matched_lines = izip_longest(calls,confs,perfs,plat_vcf,fillvalue='MISSING')
	for line_num, (call,conf,perf,plat) in enumerate(matched_lines) :

		#Check number of records from each source file matches
		if "MISSING" in (call,conf,perf,plat):
			raise Exception(err2.format(axiomdir))
      
		callval = call.strip().split('\t')
		confval = conf.strip().split('\t')
		perfval = perf.strip().split('\t')
		platval = plat.strip().split('\t')

		#Confirm all records refer to the same probeset
		if not all(x==callval[0] for x in (confval[0],perfval[0],platval[2])):
			raise Exception(err1.format(axiomdir,line_num,callFile,callval[0],confFile,confval[0],perfFile,perfval[0],bank + '/PlatformDefinitions/' + options.platform + '.vcf.gz',platval[2]))

		# SNP LEVEL: Define A and B Alleles
		field = platval[7]
		platform_info = dict(v.split("=") for v in field.split(";"))
		alleleA = platform_info['A']
		alleleB = platform_info['B']
		lookup = dict()
		lookup['-1'] = './.'
		lookup['0'] = alleleA + '/' + alleleA
		lookup['1'] = alleleA + '/' + alleleB
		lookup['2'] = alleleB + '/' + alleleB

		# Assemble Data
		data = [lookup[call] + ':' + conf for call, conf in zip(callval[1:],confval[1:])]

		#Assemble info field
		annotations = [
			annotation_string('CR',perfval[2]),
			annotation_string('HR',perfval[6]),
			annotation_string('MAC',perfval[7]),
			annotation_string('NC',perfval[8]),
			annotation_string('SP',perfval[15]),
			annotation_flag('HZ',perfval[13]),
			annotation_flag('SPB',perfval[16]),
			annotation_num('FLD',perfval[3]),
			annotation_num('HFLD',perfval[4]),
			annotation_num('HS',perfval[5])
		]

		info = [a for a in annotations if a is not None]

		#Assemble filter field
		filters = [
			filter_num('cr',perfval[2],95.0,
			filter_num('fld',perfval[3],3.6),
			filter_num('hs',perfval[5],-0.1),
			filter_num('hr',perfval[6],thresholds_homRO[perfval[8]])
		]

		fltr = [f for f in filters if f is not None]
		#Set filter to missing if nothing failed
		if len(fltr) == 0 :
			fltr.append('.')
		#Write record
		vcf.write('\t'.join(platval[0:6]) + '\t' + ';'.join(fltr) + '\t' + ';'.join(info) + '\tGT:GC\t' + '\t'.join(data) + '\n')

sorted_base = os.path.join(options.vcfpath,options.vcf + '_sorted')
sorted_vcf = sorted_base + '.vcf.bgz'
sorted_err = sorted_base + '.err'
# Ensure OUTDIR is used for TEMP output of SORT
os.environ['TMPDIR'] = outdir

#Note, the re-dir of stderr to sorted_err will only capture stderr from the last command (bgzip), not everything
cmd = "zcat {} | {} -t {} | bgzip > {} 2> {}".format(vcfFileUnsorted, bankconstants.vcf_sort, options.vcfpath, sorted_vcf, sorted_err)
os.system(cmd)

logging.info("vcf [%s] created", sorted_vcf)

# DECISION: Do not run VCF-VALIDATOR based on risk (LOW for routine process) vs cost (HIGH Long-running)

elapsed_time = time.time() - start_time
print "Elapsed Time: %.3f" % elapsed_time
