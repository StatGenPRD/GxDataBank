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
from itertools import izip_longest
import bankfunctions

#define paths
import bankconstants

#extract arguments
parser = OptionParser(description = 'usage: %prog OPTIONS')
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
#Get absolute paths for inclusion in log / error messages
outdir = os.path.abspath(options.vcfpath)
bank = os.path.abspath(options.bank)
axiomdir = os.path.abspath(options.axiom)

#validate vcf name: allow alphanum, dashes, and underscores
if options.vcf == '' :
	raise Exception('ERROR: vcf name cannot be blank')
elif re.search(r"^([\w\-])+$", options.vcf) == None :
	raise Exception('ERROR: vcf name may only contain alphanumeric characters, dash or underscore')

#validate platform name - must match an existing platform definition with both VCF and AB_RefAlt map
platformVCF_path = options.bank + '/GxDataBankPlatforms/' + options.platform + '/unsorted.vcf.gz'
if not os.access(platformVCF_path, os.R_OK) :
	raise Exception ('ERROR: Could not find or read platform definition [' + platformVCF_path + ']')
platformAB_path = options.bank + '/GxDataBankPlatforms/' + options.platform + '/AB_RefAlt_map.txt.gz'
if not os.access(platformAB_path, os.R_OK) :
	raise Exception ('ERROR: Could not find or read platform definition [' + platformAB_path + ']')

#validate axiom files - files must exist & be readable
callFile = '/AxiomGT1.calls.txt.gz'
confFile = '/AxiomGT1.confidences.txt.gz'
summFile = '/AxiomGT1.summary.txt.gz'
postFile = '/AxiomGT1.snp-posteriors.txt.gz'
perfFile = '/Ps.performance.txt.gz'
#Decided not to include intensity data in VCF as inflates size (~1 GB vs ~9 GB), requiring longer to create, sort and validate.
#Instead, will store the summary file as-is beside the VCF in the bank
#files = (callFile, confFile, summFile, postFile, perfFile)
files = (callFile, confFile, postFile, perfFile)
for file in files :
	if not os.access(options.axiom + file, os.R_OK) :
		raise Exception('ERROR: Could not find or read AxiomGT1 file [' + options.axiom + '/' + file + ']')

#validate vcf path - must be valid path to a writable directory
if not os.access(options.vcfpath, os.W_OK) :
	raise Exception('ERROR: user cannot write to vcf-path [' + options.outdir + ']')
	
#check if vcf already exists and not forcing
if os.path.isfile(options.vcfpath + '/' + options.vcf + '.vcf.gz') and not options.force :
	raise Exception('ERROR: vcf [' + options.vcf + '.vcf.gz] already exists in vcf-path [' + options.outdir + ']')
	
#create log - stand alone, not part of bank log
logging.basicConfig(filename=options.vcfpath + '/' + options.vcf + '.log', filemode='w', format='%(levelname)s\t%(asctime)s\t%(message)s', level=logging.DEBUG)

#start log with user and arguments - will validate file contents as process
if options.force :
	logging.info('[%s;%s] is overwriting vcf [%s.vcf.gz] in directory [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], options.vcf, outdir)
else :
	logging.info('[%s;%s] is creating vcf [%s.vcf.gz] in directory [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], options.vcf, outdir)
logging.info('From AxiomGT1 files here [%s]', axiomdir)
logging.info('Using definition of platform [%s] in bank [%s]', options.platform, bank)
	
vcffields = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]

#err1 = "ERROR: Probeset IDs in AxiomGT1 files in [{0}] do not match at record {1}. {2}={3}; {4}={5}; {6}={7} and {8}; {9}={10}; {11}={12}; {13}={14}; {15}={16}"
err1 = "ERROR: Probeset IDs in AxiomGT1 files in [{0}] do not match at record {1}. {2}={3}; {4}={5}; {6}={7}; {8}={9}; {10}={11}; {12}={13}"
err2 = "ERROR: AxiomGT1 files in [{0}] have different line counts."
err3 = "ERROR: Probeset [{0}] not found in platform definition [{1}]."
#err4 = "ERROR: AxiomGT1 files in [{0}] have different sample counts. {1}={2}; {3}={4}; {5}={6}"
err4 = "ERROR: AxiomGT1 files in [{0}] have different sample counts. {1}={2}; {3}={4}"
#err5 = "ERROR: AxiomGT1 files in [{0}] have different sample orders, first mismatch found at column {1}: {2}={3}; {4}={5}; {6}={7}"
err5 = "ERROR: AxiomGT1 files in [{0}] have different sample orders, first mismatch found at column {1}: {2}={3}; {4}={5}"
err6 = "ERROR: Platform [{0}] definitions in [{1}] are inconsistent at line {3}, AB map is bi-allelic {4} but VCF is multi-allelic {5}"

with gzip.open(options.axiom + callFile, 'rb') as calls, gzip.open(options.axiom + confFile, 'rb') as confs, \
	  gzip.open(options.axiom + postFile, 'rb') as posts, gzip.open(options.axiom + perfFile, 'rb') as perfs, \
	  gzip.open(platformAB_path,"rb") as orients, \
	  gzip.open(platformVCF_path,"rb") as plat_vcf, \
	  gzip.open(options.vcfpath + '/' + options.vcf + '_unsorted.vcf.gz',"wb") as vcf :
	  #gzip.open(options.axiom + summFile, 'rb') as summ
	  
	#Remove headers
	bankfunctions.read_through_headers(perfs,'#')
	bankfunctions.read_through_headers(posts,'#')
	bankfunctions.read_through_headers(plat_vcf,'##')
	orients.readline()
	#Confirm order of samples in each Axiom file matches
	call = bankfunctions.read_through_headers(calls,'#').strip().split('\t')
	conf = bankfunctions.read_through_headers(confs,'#').strip().split('\t')
	#sum = bankfunctions.read_through_headers(summ,'#').strip().split('\t')
	#if len(call) != len(conf) or len(conf) != len(sum) :
	if len(call) != len(conf) :
		raise Exception(err4.format(axiomdir,'AxiomGT1.calls.txt.gz',len(call) - 1,'AxiomGT1.confidences.txt.gz',len(conf) - 1))
	for i in range(1,len(call)) :
		#if call[i] != conf[i] or conf[i] != sum[i] :
		if call[i] != conf[i] :
			raise Exception(err5.format(axiomdir,i,'AxiomGT1.calls.txt.gz',call[i],'AxiomGT1.confidences.txt.gz',conf[i]))
		#removing .CEL extension to match File Registration value in Reflection
		#unsure if replace is case sensitive, re might be better approach in case file contains .CEL within root name
		call[i] = call[i].replace('.CEL','',1)
	#log number of samples
	logging.info('[%d] samples found in these Axiom files', len(call) - 1)
	#Start new vcf with headers
	vcf.write('##fileformat=VCFv4.1\n')
	vcf.write('##platformSource=file://us2us00013.corpnet2.com' + platformVCF_path + '\n')
	vcf.write('##axiomSource=file://us2us00013.corpnet2.com' + axiomdir + '\n')
	vcf.write('##batchName=' + options.vcf + '\n')
	vcf.write('##fileAuthor="' + pwd.getpwuid(os.getuid())[4] + '"\n')
	vcf.write('##FILTER=<ID=fld,Description="FLD less than 3.6">\n')
	vcf.write('##FILTER=<ID=hs,Description="HetSO less than -0.1">\n')
	vcf.write('##FILTER=<ID=cr,Description="Call rate less than 0.95">\n')
	vcf.write('##FILTER=<ID=hr,Description="HomRO less than 0.6|0.3|-0.9 for 1|2|3 clusters respectively">\n')
	vcf.write('##INFO=<ID=CR,Number=1,Type=Float,Description="Call rate">\n')
	vcf.write('##INFO=<ID=FLD,Number=1,Type=Float,Description="FLD">\n')
	vcf.write('##INFO=<ID=HFLD,Number=1,Type=Float,Description="HomFLD">\n')
	vcf.write('##INFO=<ID=HS,Number=1,Type=Float,Description="HetSO">\n')
	vcf.write('##INFO=<ID=HR,Number=1,Type=Float,Description="HomRO">\n')
	vcf.write('##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor allele count">\n')
	vcf.write('##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of clusters">\n')
	vcf.write('##INFO=<ID=HZ,Number=0,Type=Flag,Description="Hemizygous">\n')
	vcf.write('##INFO=<ID=SP,Number=1,Type=String,Description="SNPolisher class">\n')
	vcf.write('##INFO=<ID=SPB,Number=0,Type=Flag,Description="SNPolisher best probeset">\n')
	vcf.write('##INFO=<ID=PC,Number=.,Type=String,Description="Posterior clusters BB|AB|AA|CV. Expect a second set of values for chrX non-PAR variants to represent the haploid clusters used to call males">\n')
	vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	vcf.write('##FORMAT=<ID=GC,Number=1,Type=Float,Description="Genotype Confidence">\n')
	#vcf.write('##FORMAT=<ID=AI,Number=1,Type=Float,Description="Allele A intensity">\n')
	#vcf.write('##FORMAT=<ID=BI,Number=1,Type=Float,Description="Allele B intensity">\n')
	vcf.write('\t'.join(vcffields) + '\t' + '\t'.join(call[1:]) + '\n')
	#Posterior file will have 2 records for chrX variants, one for the hemizygous calling of males and one for diploid calling of females
	#Will use a prior/next strategy to identify these and merge as appropriate, initializing prior value here
	priorpost = posts.readline().strip().split('\t')
	matched_lines = izip_longest(calls,confs,perfs,orients,plat_vcf,fillvalue='MISSING')
	for line_num, (call,conf,perf,orient,plat) in enumerate(matched_lines) :
		#for testing
		#if line_num > 1000 :
		#	break
		"""
		#Summary file has 2 lines per probeset, one for each of the A and B signals
		sum1val = summ.readline().strip().split('\t')
		#Non-polymorphic probes only have one line and are not present in the other files (no genotypes)
		while re.match(r"^AFFX-NP",sum1val[0]) != None :
			sum1val = summ.readline().strip().split('\t')
		#Assuming the line for the B allele is immediately after the A and thus don't expect a non-polymorphic probe between
		sum2val = summ.readline().strip().split('\t')
		"""
		#Hemizygous posteriors suffixed with :1
		hemi = 0 #to track which instance is the hemizygous
		priorpostps = re.search(r"^(.+):1$",priorpost[0])
		if priorpostps != None :
			priorpostps = priorpostps.group(1)
			hemi = 1
		else :
			priorpostps = priorpost[0]
		nextpost = posts.readline().strip().split('\t')
		nextpostps = re.search(r"^(.+):1$",nextpost[0])
		if nextpostps != None :
			nextpostps = nextpostps.group(1)
			hemi = 2
		else :
			nextpostps = nextpost[0]
		#If next posterior record matches prior, merge with prior
		if nextpostps == priorpostps :
			if hemi == 2 :
				priorpost = priorpost + nextpost[1:]
			else :
				priorpost = nextpost + priorpost[1:]
			nextpost = posts.readline().strip().split('\t')
		#Check number of records from each source file matches
		#if "MISSING" in (call,conf,perf,orient,plat) or sum1val[0] == '' or sum2val[0] == '' or priorpost[0] == '':
		if "MISSING" in (call,conf,perf,orient,plat) or priorpost[0] == '':
			raise Exception(err2.format(axiomdir))
		callval = call.strip().split('\t')
		confval = conf.strip().split('\t')
		perfval = perf.strip().split('\t')
		orientval = orient.strip().split('\t')
		platval = plat.strip().split('\t')
		#Confirm all records refer to the same probeset
		if (callval[0] != confval[0] or confval[0] 
		    #!= sum1val[0][0:-2] or sum1val[0][0:-2] != sum2val[0][0:-2] or sum2val[0][0:-2] 
			!= priorpostps or
			priorpostps != perfval[0] or perfval[0] != orientval[0] or orientval[0] != platval[2]) :
			raise Exception(err1.format(axiomdir,line_num,callFile,callval[0],confFile,confval[0],
							postFile,priorpostps,metricFile,
							metval[0],perfFile,perfval[0],options.platform + '.AB_RefAlt_map.txt',orientval[0],
							bank + '/PlatformDefinitions/' + options.platform + '.vcf.gz',platval[2]))
		"""
		#Determine which summary record is from allele A and which from B
		if sum1val[0][-1] == 'A' :
			sumA = sum1val
			sumB = sum2val
		else :
			sumA = sum2val
			sumB = sum1val
		"""
		#Assemble data
		data = []
		for i in range(1,len(callval)) :
			#No call
			if callval[i] == '-1' :
				gt = './.'
			else :
				#bi-allelic
				if orientval[1] == 'REF' or orientval[2] == 'REF' :
					#Confirm not multi-allelic in platform VCF
					if len(platval[4].split(',')) > 1 :
						raise Exception(err6.format(options.platform,bank,line_num,'\t'.join(orientval),'\t'.join(platval)))
					if orientval[1] == 'REF' :
						alleleA = '0'
						alleleB = '1'
					else :
						alleleA = '1'
						alleleB = '0'
				#Multi-allelic, lookup allele index
				else :
					alleles = [platval[3]] + platval[4].split(',')
					for index,nt in enumerate(alleles) :
						if orientval[1] == nt :
							alleleA = str(index)
						if orientval[2] == nt :
							alleleB = str(index)
				#AA
				if callval[i] == '0' :
					gt = alleleA + '/' + alleleA
				#AB
				elif callval[i] == '1' :
					gt = alleleA + '/' + alleleB
				#BB
				else :
					gt = alleleB + '/' + alleleB
			#data.append(gt + ':' + confval[i] + ':' + sumA[i] + ':' + sumB[i])
			data.append(gt + ':' + confval[i])
		#Assemble info field
		info = ['CR=' + perfval[2],'HR=' + perfval[6],'MAC=' + perfval[7],
				'NC=' + perfval[8],'SP=' + perfval[15], 'PC=' + '|'.join(priorpost[1:])]
		if perfval[13] == '1' :
			info.append('HZ')
		if perfval[3] != 'NA' :
			info.append('FLD=' + perfval[3])
		if perfval[4] != 'NA' :
			info.append('HFLD=' + perfval[4])
		if perfval[5] != 'NA' :
			info.append('HS=' + perfval[5])
		if perfval[16] == '1' :
			info.append('SPB')
		#Assemble filter field
		filter = []
		if float(perfval[2]) < 0.95 :
			filter.append('cr')
		if perfval[3] != 'NA' and float(perfval[3]) < 3.6 :
			filter.append('fld')
		if perfval[5] != 'NA' and float(perfval[5]) < -0.1 :
			filter.append('hs')
		if (perfval[8] == '1' and float(perfval[6]) < 0.6) or (perfval[8] == '2' and float(perfval[6]) < 0.3) or (perfval[8] == '3' and float(perfval[6]) < -0.9) :
			filter.append('hr')
		#Set filter to missing if nothing failed
		if len(filter) == 0 :
			filter.append('.')
		#Write record
		#vcf.write('\t'.join(platval[0:6]) + '\t' + ';'.join(filter) + '\t' + ';'.join(info) + '\tGT:GC:AI:BI\t' + '\t'.join(data) + '\n')
		vcf.write('\t'.join(platval[0:6]) + '\t' + ';'.join(filter) + '\t' + ';'.join(info) + '\tGT:GC\t' + '\t'.join(data) + '\n')
		#increment priorpost
		priorpost = nextpost

#sort results - with intensity, takes ~1 hour, w/o, takes ~15 min
with open(options.vcfpath + '/' + options.vcf + '_sorted.vcf',"w") as sortedvcf, open(options.vcfpath + '/sort.err', "w") as sorterr :
	subprocess.call(["perl", "-I", bankconstants.vcftools_base + '/perl', bankconstants.vcf_sort, "-t", options.vcfpath, options.vcfpath + '/' + options.vcf + '_unsorted.vcf.gz', ],stdout=sortedvcf, stderr=sorterr)
#check any sort err
with open(options.vcfpath + '/sort.err', "r") as sorterr :
	error = sorterr.read().strip()
	if len(error) > 0 :
		raise Exception('ERROR: sorting failed: %s',error)

#Try to bgzip
try:
	subprocess.check_output(["bgzip", "-f", options.vcfpath + '/' + options.vcf + '_sorted.vcf'],stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
	raise Exception("ERROR: bgzip'ing the sorted file failed: %s", e.output)
os.rename(options.vcfpath + '/' + options.vcf + '_sorted.vcf.gz', options.vcfpath + '/' + options.vcf + '.vcf.gz')

#Restrict to user access since don't know what group is or where user is working so allowing group access could make data accessible to those not authorized
os.chmod(options.vcfpath + '/' + options.vcf + '.vcf.gz', 0600)
#Remove sort error file and pre-sorted vcf
os.remove(options.vcfpath + '/' + options.vcf + '_unsorted.vcf.gz')
os.remove(options.vcfpath + '/sort.err')

logging.info("vcf [%s] created", outdir + '/' + options.vcf + '.vcf.gz')

#Run vcf-validator? It is an o/n job and will be done as part of bank load so may not be appropriate.
#With intensity, takes ~10.5 hours, without ~7.25
#Convert chrY, chrM, and male chrX to haploid ? Unclear there will be much added value

#Checked the first PGx6712 clusterset against the plink export from Genotyping Console provided by BSA.  
#The NORM probesets were exported from GTC as missing (variant annotations including alleles likely missing from the config files used by GTC) so could not be compared.
#All other probesets were 100% concordant using plink merge-mode 6 (including NoCalls in comparison)