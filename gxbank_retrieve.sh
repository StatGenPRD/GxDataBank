#! /bin/bash

#######
#Input parameters
#######
#Raw data in bank - assumes all data in single clusterset
RAW=$1
#Study ID used to name files
STUDY=$2
#Output directory
OUTDIR=$3
#paths to programs
vcftools=/GWD/appbase/projects/GXapp/vcftools/vcftools_0.1.12b/bin/vcftools
plink64=/GWD/bioinfo/apps/bin/plink64
plink2=/GWD/appbase/projects/GXapp/plink/plink-1.9.0_Jan-12-2015/plink
#Platform definition - make sure latest from git repo
PLATFORM=/GWD/appbase/projects/statgen4/GxBank/Prod/GxDataBankPlatforms/$4/unsorted.vcf.gz
#File containing CEL (File Registration) to USUBJ ID mapping of only samples to keep (i.e. best call rate for duplicates, no QC failures) in plink update format (FID IID FID IID)
CELtoUSUBJ=$5

#######
#Execution
#######
#create output directory as necessary
mkdir -p $OUTDIR

#Initiate log
#Note, log counts are just what is readily countable - not necessarily what is most relevant/useful/comprehensible
echo "Reference genome version: `zgrep '^##reference=file:' $PLATFORM`" >$OUTDIR/RETRIEVE.log
echo "Assays in bank platform: `zgrep -v '^#' $PLATFORM | wc -l`" >>$OUTDIR/RETRIEVE.log

#Get passing SNPolisher classes from raw data
#Note, using these 3 classes, it is not possible to fail any of the clustering metrics so the FILTER field is irrelevant
zgrep -e '^#|SP=PolyHighResolution|SP=MonoHighResolution|SP=NoMinorHom' $RAW | bgzip >$OUTDIR/SPClass.vcf.bgz

#Get passing probeset from platform definition - make sure latest version from git repo
#Note, although this is an invalid VCF file as it is unsorted, this specific FILTER based operation seems to work as expected
$vcftools --gzvcf $PLATFORM --remove-filtered-all --recode --stdout --out filter_platform | grep -v '^#' | cut -f 3 >$OUTDIR/platform_pass.ps
echo "Assays considered usable from platform a priori: `cat $OUTDIR/platform_pass.ps | wc -l`" >>$OUTDIR/RETRIEVE.log
echo "Assays in bank with data: `zgrep -v '^#' $RAW | cut -f 1 | wc -l`" >>$OUTDIR/RETRIEVE.log
echo "Assays passing lab QC: `zgrep -v '^#' $OUTDIR/SPClass.vcf.bgz | cut -f 1 | wc -l`" >>$OUTDIR/RETRIEVE.log

#Determine non-best probesets where replicates per SNPolisher results in data
#Note, for non-replicate probesets, will still have SPB flag indicating the sole probeset is the best
zgrep -v 'SPB' $OUTDIR/SPClass.vcf.bgz | cut -f 3 >$OUTDIR/not_SPB.ps
echo "Poorer performing of replicate assays excluded: `cat $OUTDIR/not_SPB.ps | wc -l`" >>$OUTDIR/RETRIEVE.log

#Bug in plink2 that --keep happens BEFORE --update-ids. Reported to developers but workaround to use pre-update IDs
cut -f 1-2 $CELtoUSUBJ >$OUTDIR/CEL.keep
echo "Subjects requested: `cat $OUTDIR/CEL.keep | wc -l`" >>$OUTDIR/RETRIEVE.log

$plink2 --vcf $OUTDIR/SPClass.vcf.bgz --exclude $OUTDIR/not_SPB.ps --extract $OUTDIR/platform_pass.ps --double-id --make-bed --update-ids $OUTDIR/CEL_to_USUBJ.txt --out $OUTDIR/$STUDY --keep CEL.keep --split-x hg19 no-fail
echo "Subjects retrieved: `cat $OUTDIR/$STUDY.fam | wc -l`" >>$OUTDIR/RETRIEVE.log
echo "Assays retrieved: `cat $OUTDIR/$STUDY.bim | wc -l`" >>$OUTDIR/RETRIEVE.log



######
#Alternative approach to use vcftools instead of plink2
#Takes ~16 min for a clusterset instead of ~2 min for plink2
#Unlike plink2, removes 9 ps that interrogate 2 non-ref alleles as non-bi-allelic even though the ref allele is not interrogated by a separate probeset (so from GSKBB2 perspective, it is bi-allelic and thus not filtered as oa)
######
#$vcftools --gzvcf $OUTDIR/SPClass.vcf.bgz --keep-INFO SPB --snps $OUTDIR/platform_pass.ps --plink-tped --out $OUTDIR/BestPS
#cut -f 3-4 $OUTDIR/CEL_to_USUBJ.txt >$OUTDIR/USUBJ.keep
#$plink64 --tfile BestPS --make-bed --update-ids $OUTDIR/CEL_to_USUBJ.txt --keep $OUTDIR/USUBJ.keep  --out $OUTDIR/$STUDY