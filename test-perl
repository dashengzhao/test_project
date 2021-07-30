#!/usr/bin/perl
#2013.07.29

use strict;
use warnings;

my (@line,%para,@tmp,@samples,@samplesName);
my ($numOfSamples,@refDir,$refSeqName,$isPairedEnd);

open IN,"configuration" or die("\nError: Can't open the file: configuration\n");
while (<IN>) {
	if(! /^\s*$|^\s*\#/){
		$_=~/(\w+)=([\w\/\.]+)/;
		$para{$1}=$2;
	}
}
close IN;

##Check parameters
if($para{'numOfThreads'}<4){
	$para{'numOfThreads'}=4;
	print "\nthe number of used threads has been changed to 4 automatically.\n";
}

print "\n-----------------\n";
print "Custom parameters\n";
print "-----------------\n";
print sprintf("%-25s","numOfThreads").$para{'numOfThreads'}."\n";
print "\n";

##check input folder
if(! -e "input"){die("\nError: you should prepare the resequencing data in fastq format.\n");}
$numOfSamples=0;
@samples=glob("input/*.fq.gz");
@samples=sort @samples;
if($samples[0]=~/_[FR]/){
	$isPairedEnd=1;
	$numOfSamples=@samples/2;
}
else{
	$isPairedEnd=0;
	$numOfSamples=@samples;
}
if($numOfSamples<1){
	if($para{'alignReads'}*$para{'createFinalBam'}!=0){
		die("\nThere is not any input file in the {input} directory OR the format of the input file's names doesn't match the standard (like 'g001_F.fq.gz').\n");
	}
}
else{
	&getSamplesName();
	print sprintf("%-25s","numOfSamples").$numOfSamples."\n";
}

##check ref folder
if(! -e "ref"){die("\nError: you should prepare the genome reference sequence in fasta format, named as [ref].fa.\n");}
@refDir=glob("ref/*.fa");
if(@refDir<1){
	die("\nError: you should prepare the genome reference sequence in fasta format, named as [ref].fa.\n");
}
elsif(@refDir==1){
	@tmp=split /\//,$refDir[0];
	$refSeqName=$tmp[1];
	print sprintf("%-25s","refSeqName").$refSeqName."\n";
}
else{
	die("\nError: there are more than one ref file.\n");
}

print "\n";

if(! -e 'log'){system('mkdir log');}
open LOG,">log/runPipeline.log";

#===========================
#[Index reference sequence]
#===========================

print "\n#==========================\n";
print "#[Index reference sequence]\n";
print "#==========================\n";
print LOG sprintf("%-60s","Index reference sequence...Begin")."\t".localtime(time)."\n";

if($para{'indexRefSeq'}==1){	
	print "bwaIndex...\n";
	if(! -e "ref/".$refSeqName.".bwt"){
		chdir("ref");
		system("bwa index -a bwtsw $refSeqName");
		print "samtoolsIndex...\n";
		system("samtools faidx $refSeqName");
		print "picardIndex...\n";
		system("java -Xmx30g -jar /usr/biosoft/picard/picard-tools-1.91/CreateSequenceDictionary.jar REFERENCE=".$refSeqName." OUTPUT=".substr($refSeqName,0,-3).".dict");
		chdir("..");
		print "Done!\n";
	}
	else{
		print "The indexed file has existed!\n";
	}
}
else{
	print "Skipped!\n";
}

#========================
#[Align the reads by bwa]
#========================

print "\n#========================\n";
print "#[Align the reads by bwa]\n";
print "#========================\n";
print LOG sprintf("%-60s","Align the reads by bwa...Begin")."\t".localtime(time)."\n";

if($para{'alignReads'}==1){
	if(! -e 'output'){system('mkdir output');}
	if(! -e 'output/bwaOut'){system('mkdir output/bwaOut');}

	print "Align the reads by bwa...\n";
	if($isPairedEnd){
		system("perl bin/runBwaPE.pl ".join('#',@samplesName)." $para{'numOfThreads'} $refSeqName");
	}
	else{
		system("perl bin/runBwaSE.pl ".join('#',@samplesName)." $para{'numOfThreads'} $refSeqName");
	}
	print "Extract alignment results...\n";
	if($isPairedEnd){
		system("perl bin/extractSamPE.pl ".join('#',@samplesName)." $para{'numOfThreads'}");
	}
	else{
		system("perl bin/extractSamSE.pl ".join('#',@samplesName)." $para{'numOfThreads'}");
	}
	print "Done!\n";
}
else{
	print "Skipped!\n";
}

#==================================
#[Convert initial bam to final bam]
#==================================

print "\n#==================================\n";
print "#[Convert initial bam to final bam]\n";
print "#==================================\n";
print LOG sprintf("%-60s","Convert initial bam to final bam...Begin")."\t".localtime(time)."\n";

if($para{'createFinalBam'}==1){
	if(! -e 'output/sam2BamTemp'){system('mkdir output/sam2BamTemp');}
	if(! -e 'output/finalBam'){system('mkdir output/finalBam');}

	print "Create final bam files...\n";
	system("perl bin/createFinalBam.pl ".join('#',@samplesName)." $para{'numOfThreads'} $refSeqName");

	print "Done!\n";
}
else{
	print "Skipped!\n";
}


#===================================
#[Call variants using multi-samples]
#===================================

print "\n#===================================\n";
print "#[Call variants using multi-samples]\n";
print "#===================================\n";
print LOG sprintf("%-60s","Call variants using multi-samples...Begin")."\t".localtime(time)."\n";

if($para{'CallVariantsMultiSamples'}==1){
	if(! -e 'output/geneticVariants'){system('mkdir output/geneticVariants');}

	print "call variants using multi-samples...\n";
	system("perl bin/callVariantsMultiSamples.pl $para{'numOfThreads'} $refSeqName");

	print "Done!\n";
}
else{
	print "Skipped!\n";
}

#===============================================
#[Filter variants for allSamples.gatk.final.vcf]
#===============================================

print "\n#===============================================\n";
print "#[Filter variants for allSamples.gatk.final.vcf]\n";
print "#===============================================\n";
print LOG sprintf("%-60s","Filter variants for allSamples.gatk.final.vcf...Begin")."\t".localtime(time)."\n";
if($para{'filterVariants'}==1){
	print sprintf("%-25s","MinAF").$para{'MinAF'}."\n";
	print sprintf("%-25s","MaxAF").$para{'MaxAF'}."\n";
	print sprintf("%-25s","missingDataFreq").$para{'missingDataFreq'}."\n";
	print sprintf("%-25s","partialHeterozygosity").$para{'partialHeterozygosity'}."\n\n";

	print "Filter variants for allSamples.gatk.final.vcf...\n";
	system("perl bin/filterVcf.gatk.pl $para{'MinAF'} $para{'MaxAF'} $para{'missingDataFreq'} $para{'partialHeterozygosity'}");

	print "Done!\n";
}
else{
	print "Skipped!\n";
}

print LOG sprintf("%-60s","Pipeline finished...End")."\t".localtime(time)."\n";



##############################################################################################

sub getSamplesName{
	my ($si);

	if($isPairedEnd){
		for($si=0;$si<@samples;$si+=2){
			$samplesName[$si/2]=substr($samples[$si],6,-8);
		}
	}
	else{
		for($si=0;$si<@samples;$si++){
			$samplesName[$si]=substr($samples[$si],6,-8);
		}
	}
}
