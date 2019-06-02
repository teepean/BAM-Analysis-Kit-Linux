#!/bin/bash

/bin/echo -ne "\033]2;BAM Analysis Kit 2.6 Linux\007"

#     The MIT License (MIT)
#     Copyright © 2013-2015 Felix Immanuel
#     http://www.y-str.org
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy
#     of this software and associated documentation files (the Software), to deal
#     in the Software without restriction, including without limitation the rights
#     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#     copies of the Software, and to permit persons to whom the Software is furnished
#     to do so, subject to the following conditions: The above copyright notice and
#     this permission notice shall be included in all copies or substantial portions
#     of the Software. THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY
#     KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
#     EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
#     OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
#     ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#     OTHER DEALINGS IN THE SOFTWARE.

/bin/echo
/bin/echo "*** BAM Analysis Kit 2.6 ***"
/bin/echo
/bin/echo "Project Website: http://www.y-str.org/2014/04/bam-analysis-kit.html"
/bin/echo "Tools Used: SAMTools, picard, bamtools, Cygwin, GATK, Java, Yleaf v2, haplogrep"
/bin/echo "Script Developer: Felix Immanuel ^<i@fi.id.au^>"
/bin/echo

if [ $# = 0 ]; then
	/bin/echo
	/bin/echo "Syntax:"
	/bin/echo "   console_bam.sh <bam-file>"
	/bin/echo
	read key
	exit 1
fi

# - start reporting versions..
bedtools --version
samtools --version
/bin/echo "Kernel Version: $(uname -r)"
java -version 2>&1 |grep -i version
/bin/echo -n 'Picard Version: '; java -jar bin/picard/picard.jar AddOrReplaceReadGroups --version
# - end reporting versions..

# - saving old processing if it wasnt saved and accidentally started new processing
if [ -f out.old/genome_complete.txt ]; then
	rm -fr out.old
fi

if [ -f out/genome_complete.txt ]; then
	mv out out.old
	mkdir out
fi

# source bamkit.config

# These parameters define what to process
#
CHR_1=no
CHR_2=no
CHR_3=no
CHR_4=no
CHR_5=no
CHR_6=no
CHR_7=no
CHR_8=no
CHR_9=no
CHR_10=no
CHR_11=no
CHR_12=no
CHR_13=no
CHR_14=no
CHR_15=no
CHR_16=no
CHR_17=no
CHR_18=no
CHR_19=no
CHR_20=no
CHR_21=no
CHR_22=no
CHR_X=no
CHR_Y=no
CHR_M=no
#
# Execution Configuration
BAMKIT_JVM=-Xmx2048m
BAMKIT_THREADS=8
DEL_VCF=no


/bin/echo
/bin/echo "Configuration from bamkit.config"
/bin/echo
/bin/echo "Process chr 1 set to $CHR_1"

/bin/echo "Process chr 2 set to $CHR_2"
/bin/echo "Process chr 3 set to $CHR_3"
/bin/echo "Process chr 4 set to $CHR_4"
/bin/echo "Process chr 5 set to $CHR_5"
/bin/echo "Process chr 6 set to $CHR_6"
/bin/echo "Process chr 7 set to $CHR_7"
/bin/echo "Process chr 8 set to $CHR_8"
/bin/echo "Process chr 9 set to $CHR_9"
/bin/echo "Process chr 10 set to $CHR_10"
/bin/echo "Process chr 11 set to $CHR_11"
/bin/echo "Process chr 12 set to $CHR_12"
/bin/echo "Process chr 13 set to $CHR_13"
/bin/echo "Process chr 14 set to $CHR_14"
/bin/echo "Process chr 15 set to $CHR_15"
/bin/echo "Process chr 16 set to $CHR_16"
/bin/echo "Process chr 17 set to $CHR_17"
/bin/echo "Process chr 18 set to $CHR_18"
/bin/echo "Process chr 19 set to $CHR_19"
/bin/echo "Process chr 20 set to $CHR_20"
/bin/echo "Process chr 21 set to $CHR_21"
/bin/echo "Process chr 22 set to $CHR_22"
/bin/echo "Process chr X set to $CHR_X"
/bin/echo "Process chr Y set to $CHR_Y"
/bin/echo "Process chr M set to $CHR_M"
/bin/echo "JVM value is $BAMKIT_JVM"
/bin/echo "No of parallel threads is $BAMKIT_THREADS"
/bin/echo "Delete VCF after processing set to $DEL_VCF"

/bin/echo
/bin/echo
/bin/echo "Input BAM : $1"
/bin/echo

/bin/echo "Pre-execution Cleanup ..."

rm -f chrY_1.tab \
	chrY_1.tmp \
	chrM_2.tmp \
	mtdna_max_pos \
	dbsnp_chrM.tab \
	snps.sorted \
	ref.sorted \
	snps.tmp \
	chrY.tmp \
	inchr.bam \
	bam_wh_tmp.bam \
	ref.fa \
	ref.fa.fai \
	ref.dict \
	bam_sorted.bam \
	bam_sorted.bam.bai \
	bam_sorted_realigned.bam \
	header \
	header01 \
	header02 \
	inchr.sam \
	tmp.sam \
	chr*.bam \
	reads.bam \
	bam.intervals \
	bam_sorted_realigned.bam.bai \
	bam_sorted_realigned.bai \
	bam_wh.bam \
	bam_out.vcf \
	bam_out.vcf.idx \
	chr \
	bam_chr*.vcf \
	bam_chr*.vcf.gz \
	bam_chr*.vcf.gz.tbi \
	bam_chr*.vcf.idx \
	chr*.tab \
	chrM.seq \
	chrM.tmp.tab \
	chrM.tmp \
	genotype.txt \
	snps.txt \
	bam_complete_sorted.bam \
	bam_complete_sorted.bam.bai \
	ystr.filters \
	bed.a \
	bam_strs.aligned.stats \
	bam_strs_sorted.bam \
	bam_strs_sorted.bam.bai \
	bam_ystrs.allelotype.stats \
	bam_ystrs.vcf \
	bam_strs.aligned.bam \
	bam_out_variants.vcf
/bin/echo

/bin/echo "Sorting ..."

samtools sort -@ $BAMKIT_THREADS "$1" -o bam_complete_sorted.bam

/bin/echo
/bin/echo "Indexing the sorted BAM file ..."
samtools index -@ $BAMKIT_THREADS bam_complete_sorted.bam

/bin/echo -e '# rsid\tchr\tpos\tgenotype' > out/genome_full_snps.txt
/bin/echo -e '# chr\tpos\tgenotype\trsid' > out/genome_complete.txt
/bin/echo -e 'RSID,CHROMOSOME,POSITION,RESULT' > out/filtered-autosomal-o37-results.csv
/bin/echo -e 'RSID,CHROMOSOME,POSITION,RESULT' > out/filtered-x-chromosome-o37-results.csv

# --- Execute genome tools
for A in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M ; do
  mulkku=CHR_$A
  
	if [[ "${!mulkku}" == yes ]]; then

		/bin/echo
		/bin/echo "******* Processing Chromosome $A *******"
		/bin/echo

		/bin/echo "Splitting and preparing Chr $A ..."
		# it can be chr22 or 22 -- must be a way to detect.
		samtools view -H bam_complete_sorted.bam |cut -f2 |/bin/grep SN |/bin/grep $A |cut -d: -f2 |head -1 > chr
		for C in $(/usr/bin/gawk 'NF {print $1}' chr); do
			samtools view -@ $BAMKIT_THREADS -b bam_complete_sorted.bam $C > chr$A.bam
		done
		
		cp ref/chr$A.fa.gz ref.fa.gz >/dev/null

		gzip -d ref.fa.gz
        
        if [ "$A" = "M" ]; then
		if test -f "chrchrM.bam"; then
			mv chrchrM.bam chrM.bam > /dev/null
		fi
		if test -f "chrMT.bam"; then
			mv chrMT.bam chrM.bam > /dev/null
		fi
		fi
		
		/bin/echo "Checking and fixing ..."
		rm -f inchr.bam
		cp chr$A.bam inchr.bam >/dev/null

		/bin/echo
		/bin/echo "Indexing Human Reference Genome ..."
		samtools faidx ref.fa
		/bin/echo
		/bin/echo "Creating Sequence Dictionary for Human Reference Genome ..."
		java $BAMKIT_JVM -jar bin/picard/picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict

		/bin/echo
		/bin/echo "Preparing BAM ..."
		samtools view  -@ $BAMKIT_THREADS inchr.bam |/bin/sed 's/\t/\tchr/2' > tmp.sam
		/bin/cat tmp.sam |/bin/sed 's/\tchrchr/\tchr/' > inchr.sam

		if [ "$A" = "M" ]; then
		echo "kopiooraan aamm"
			/bin/cat inchr.sam |/bin/sed 's/\tchrMT/\tchrM/' > inchr_tmp.sam
			rm -f inchr.sam
			cp inchr_tmp.sam inchr.sam > /dev/null
			rm -f inchr_tmp.sam
		fi

		samtools view -@ $BAMKIT_THREADS -bT ref.fa inchr.sam > reads.bam
		samtools view -@ $BAMKIT_THREADS -H inchr.bam |/bin/grep -v SN > header01
		samtools view -@ $BAMKIT_THREADS -H reads.bam > header02
		/bin/cat header01 header02 > header
		samtools reheader header reads.bam > bam_wh_tmp.bam

		/bin/echo "Adding or Replace Read Group Header ..."
		java $BAMKIT_JVM -jar bin/picard/picard.jar AddOrReplaceReadGroups INPUT=bam_wh_tmp.bam OUTPUT=bam_wh.bam SORT_ORDER=coordinate RGID=rgid RGLB=rglib RGPL=illumina RGPU=rgpu RGSM=sample VALIDATION_STRINGENCY=SILENT

		/bin/echo "Sorting ..."
		samtools sort -@ $BAMKIT_THREADS bam_wh.bam -o bam_sorted.bam

		/bin/echo
		/bin/echo "Indexing the sorted BAM file ..."
		samtools index -@ $BAMKIT_THREADS bam_sorted.bam

		/bin/echo
		/bin/echo "Realignment of the sorted and indexed BAM file ..."
		if [ "$A" = "M" ]; then
			java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator   -window 3  -minReads 1 -R ref.fa -I bam_sorted.bam  -o bam.intervals
			java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -T IndelRealigner  -maxPosMove 10 -R ref.fa -I bam_sorted.bam -targetIntervals bam.intervals -o bam_sorted_realigned.bam
		else
			java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref.fa -I bam_sorted.bam  -o bam.intervals
			java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R ref.fa -I bam_sorted.bam -targetIntervals bam.intervals -o bam_sorted_realigned.bam
		fi

		/bin/echo
		/bin/echo "Indexing the realigned BAM file ..."
		samtools index -@ $BAMKIT_THREADS bam_sorted_realigned.bam

		/bin/echo
		/bin/echo "Invoke the variant caller ..."

		unset HAPLOID
		[ "$A" = "M" -o "$A" = "Y" ] && HAPLOID=1
		
		if [ -n "$HAPLOID" ]; then
		
			if [ "$A" = "M" ]; then
				
				java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -l INFO -R ref.fa -T HaplotypeCaller -I bam_sorted_realigned.bam -nct $BAMKIT_THREADS -o bam_chr$A.vcf
				bgzip bam_chr$A.vcf
				tabix bam_chr$A.vcf.gz

				echo "Extracting mtDNA markers .."
				bcftools query -f '%CHROM\t%POS\t[%IUPACGT]\n' bam_chr$A.vcf.gz |sed 's/chr//g' |sed 's/\///g' > chr$A.tab
				bcftools query -f '%CHROM\t%POS\t[%IUPACGT]\n' bam_chr$A.vcf.gz |sed 's/chr//g' |sed 's/\///g' |cut -f2,3 |/usr/bin/gawk '{if($2 ~ /^[ATGC]$/) printf $1  $2 " "}' |sed 's/\s$//g' |sed 's/\s/, /g' > out/rCRS_mtDNA.txt

				# - realign with RSRS to get accurate mtDNA markers in RSRS

				rm -f \
					bam_chr$A.* \
					ref.fa \
					ref.dict \
					ref.fa.fai
				cp ref/chrM.RSRS.fasta ref.fa
				
				samtools faidx ref.fa
				java $BAMKIT_JVM -jar bin/picard/picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict
				rm -f bam.intervals
				java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator   -window 3  -minReads 1 -R ref.fa -I bam_sorted.bam  -o bam.intervals
				rm -f bam_sorted_realigned.bam
				java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -T IndelRealigner  -maxPosMove 10 -R ref.fa -I bam_sorted.bam -targetIntervals bam.intervals -o bam_sorted_realigned.bam
				java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -l INFO -R ref.fa -T HaplotypeCaller -I bam_sorted_realigned.bam -nct $BAMKIT_THREADS -o bam_chr$A.vcf
				bgzip bam_chr$A.vcf
				tabix bam_chr$A.vcf.gz
				java -jar bin/haplogrep/haplogrep-2.1.20.jar --in bam_chrM.vcf.gz --format vcf --out out/mtDNA-haplogroup-haplogrep.txt
				bcftools query -f '%CHROM\t%POS\t[%IUPACGT]\n' bam_chr$A.vcf.gz |sed 's/chr//g' |sed 's/\///g' > chr$A.tab	
				bcftools query -f '%REF\t%POS\t[%IUPACGT]\n' bam_chr$A.vcf.gz |sed 's/chr//g' |sed 's/\///g' |/usr/bin/gawk '{if($3 ~ /^[ATGC]$/) printf $1 $2 $3 " "}' |sed 's/\s$//g' |sed 's/\s/, /g' > out/RSRS_mtDNA.txt

			fi
		
			if [ "$A" = "Y" ]; then
				java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -l INFO -R ref.fa -T UnifiedGenotyper -glm SNP -I bam_sorted_realigned.bam -rf BadCigar -nct $BAMKIT_THREADS -o bam_chr$A.vcf --output_mode EMIT_ALL_CONFIDENT_SITES
				bgzip bam_chr$A.vcf
				tabix bam_chr$A.vcf.gz
				echo "Extracting Y-SNP markers .."
				bcftools query -f '%CHROM\t%POS\t[%IUPACGT]\n' bam_chr$A.vcf.gz | sed 's/chr//g' > chr$A.tab
				python3 bin/yleaf/Yleaf.py -bam bam_sorted_realigned.bam -ref hg19 -out out -q 10 -b 40 -t 8 -r 1
				mv out/out.hg out/y-haplogroup.txt
				mv out/bam_sorted_realigned/bam_sorted_realigned.chr out
				mv out/bam_sorted_realigned/bam_sorted_realigned.fmf out
				mv out/bam_sorted_realigned/bam_sorted_realigned.log out
				mv out/bam_sorted_realigned/bam_sorted_realigned.out out
				rm -rf out/bam_sorted_realigned
			fi
		else
			# not defined haploids
			java $BAMKIT_JVM -jar bin/gatk/GenomeAnalysisTK.jar -l INFO -R ref.fa -T UnifiedGenotyper -glm SNP -I bam_sorted_realigned.bam -rf BadCigar -nct $BAMKIT_THREADS -o bam_chr$A.vcf --output_mode EMIT_ALL_CONFIDENT_SITES
			bgzip bam_chr$A.vcf
			tabix bam_chr$A.vcf.gz
			bcftools query -f '%CHROM\t%POS\t[%TGT]\n' bam_chr$A.vcf.gz |/bin/sed 's/chr//g' > chr$A.tab
		fi

		/usr/bin/gawk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{ print a[$1,$2],$1,$2,$3}' <(/bin/gzip -dc ref/dbsnp_chr$A.tab.gz) chr$A.tab|/bin/sed 's/\\s/\t/g' > snps.tmp
		/bin/cat snps.tmp | /bin/sed 's/\///g' >> out/genome_full_snps.txt
		/usr/bin/gawk 'NR==FNR{a[$1,$2]=$3;next} { print $1,$2,$3,a[$1,$2]}' <(/bin/gzip -dc ref/dbsnp_chr$A.tab.gz) chr$A.tab|/bin/sed 's/\\s/\t/g' >> out/genome_complete.txt

		/bin/cat snps.tmp |sort -k 1 > snps.sorted

		if [ "$A" = "X" ]; then
			/bin/cat ref/snps_filtered/chr$A |sort -k 1 > ref.sorted
			join snps.sorted ref.sorted |/bin/sed 's/\///g' |sort -n -k 3 |/usr/bin/gawk '{print "\042" $1 "\042," "\042" $2 "\042," "\042" $3 "\042," "\042" $4 "\042"}' >> out/filtered-x-chromosome-o37-results.csv
		else
			if [ "$A" != "Y" -a "$A" != "M" ]; then
				/bin/cat ref/snps_filtered/chr$A |sort -k 1 > ref.sorted
				join snps.sorted ref.sorted |/bin/sed 's/\///g' |sort -n -k 3 |/usr/bin/gawk '{print "\042" $1 "\042," "\042" $2 "\042," "\042" $3 "\042," "\042" $4 "\042"}' >> out/filtered-autosomal-o37-results.csv
			fi
		fi

		# -- final cleanup
		rm -f \
			chrY_1.tab  \
			chrY_1.tmp  \
			chrM_2.tmp  \
			mtdna_max_pos  \
			dbsnp_chrM.tab \
			snps.sorted \
			ref.sorted \
			snps.tmp \
			chrY.tmp \
			lobSTR_CODIS.out \
			lobSTR_Y-STR.out \
			inchr.bam \
			bam_wh_tmp.bam \
			ref.fa \
			ref.fa.fai \
			ref.dict \
			bam_sorted.bam \
			bam_sorted.bam.bai \
			bam_sorted_realigned.bam \
			header \
			header01 \
			header02 \
			inchr.sam \
			tmp.sam \
			chr$A.bam \
			reads.bam \
			bam.intervals \
			bam_sorted_realigned.bam.bai \
			bam_sorted_realigned.bai \
			bam_wh.bam \
			bam_out.vcf \
			bam_out.vcf.idx \
			chr \
			bam_chr$A.vcf  \
			bam_chr$A.vcf.gz.tbi \
			bam_chr$A.vcf.idx \
			chr$A.tab \
			chrM.seq \
			chrM.tmp.tab \
			chrM.tmp
			
		if [ "$DEL_VCF" = "yes" ]; then
			rm -f bam_chr$A.vcf.gz
		else
			[ -f "bam_chr$A.vcf.gz" ] && mv bam_chr$A.vcf.gz out
		fi
	fi
done

rm -f \
	genotype.txt \
	snps.txt \
	bam_complete_sorted.bam \
	bam_complete_sorted.bam.bai \
	ystr.filters \
	bed.a \
	bam_strs.aligned.stats \
	bam_strs_sorted.bam \
	bam_strs_sorted.bam.bai \
	bam_ystrs.allelotype.stats \
	bam_ystrs.vcf \
	bam_strs.aligned.bam \
	bam_out_variants.vcf

# xdg-open out
/bin/echo
/bin/echo "All Tasks Completed. Please find results in out subfolder."
/bin/echo "Also check the logs/info in this window for errors (if any)."
# read key
