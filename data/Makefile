#################################################
#########   HUMAN DATA DOWNLOAD   ###############
#################################################
#################################################
################   SOURCES  #####################
#################################################
#####    ftp://ftp.1000genomes.ebi.ac.uk     ####
#####    ftp://ftp.sra.ebi.ac.uk             ####
#####    ftp://ftp-trace.ncbi.nih.gov        ####
#################################################
#################################################
# use following command to download all data    #
# used for our exmeriments                      #
#                                               #
# make                                          #
#                                               #
#################################################

all_for_experiments: hs37d5.fa human.vcf reads.fq reads2.fq


clean:
	rm -f hs37d5.fa human.vcf reads.fq reads2.fq all_vcf_downloaded hs37d5.fa.fai
	rm -f -r gz

gz_dir:
	mkdir -p gz

all_vcf_downloaded: gz_dir
	#DOWNLOAD CHROM 1 - 22 + X + Y
	for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do \
	wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$$i.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -O gz/chrom$$i.vcf.gz &  \
	done; \
	wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz -O gz/chromX.vcf.gz & \
	wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz -O gz/chromY.vcf.gz & \
	echo waiting for download... ; \
	wait
	touch all_vcf_downloaded

human.vcf: all_vcf_downloaded
	# UNZIP AND CONCATENATE VCF FILES
	# -dkc: d = decompress, k = keep original, c = write to stdio
	echo unzip chrom1
	gzip -dkc gz/chrom1.vcf.gz | awk -v OFS='\t' '{print $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8}' >  human.vcf &
	for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do \
	echo unzip chrom$$i; \
	gzip -dkc gz/chrom$$i.vcf.gz | awk -v OFS='\t' '{print $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8}' | grep -v '##' >  tmp$$i.vcf & \
	done; \
	wait
	for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do \
	cat tmp$$i.vcf >> human.vcf; rm tmp$$i.vcf; done


hs37d5.fa: gz_dir
	#GET hs37d5.fa FROM 1000 GENOMES PROJEKT PHASE 2
	curl ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz -o gz/hs37d5.fa.gz
	#awk '{print $1}': only print '>1' not '>1 dna:chromosome chromosome:GRCh37:1:1:249250621:1' as sequence name. ONLY CHROM 1 - 22 + X + Y
	gunzip -d -c gz/hs37d5.fa.gz | awk  'NR==1, NR==51594656 {print $$1}' > hs37d5.fa
	#rm -f hs37d5.fa.gz

################ FOR EXPERIMENTS WITH ONLY CHROM 1 OR 22 #########################
chrom1 : chrom1.fa  chrom1.vcf
chrom22: chrom22.fa chrom22.vcf

chrom22.fa: hs37d5.fa
	#CUT CHROMOSOM 22 FROM hs37d5
	awk 'NR==47162177, NR==48017254 {print $$0}'  hs37d5.fa > chrom22.fa

gz/chrom22.vcf.gz: gz_dir
	#GET DATA FROM 1000 GENOMES PROJEKT PHASE 3
	curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -o gz/chrom22.vcf.gz

chrom22.vcf: gz/chrom22.vcf.gz
	gzip -dkc gz/chrom22.vcf.gz | awk -v OFS='\t' '{print $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8}' >>  chrom22.vcf

chrom1.fa: hs37d5.fa
	#CUT CHROMOSOM 1 FROM hs37d5
	awk 'NR==1, NR==4154179 {print $$0}'  hs37d5.fa > chrom1.fa

gz/chrom1.vcf.gz: gz_dir
	#GET DATA FROM 1000 GENOMES PROJEKT PHASE 3
	curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -o gz/chrom1.vcf.gz

chrom1.vcf: gz/chrom1.vcf.gz
	gzip -dkc gz/chrom1.vcf.gz | awk -v OFS='\t' '{print $$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8}' >>  chrom1.vcf


################################ REAL READS ##################################
reads.fq:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR769/SRR769545/SRR769545_2.fastq.gz
	gzip -d SRR769545_2.fastq.gz
	head -n 40000000 SRR769545_2.fastq > reads.fq
	rm SRR769545_2.fastq

reads2.fq:
	wget ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
	gzip -d NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
	head -n 40000000 NIST7035_TAAGGCGA_L001_R1_001.fastq > reads2.fq
	rm NIST7035_TAAGGCGA_L001_R1_001.fastq
