include ../PATHS
#################################################
################# EXPERIMENT ####################
#################################################
# use following commands to run experiments of  #
# the paper (please ensure you downloaded the   #
# data before)                                  #
#                                               #
# make                                          #
#                                               #
#################################################
GEDMAP_DIR  = gedmap
HISAT_DIR   = hisat2
GIRAFFE_DIR = vg_giraffe
#################################################

all_expermiments:   $(GEDMAP_DIR)/index.out samples/all $(GEDMAP_DIR)/eval.out $(HISAT_DIR)/index.out  $(HISAT_DIR)/eval.out  $(GIRAFFE_DIR)/index.out $(GIRAFFE_DIR)/eval.out  compare_sam_ged_hisat2
clean: $(GEDMAP_DIR)/clean $(HISAT_DIR)/clean $(GIRAFFE_DIR)/clean samples/clean
	rm -f compare_sam_ged_hisat2

################ INPUT FILES  ###################

DATA   = ../data
GZ_DIR = ../data/gz

FA   = $(DATA)/hs37d5.fa
VCF  = $(DATA)/human.vcf
READS = $(DATA)/reads.fq
READS2 = $(DATA)/reads2.fq
VCF_GZ  = $(GZ_DIR)/chrom1.vcf.gz \
-v $(GZ_DIR)/chrom2.vcf.gz \
-v $(GZ_DIR)/chrom3.vcf.gz \
-v $(GZ_DIR)/chrom4.vcf.gz \
-v $(GZ_DIR)/chrom5.vcf.gz \
-v $(GZ_DIR)/chrom6.vcf.gz \
-v $(GZ_DIR)/chrom7.vcf.gz \
-v $(GZ_DIR)/chrom8.vcf.gz \
-v $(GZ_DIR)/chrom10.vcf.gz \
-v $(GZ_DIR)/chrom11.vcf.gz \
-v $(GZ_DIR)/chrom12.vcf.gz \
-v $(GZ_DIR)/chrom13.vcf.gz \
-v $(GZ_DIR)/chrom14.vcf.gz \
-v $(GZ_DIR)/chrom15.vcf.gz \
-v $(GZ_DIR)/chrom16.vcf.gz \
-v $(GZ_DIR)/chrom17.vcf.gz \
-v $(GZ_DIR)/chrom18.vcf.gz \
-v $(GZ_DIR)/chrom19.vcf.gz \
-v $(GZ_DIR)/chrom20.vcf.gz \
-v $(GZ_DIR)/chrom21.vcf.gz \
-v $(GZ_DIR)/chrom22.vcf.gz \
-v $(GZ_DIR)/chromX.vcf.gz \
-v $(GZ_DIR)/chromY.vcf.gz


#################  PARAMETERS ####################

GEDMAP_K = 15
GEDMAP_W = 7

GEDMAP_PARAMS = -rc -tc 64
GEDMAP_PARAMS_SAM0 =  $(GEDMAP_PARAMS) -mac 50 -mat 100
GEDMAP_PARAMS_SAM1 =  $(GEDMAP_PARAMS) -mac 50 -mat 100
GEDMAP_PARAMS_SAM3 =  $(GEDMAP_PARAMS) -mac 50 -mat 100
GEDMAP_PARAMS_SAM8 =  $(GEDMAP_PARAMS) -d 20
GEDMAP_PARAMS_SAM9 =  $(GEDMAP_PARAMS) -d 20
GEDMAP_PARAMS_SAM11 = $(GEDMAP_PARAMS)
GEDMAP_PARAMS_READ =  $(GEDMAP_PARAMS)

HISAT_PARAMS_SAM0 = -p 64 --ignore-quals --score-min L,0,-0.6 --no-softclip --mp 1,1 --rdg 0,1 --rfg 0,1 --np 0 -k 1 -t
HISAT_PARAMS_SAM1 = $(HISAT_PARAMS_SAM0)
HISAT_PARAMS_SAM2 = $(HISAT_PARAMS_SAM0)
HISAT_PARAMS_SAM3 = $(HISAT_PARAMS_SAM0) --score-min L,0,-0.8
HISAT_PARAMS_SAM8 = $(HISAT_PARAMS_SAM0)
HISAT_PARAMS_SAM9 = $(HISAT_PARAMS_SAM0)
HISAT_PARAMS_SAM10 = $(HISAT_PARAMS_SAM0)
HISAT_PARAMS_SAM11 = $(HISAT_PARAMS_SAM0) --score-min L,0,-0.8

GIRAFFE_PARAMS = -t 64 -o SAM

################### SAMPLES #####################
SAMPLE_DIR = samples

SAMPLE0 = 100_none
SAMPLE1 = 100_low
SAMPLE3 = 100_high

SAMPLE8 = 500_none
SAMPLE9 = 500_low
SAMPLE11 = 500_high

SAMPLE_OE = 100_low_edge

SAMPLES = $(SAMPLE0) $(SAMPLE8) $(SAMPLE1) $(SAMPLE9) $(SAMPLE3) $(SAMPLE11) $(SAMPLE_OE)

common_sample_params = $(GEDS) -c 10000000 -2fa $(2FA) -a $(ADJ)
zero_error = -e 0   -e_d 0    -e_i 0
low_error  = -e 100 -e_d 1000 -e_i 1000
mid_error  = -e 50  -e_d 500  -e_i 500
high_error = -e 25  -e_d 250  -e_i 250

$(SAMPLE_DIR)/$(SAMPLE0).fq:  $(GEDS)
	mkdir -p $(SAMPLE_DIR); $(gedmap) sample $(common_sample_params) -o $(SAMPLE_DIR)/$(SAMPLE0).fq -l 100 $(zero_error)
$(SAMPLE_DIR)/$(SAMPLE1).fq:  $(GEDS)
	mkdir -p $(SAMPLE_DIR); $(gedmap) sample $(common_sample_params) -o $(SAMPLE_DIR)/$(SAMPLE1).fq -l 100 $(low_error)
$(SAMPLE_DIR)/$(SAMPLE3).fq:  $(GEDS)
	mkdir -p $(SAMPLE_DIR); $(gedmap) sample $(common_sample_params) -o $(SAMPLE_DIR)/$(SAMPLE3).fq -l 100 $(high_error)
$(SAMPLE_DIR)/$(SAMPLE8).fq:  $(GEDS)
	mkdir -p $(SAMPLE_DIR); $(gedmap) sample $(common_sample_params) -o $(SAMPLE_DIR)/$(SAMPLE8).fq  -l 500 $(zero_error)
$(SAMPLE_DIR)/$(SAMPLE9).fq:  $(GEDS)
	mkdir -p $(SAMPLE_DIR); $(gedmap) sample $(common_sample_params) -o $(SAMPLE_DIR)/$(SAMPLE9).fq  -l 500 $(low_error)
$(SAMPLE_DIR)/$(SAMPLE11).fq:  $(GEDS)
	mkdir -p $(SAMPLE_DIR); $(gedmap) sample $(common_sample_params) -o $(SAMPLE_DIR)/$(SAMPLE11).fq -l 500 $(high_error)
$(SAMPLE_DIR)/$(SAMPLE_OE).fq:  $(GEDS)
	mkdir -p $(SAMPLE_DIR); $(gedmap) sample $(common_sample_params) -ae $(ADJ) -o $(SAMPLE_DIR)/$(SAMPLE_OE).fq -l 100 $(low_error)

samples/all: $(GEDS)
	for target in $(SAMPLES); do make $(SAMPLE_DIR)/$${target}.fq ; done

samples/clean:
	rm -r -f $(SAMPLE_DIR)
##################### GEDMAP ######################
gedmap     = ../gedmap

$(gedmap):
	cd .. && make gedmap

GEDS = $(GEDMAP_DIR)/human.geds
MINI = $(GEDMAP_DIR)/human.$(GEDMAP_W).$(GEDMAP_K).min
ADJ = $(GEDS).adj
2FA = $(GEDS).2fa

$(GEDS): $(FA) $(VCF)
	mkdir -p $(GEDMAP_DIR)
	$(gedmap) parse $(FA) $(VCF) $(GEDS)

$(MINI): $(GEDS)
	mkdir -p $(GEDMAP_DIR)
	$(gedmap) index  $(GEDS)  -k $(GEDMAP_K) -w $(GEDMAP_W) -o $(MINI) -a $(ADJ) -2fa $(2FA) -trim 1000

common_align_params = $(MINI) -2fa $(2FA) -a $(ADJ)

$(GEDMAP_DIR)/$(SAMPLE0).sam: $(MINI) $(SAMPLE_DIR)/$(SAMPLE0).fq
	$(gedmap) align  $(GEDS) $(SAMPLE_DIR)/$(SAMPLE0).fq  $(common_align_params) -o  $(GEDMAP_DIR)/$(SAMPLE0).sam   $(GEDMAP_PARAMS_SAM0)
$(GEDMAP_DIR)/$(SAMPLE1).sam: $(MINI) $(SAMPLE_DIR)/$(SAMPLE1).fq
	$(gedmap) align  $(GEDS) $(SAMPLE_DIR)/$(SAMPLE1).fq  $(common_align_params) -o  $(GEDMAP_DIR)/$(SAMPLE1).sam   $(GEDMAP_PARAMS_SAM1)
$(GEDMAP_DIR)/$(SAMPLE3).sam: $(MINI) $(SAMPLE_DIR)/$(SAMPLE3).fq
	$(gedmap) align  $(GEDS) $(SAMPLE_DIR)/$(SAMPLE3).fq  $(common_align_params) -o  $(GEDMAP_DIR)/$(SAMPLE3).sam   $(GEDMAP_PARAMS_SAM3)
$(GEDMAP_DIR)/$(SAMPLE8).sam: $(MINI) $(SAMPLE_DIR)/$(SAMPLE8).fq
	$(gedmap) align  $(GEDS) $(SAMPLE_DIR)/$(SAMPLE8).fq  $(common_align_params) -o  $(GEDMAP_DIR)/$(SAMPLE8).sam   $(GEDMAP_PARAMS_SAM8)
$(GEDMAP_DIR)/$(SAMPLE9).sam: $(MINI) $(SAMPLE_DIR)/$(SAMPLE9).fq
	$(gedmap) align  $(GEDS) $(SAMPLE_DIR)/$(SAMPLE9).fq  $(common_align_params) -o  $(GEDMAP_DIR)/$(SAMPLE9).sam   $(GEDMAP_PARAMS_SAM9)
$(GEDMAP_DIR)/$(SAMPLE11).sam: $(MINI) $(SAMPLE_DIR)/$(SAMPLE11).fq
	$(gedmap) align  $(GEDS) $(SAMPLE_DIR)/$(SAMPLE11).fq  $(common_align_params) -o  $(GEDMAP_DIR)/$(SAMPLE11).sam   $(GEDMAP_PARAMS_SAM11)
$(GEDMAP_DIR)/$(SAMPLE_OE).sam: $(MINI) $(SAMPLE_DIR)/$(SAMPLE_OE).fq
	$(gedmap) align  $(GEDS) $(SAMPLE_DIR)/$(SAMPLE_OE).fq  $(common_align_params) -o  $(GEDMAP_DIR)/$(SAMPLE_OE).sam   $(GEDMAP_PARAMS_SAM1)
$(GEDMAP_DIR)/reads.sam: $(MINI) $(READS)
	$(gedmap) align  $(GEDS) $(READS)  $(common_align_params) -o  $(GEDMAP_DIR)/reads.sam   $(GEDMAP_PARAMS_READ)
$(GEDMAP_DIR)/reads2.sam: $(MINI) $(READS2)
	$(gedmap) align  $(GEDS) $(READS2)  $(common_align_params) -o  $(GEDMAP_DIR)/reads2.sam   $(GEDMAP_PARAMS_READ)
#### in order to compare alignment-wise #####
$(GEDMAP_DIR)/reads.in-order.sam: $(MINI) $(READS)
	$(gedmap) align  $(GEDS) $(READS)  $(common_align_params) -o  $(GEDMAP_DIR)/reads.in-order.sam   $(GEDMAP_PARAMS_READ) -io
$(GEDMAP_DIR)/reads2.in-order.sam: $(MINI) $(READS2)
	$(gedmap) align  $(GEDS) $(READS2)  $(common_align_params) -o  $(GEDMAP_DIR)/reads2.in-order.sam   $(GEDMAP_PARAMS_READ) -io

################## EVALUATE INDEX WITH DIFFERENT PARAMETER
IDX_EVAL_DIR  = $(GEDMAP_DIR)/index_eval


eval_index_parameter: $(IDX_EVAL_DIR)/samples.1000000.fq $(IDX_EVAL_DIR)/samples.10000.fq $(IDX_EVAL_DIR)/samples.100.fq $(IDX_EVAL_DIR)/indexes $(IDX_EVAL_DIR)/eval_idx $(IDX_EVAL_DIR)/eval_map $(IDX_EVAL_DIR)/eval_qual

$(IDX_EVAL_DIR)/indexes:
	mkdir -p $(IDX_EVAL_DIR)/indexes
	rm $(IDX_EVAL_DIR)/indexes/index.out -f
	for t in 0 1000; do \
		for w in 4 7 15 30; do \
			for k in 5 10 13 15 16 17; do \
				$(gedmap) index  $(GEDS)  -k $$k -w $$w  -o $(IDX_EVAL_DIR)/indexes/human.t$$t.$$w.$$k.min -a $(ADJ) -2fa $(2FA) -trim $$t  | tee -a $(IDX_EVAL_DIR)/indexes/index.out ; \
			done; \
		done; \
	done

$(IDX_EVAL_DIR)/samples.1000000.fq:
	mkdir -p $(IDX_EVAL_DIR)/
	$(gedmap) sample $(GEDS) -c 1000000 -o $(IDX_EVAL_DIR)/samples.1000000.fq -l 100 -e 0 -e_d 0 -e_i 0 -rc 0


$(IDX_EVAL_DIR)/samples.10000.fq: $(IDX_EVAL_DIR)/samples.1000000.fq
	head -n 40000 $(IDX_EVAL_DIR)/samples.1000000.fq > $(IDX_EVAL_DIR)/samples.10000.fq


$(IDX_EVAL_DIR)/samples.100.fq: $(IDX_EVAL_DIR)/samples.1000000.fq
	head -n 400 $(IDX_EVAL_DIR)/samples.1000000.fq > $(IDX_EVAL_DIR)/samples.100.fq


$(IDX_EVAL_DIR)/eval_idx:
	mkdir -p $(IDX_EVAL_DIR)/eval
	for t in 0 1000; do \
		for w in 4 7 15 30; do \
			../eval_idx $(IDX_EVAL_DIR)/indexes/human.t1000.30.15.min $(IDX_EVAL_DIR)/samples.100.fq | head -n 1 > $(IDX_EVAL_DIR)/eval/eval_idx.t$$t.w$$w.dat; \
			for k in 5; do \
				../eval_idx $(IDX_EVAL_DIR)/indexes/human.t$$t.$$w.$$k.min $(IDX_EVAL_DIR)/samples.100.fq | tail -n 1 >> $(IDX_EVAL_DIR)/eval/eval_idx.t$$t.w$$w.dat;  \
			done; \
			for k in 10 13 15 16 17; do \
				../eval_idx $(IDX_EVAL_DIR)/indexes/human.t$$t.$$w.$$k.min $(IDX_EVAL_DIR)/samples.10000.fq | tail -n 1 >> $(IDX_EVAL_DIR)/eval/eval_idx.t$$t.w$$w.dat; \
			done; \
		done; \
	done

$(IDX_EVAL_DIR)/eval_map:
	mkdir -p $(IDX_EVAL_DIR)/map
	mkdir -p $(IDX_EVAL_DIR)/eval
	for t in 0 1000; do \
		for w in 4 7 15 30; do \
			rm -f $(IDX_EVAL_DIR)/map/map_t$$t.w$$w.out; \
			for k in 10 13 15 16 17; do \
				$(gedmap) align  $(GEDS) $(IDX_EVAL_DIR)/samples.1000000.fq  $(IDX_EVAL_DIR)/indexes/human.t$$t.$$w.$$k.min -a $(ADJ) -o  $(IDX_EVAL_DIR)/map/t$$t.w$$w.k$$k.sam   $(GEDMAP_PARAMS_SAM0) | tee -a $(IDX_EVAL_DIR)/map/map_t$$t.w$$w.out;\
			done; \
			grep mapping $(IDX_EVAL_DIR)/map/map_t$$t.w$$w.out | awk 'BEGIN{OFS="\t"; print "k", "mapping_time";a[1]=10;a[2]=13;a[3]=15;a[4]=16;a[5]=17}{print  a[NR], $$4 }' > $(IDX_EVAL_DIR)/eval/eval_map.t$$t.w$$w.dat; \
		done; \
	done

$(IDX_EVAL_DIR)/eval_qual:
	mkdir -p $(IDX_EVAL_DIR)/eval
	for t in 0 1000; do \
		for w in 4 7 15 30; do \
			rm -f $(IDX_EVAL_DIR)/eval/eval_qual.t$$t.w$$w.dat ; \
			../eval_sam /dev/null | head -n  1  | tee -a $(IDX_EVAL_DIR)/eval/eval_qual.t$$t.w$$w.dat; \
			for k in 10 13 15 16 17; do \
				../eval_sam  $(IDX_EVAL_DIR)/map/t$$t.w$$w.k$$k.sam   | tail -n 1  | tee -a $(IDX_EVAL_DIR)/eval/eval_qual.t$$t.w$$w.dat; \
			done; \
			awk 'BEGIN{FS=OFS="\t";a[1]="k";a[2]=10;a[3]=13;a[4]=15;a[5]=16;a[6]=17}{print  $$0,  a[NR] }' $(IDX_EVAL_DIR)/eval/eval_qual.t$$t.w$$w.dat > $(IDX_EVAL_DIR)/eval/eval_qual.t$$t.w$$w.dat.tmp; \
			mv $(IDX_EVAL_DIR)/eval/eval_qual.t$$t.w$$w.dat.tmp $(IDX_EVAL_DIR)/eval/eval_qual.t$$t.w$$w.dat; \
		done; \
	done

##################  GENERATE INDEX AND MAP AND EVAL

$(GEDMAP_DIR)/index.out: $(gedmap)
	mkdir -p $(GEDMAP_DIR)
	############ GEDMAP INDEX ##################
	$(time) -v -o $(GEDMAP_DIR)/index.time.out \
	make $(MINI) | tee $(GEDMAP_DIR)/index.out

$(GEDMAP_DIR)/map.out: $(gedmap)
	mkdir -p $(GEDMAP_DIR)
	############ GEDMAP MAP # ##################
	rm -f $(GEDMAP_DIR)/map.time.out
	for target in $(SAMPLES); \
	do \
	$(time) -v -o $(GEDMAP_DIR)/map.time.out -a \
	make $(GEDMAP_DIR)/$${target}.sam | tee -a $(GEDMAP_DIR)/map.out ; \
	done
	$(time) -v -o $(GEDMAP_DIR)/map.time.out -a \
	make $(GEDMAP_DIR)/reads.sam | tee -a $(GEDMAP_DIR)/map.out
	$(time) -v -o $(GEDMAP_DIR)/map.time.out -a \
	make $(GEDMAP_DIR)/reads2.sam | tee -a $(GEDMAP_DIR)/map.out

$(GEDMAP_DIR)/eval.out: $(GEDMAP_DIR)/map.out
	mkdir -p $(GEDMAP_DIR)
	../eval_sam /dev/null | head -n  1 > $(GEDMAP_DIR)/eval.out
	for target in $(SAMPLES); do ../eval_sam $(GEDMAP_DIR)/$${target}.sam | tail -n 1 >> $(GEDMAP_DIR)/eval.out; done
	../eval_sam $(GEDMAP_DIR)/reads.sam | tail -n 1 >> $(GEDMAP_DIR)/eval.out
	../eval_sam $(GEDMAP_DIR)/reads2.sam | tail -n 1 >> $(GEDMAP_DIR)/eval.out

$(GEDMAP_DIR)/clean:
	rm -f -r $(GEDMAP_DIR)

$(GEDMAP_DIR)/all: $(GEDMAP_DIR)/index.out $(GEDMAP_DIR)/map.out $(GEDMAP_DIR)/eval.out

gedmap_all: $(GEDMAP_DIR)/index.out samples/all $(GEDMAP_DIR)/map.out $(GEDMAP_DIR)/eval.out

##################### HISAT #####################


ht2-extr = $(HISAT_DIR)/human.extr
ht2-idx  = $(HISAT_DIR)/human.idx

$(HISAT_DIR)/hs_extract_build: $(VCF) $(FA)
	mkdir -p $(HISAT_DIR)
	$(hs_dir)/hisat2_extract_snps_haplotypes_VCF.py $(FA) $(VCF) $(ht2-extr)
	$(hs_dir)/hisat2-build -f --snp $(ht2-extr).snp --haplotype $(ht2-extr).haplotype $(FA) $(ht2-idx)

$(HISAT_DIR)/$(SAMPLE0).sam: $(SAMPLE_DIR)/$(SAMPLE0).fq $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(SAMPLE_DIR)/$(SAMPLE0).fq -S $(HISAT_DIR)/$(SAMPLE0).sam  $(HISAT_PARAMS_SAM0)
$(HISAT_DIR)/$(SAMPLE1).sam: $(SAMPLE_DIR)/$(SAMPLE1).fq $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(SAMPLE_DIR)/$(SAMPLE1).fq -S $(HISAT_DIR)/$(SAMPLE1).sam  $(HISAT_PARAMS_SAM1)
$(HISAT_DIR)/$(SAMPLE3).sam: $(SAMPLE_DIR)/$(SAMPLE3).fq $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(SAMPLE_DIR)/$(SAMPLE3).fq -S $(HISAT_DIR)/$(SAMPLE3).sam  $(HISAT_PARAMS_SAM3)
$(HISAT_DIR)/$(SAMPLE8).sam: $(SAMPLE_DIR)/$(SAMPLE8).fq $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(SAMPLE_DIR)/$(SAMPLE8).fq -S $(HISAT_DIR)/$(SAMPLE8).sam  $(HISAT_PARAMS_SAM8)
$(HISAT_DIR)/$(SAMPLE9).sam: $(SAMPLE_DIR)/$(SAMPLE9).fq $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(SAMPLE_DIR)/$(SAMPLE9).fq -S $(HISAT_DIR)/$(SAMPLE9).sam  $(HISAT_PARAMS_SAM9)
$(HISAT_DIR)/$(SAMPLE11).sam: $(SAMPLE_DIR)/$(SAMPLE11).fq $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(SAMPLE_DIR)/$(SAMPLE11).fq -S $(HISAT_DIR)/$(SAMPLE11).sam  $(HISAT_PARAMS_SAM11)
$(HISAT_DIR)/$(SAMPLE_OE).sam: $(SAMPLE_DIR)/$(SAMPLE_OE).fq $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(SAMPLE_DIR)/$(SAMPLE_OE).fq -S $(HISAT_DIR)/$(SAMPLE_OE).sam  $(HISAT_PARAMS_SAM1)
$(HISAT_DIR)/reads.sam: $(READS) $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(READS) -S $(HISAT_DIR)/reads.sam  $(HISAT_PARAMS_SAM1)
$(HISAT_DIR)/reads2.sam: $(READS) $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(READS2) -S $(HISAT_DIR)/reads2.sam  $(HISAT_PARAMS_SAM1)
#### in order to compare alignment-wise #####
$(HISAT_DIR)/reads.in-order.sam: $(READS) $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(READS) -S $(HISAT_DIR)/reads.in-order.sam  $(HISAT_PARAMS_SAM1) --reorder
$(HISAT_DIR)/reads2.in-order.sam: $(READS) $(ht2-idx).1.ht2
	$(hs_dir)/hisat2 -q -x $(ht2-idx) -U $(READS2) -S $(HISAT_DIR)/reads2.in-order.sam  $(HISAT_PARAMS_SAM1) --reorder

$(HISAT_DIR)/index.out:
	mkdir -p $(HISAT_DIR)
	########### HISAT INDEX ##################
	$(time) -v -o $(HISAT_DIR)/index.time.out \
	make $(HISAT_DIR)/hs_extract_build | tee $(HISAT_DIR)/index.out

$(HISAT_DIR)/map.out: $(gedmap)
	mkdir -p $(HISAT_DIR)
	############ HISAT MAP ###################
	rm -f $(HISAT_DIR)/map.time.out
	for target in $(SAMPLES); \
	do \
	$(time) -v -o $(HISAT_DIR)/map.time.out -a \
	make $(HISAT_DIR)/$${target}.sam | tee -a $(HISAT_DIR)/map.out ;\
	done
	$(time) -v -o $(HISAT_DIR)/map.time.out -a \
	make $(HISAT_DIR)/reads.sam | tee -a $(HISAT_DIR)/map.out
	$(time) -v -o $(HISAT_DIR)/map.time.out -a \
	make $(HISAT_DIR)/reads2.sam | tee -a $(HISAT_DIR)/map.out

$(HISAT_DIR)/eval.out: $(HISAT_DIR)/map.out
	mkdir -p $(HISAT_DIR)
	../eval_sam /dev/null | head -n  1 > $(HISAT_DIR)/eval.out
	for target in $(SAMPLES); do ../eval_sam $(HISAT_DIR)/$${target}.sam | tail -n 1 >> $(HISAT_DIR)/eval.out; done
	../eval_sam $(HISAT_DIR)/reads.sam | tail -n 1 >> $(HISAT_DIR)/eval.out
	../eval_sam $(HISAT_DIR)/reads2.sam | tail -n 1 >> $(HISAT_DIR)/eval.out

$(HISAT_DIR)/clean:
	rm -f -r $(HISAT_DIR)

$(HISAT_DIR)/all: $(HISAT_DIR)/index.out $(HISAT_DIR)/map.out $(HISAT_DIR)/eval.out

##################### VG GIRFAFFE #################################

vg = $(vg_dir)/bin/vg

vg_prefix    = $(GIRAFFE_DIR)/human
vg_graph  = $(vg_prefix).vg
xg_graph  = $(vg_prefix).xg
gfa_graph = $(vg_prefix).gfa
gbz       = $(vg_prefix).giraffe.gbz

$(GZ_DIR)/clean_all_vcf_tabix:
	rm -f $(GZ_DIR)/all_vcf_tabix
	rm -f $(GZ_DIR)/*.tbi

$(GZ_DIR)/all_vcf_tabix: $(GZ_DIR)/all_vcf_downloaded
	#DOWNLOAD CHROM 1 - 22 + X + Y
	for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do \
	tabix -p vcf $(GZ_DIR)/chrom$$i.vcf.gz &  \
	done; \
	tabix -p vcf $(GZ_DIR)/chromX.vcf.gz & \
	tabix -p vcf $(GZ_DIR)/chromY.vcf.gz & \
	wait
	touch $(GZ_DIR)/all_vcf_tabix

$(vg_graph): $(FA) $(GZ_DIR)/all_vcf_tabix
	$(vg) construct -r  $(FA) -v $(VCF_GZ) > $(vg_graph)

$(xg_graph): $(vg_graph)
	$(vg) index -x $(xg_graph) $(vg_graph)

$(gfa_graph): $(vg_graph)
	$(vg) view $(vg_graph) > $(gfa_graph)

$(gbz): $(gfa_graph) $(xg_graph)
	$(vg) autoindex --workflow giraffe -g $(gfa_graph) -p $(vg_prefix)

$(GIRAFFE_DIR)/$(SAMPLE0).sam:  $(SAMPLE_DIR)/$(SAMPLE0).fq $(gbz)
	$(vg) giraffe -Z $(gbz) -f $(SAMPLE_DIR)/$(SAMPLE0).fq $(GIRAFFE_PARAMS) > $(GIRAFFE_DIR)/$(SAMPLE0).sam
$(GIRAFFE_DIR)/$(SAMPLE1).sam:  $(SAMPLE_DIR)/$(SAMPLE1).fq $(gbz)
	$(vg) giraffe -Z $(gbz) -f $(SAMPLE_DIR)/$(SAMPLE1).fq $(GIRAFFE_PARAMS) > $(GIRAFFE_DIR)/$(SAMPLE1).sam
$(GIRAFFE_DIR)/$(SAMPLE3).sam:  $(SAMPLE_DIR)/$(SAMPLE3).fq $(gbz)
	$(vg) giraffe -Z $(gbz) -f $(SAMPLE_DIR)/$(SAMPLE3).fq $(GIRAFFE_PARAMS) > $(GIRAFFE_DIR)/$(SAMPLE3).sam
$(GIRAFFE_DIR)/$(SAMPLE8).sam:  $(SAMPLE_DIR)/$(SAMPLE8).fq $(gbz)
	$(vg) giraffe -Z $(gbz) -f $(SAMPLE_DIR)/$(SAMPLE8).fq $(GIRAFFE_PARAMS) > $(GIRAFFE_DIR)/$(SAMPLE8).sam
$(GIRAFFE_DIR)/$(SAMPLE9).sam:  $(SAMPLE_DIR)/$(SAMPLE9).fq $(gbz)
	$(vg) giraffe -Z $(gbz) -f $(SAMPLE_DIR)/$(SAMPLE9).fq $(GIRAFFE_PARAMS) > $(GIRAFFE_DIR)/$(SAMPLE9).sam
$(GIRAFFE_DIR)/$(SAMPLE11).sam:  $(SAMPLE_DIR)/$(SAMPLE11).fq $(gbz)
	$(vg) giraffe -Z $(gbz) -f $(SAMPLE_DIR)/$(SAMPLE11).fq $(GIRAFFE_PARAMS) >$(GIRAFFE_DIR)/$(SAMPLE11).sam
$(GIRAFFE_DIR)/$(SAMPLE_OE).sam:  $(SAMPLE_DIR)/$(SAMPLE_OE).fq $(gbz)
	$(vg) giraffe -Z $(gbz) -f $(SAMPLE_DIR)/$(SAMPLE_OE).fq $(GIRAFFE_PARAMS) > $(GIRAFFE_DIR)/$(SAMPLE_OE).sam
$(GIRAFFE_DIR)/reads.sam:  $(READS) $(gbz)
	$(vg) giraffe -Z $(gbz) -f $(READS) $(GIRAFFE_PARAMS) > $(GIRAFFE_DIR)/reads.sam
$(GIRAFFE_DIR)/reads2.sam:  $(READS) $(gbz)
	$(vg) giraffe -Z $(gbz) -f $(READS2) $(GIRAFFE_PARAMS) > $(GIRAFFE_DIR)/reads2.sam

$(GIRAFFE_DIR)/index.out:
	mkdir -p $(GIRAFFE_DIR)
	########### GIRAFFE_DIR INDEX ##################
	$(time) -v -o $(GIRAFFE_DIR)/index.time.out \
	make $(gbz) | tee $(GIRAFFE_DIR)/index.out

$(GIRAFFE_DIR)/map.out:
	mkdir -p $(GIRAFFE_DIR)
	############ GIRFAFFE MAP ###################
	rm -f $(GIRAFFE_DIR)/map.time.out
	for target in $(SAMPLES); \
	do \
	$(time) -v -o $(GIRAFFE_DIR)/map.time.out -a \
	make $(GIRAFFE_DIR)/$${target}.sam | tee -a $(GIRAFFE_DIR)/map.out ;\
	done
	$(time) -v -o $(GIRAFFE_DIR)/map.time.out -a \
	make $(GIRAFFE_DIR)/reads.sam | tee -a $(GIRAFFE_DIR)/map.out
	$(time) -v -o $(GIRAFFE_DIR)/map.time.out -a \
	make $(GIRAFFE_DIR)/reads2.sam | tee -a $(GIRAFFE_DIR)/map.out

$(GIRAFFE_DIR)/eval.out: $(GIRAFFE_DIR)/map.out
	mkdir -p $(GIRAFFE_DIR)
	../eval_sam /dev/null | head -n  1 > $(GIRAFFE_DIR)/eval.out
	for target in $(SAMPLES); do ../eval_sam $(GIRAFFE_DIR)/$${target}.sam | tail -n 1 >> $(GIRAFFE_DIR)/eval.out; done
	../eval_sam $(GIRAFFE_DIR)/reads.sam | tail -n 1 >> $(GIRAFFE_DIR)/eval.out
	../eval_sam $(GIRAFFE_DIR)/reads2.sam | tail -n 1 >> $(GIRAFFE_DIR)/eval.out

$(GIRAFFE_DIR)/all: $(GIRAFFE_DIR)/index.out $(GIRAFFE_DIR)/map.out $(GIRAFFE_DIR)/eval.out

$(GIRAFFE_DIR)/clean:	$(GZ_DIR)/clean_all_vcf_tabix
	rm -f -r $(GIRAFFE_DIR)
	rm -f $(FA).fai
	
##################### 
../comp_sams:
	cd .. && make comp_sams
	
compare_sam_ged_hisat2: $(HISAT_DIR)/reads.in-order.sam $(HISAT_DIR)/reads2.in-order.sam $(GEDMAP_DIR)/reads.in-order.sam $(GEDMAP_DIR)/reads2.in-order.sam ../comp_sams
	echo gedmap vs hisat read set 1 > compare_sam_ged_hisat2
	../comp_sams ../data/reads.fq gedmap/reads.in-order.sam hisat2/reads.in-order.sam >> compare_sam_ged_hisat2
	echo gedmap vs hisat read set 2 >> compare_sam_ged_hisat2
	../comp_sams ../data/reads2.fq gedmap/reads2.in-order.sam hisat2/reads2.in-order.sam >> compare_sam_ged_hisat2
	cat compare_sam_ged_hisat2
	
#####################
