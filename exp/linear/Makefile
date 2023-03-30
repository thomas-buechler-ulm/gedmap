include ../../PATHS
#################################################
################# EXPERIMENT ####################
#################################################
# use following commands to run experiments of  #
# the paper (please ensure you downloaded the   #
# data before)                                  #
#                                               #
# make                                          #
#                                               #
# to only run the experiments for gedmap use    #
#                                               #
# make gedmap/all                               #
#                                               #
#################################################
# the evaluation of mapping rate and accuracy   #
# is stored in subdir/map.eval                  #
# the evaluation of time and memory is stored   #
# in subdir/map.out                             #
#################################################

gedmap   = ../..

#OUTPUT DIRECTORIES
GEDMAP_DIR  = gedmap
HISAT_DIR   = hisat2
GIRAFFE_DIR = giraffe
MINIMAP_DIR = minimap2
samples     = ../samples

#INPUT DATA
DATA = $(gedmap)/data
FA   = $(DATA)/hs37d5.fa
VCF  = /dev/null
GZ_DIR =$(DATA)/gz


all: $(GEDMAP_DIR)/all $(GIRAFFE_DIR)/all $(HISAT_DIR)/all $(MINIMAP_DIR)/all
$(GEDMAP_DIR)/all:   $(GEDMAP_DIR)/index.out  $(GEDMAP_DIR)/map.eval   $(GEDMAP_DIR)/map.out 
$(HISAT_DIR)/all:    $(HISAT_DIR)/index.out   $(HISAT_DIR)/map.eval    $(HISAT_DIR)/map.out 
$(GIRAFFE_DIR)/all:  $(GIRAFFE_DIR)/index.out $(GIRAFFE_DIR)/map.eval  $(GIRAFFE_DIR)/map.out 
$(MINIMAP_DIR)/all:  $(MINIMAP_DIR)/index.out $(MINIMAP_DIR)/map.eval  $(MINIMAP_DIR)/map.out 


############################ SAMPLES ############################

../samples: 
	cd .. && make samples

############################  GEDMAP   ############################
#BUILD

geds = $(GEDMAP_DIR)/human.geds

$(GEDMAP_DIR)/index.out: $(FA) $(VCF)
	mkdir -p $(GEDMAP_DIR)
	{ /usr/bin/time -v  $(gedmap)/gedmap parse $(FA) $(VCF) -o $(geds)                       ; } 2>&1 | tee    $(GEDMAP_DIR)/index.out
	{ /usr/bin/time -v $(gedmap)/gedmap index  $(geds) -2fa $(geds).2fa -t 64                ; } 2>&1 | tee -a $(GEDMAP_DIR)/index.out

#MAP
$(GEDMAP_DIR)/map.out: ../samples
	{ /usr/bin/time -v $(gedmap)/gedmap align $(samples)/set_100.fq  $(geds) $(geds).min -2fa $(geds).2fa -o /dev/null -rc -t 8                 ; } 2>&1 | tee $(GEDMAP_DIR)/maps1.out
	{ /usr/bin/time -v $(gedmap)/gedmap align $(samples)/set_250.fq  $(geds) $(geds).min -2fa $(geds).2fa -o /dev/null -rc -t 8                 ; } 2>&1 | tee $(GEDMAP_DIR)/maps2.out
	{ /usr/bin/time -v $(gedmap)/gedmap align $(samples)/set_cnv.fq  $(geds) $(geds).min -2fa $(geds).2fa -o /dev/null -rc -t 8                 ; } 2>&1 | tee $(GEDMAP_DIR)/mapcn.out
	{ /usr/bin/time -v $(gedmap)/gedmap align $(DATA)/reads1.fq      $(geds) $(geds).min -2fa $(geds).2fa -o /dev/null -rc -t 8 -d 100          ; } 2>&1 | tee $(GEDMAP_DIR)/mapr1.out
	{ /usr/bin/time -v $(gedmap)/gedmap align $(DATA)/reads2.fq      $(geds) $(geds).min -2fa $(geds).2fa -o /dev/null -rc -t 8 -d 100          ; } 2>&1 | tee $(GEDMAP_DIR)/mapr2.out

	$(call get_time_and_space,  $(GEDMAP_DIR))


#EVALUATE
$(GEDMAP_DIR)/map.eval: ../samples $(DATA)/reads1a.fq $(DATA)/reads2a.fq
	$(gedmap)/gedmap align $(samples)/set_100a.fq  $(geds)  $(geds).min -2fa $(geds).2fa -o $(GEDMAP_DIR)/tmps1.sam -rc -t 8
	$(gedmap)/gedmap align $(samples)/set_250a.fq  $(geds)  $(geds).min -2fa $(geds).2fa -o $(GEDMAP_DIR)/tmps2.sam -rc -t 8
	$(gedmap)/gedmap align $(samples)/set_cnva.fq  $(geds)  $(geds).min -2fa $(geds).2fa -o $(GEDMAP_DIR)/tmpcn.sam -rc -t 8
	$(gedmap)/gedmap align $(DATA)/reads1a.fq      $(geds)  $(geds).min -2fa $(geds).2fa -o $(GEDMAP_DIR)/tmpr1.sam -rc -t 8 -d 100
	$(gedmap)/gedmap align $(DATA)/reads2a.fq      $(geds)  $(geds).min -2fa $(geds).2fa -o $(GEDMAP_DIR)/tmpr2.sam -rc -t 8 -d 100

	$(call all_my_eval, GEDMAP, $(GEDMAP_DIR))


############################ HISAT2 ############################
#BUILD

hs_idx=$(HISAT_DIR)/human.idx

$(HISAT_DIR)/index.out:  $(FA)
	mkdir -p $(HISAT_DIR)
	{ /usr/bin/time -v $(hisat2_dir)/hisat2-build -f --large-index --bmax 10000000 --dcv 4096 --threads 64 $(FA) $(HISAT_DIR)/human.idx ; } 2>&1 | tee -a $(HISAT_DIR)/index.out


#MAP
$(HISAT_DIR)/map.out: ../samples
	{ /usr/bin/time -v $(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(samples)/set_100.fq -S /dev/null --threads 8 --sensitive --np 0 -k 1 ; } 2>&1 | tee $(HISAT_DIR)/maps1.out
	{ /usr/bin/time -v $(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(samples)/set_250.fq -S /dev/null --threads 8 --sensitive --np 0 -k 1 ; } 2>&1 | tee $(HISAT_DIR)/maps2.out
	{ /usr/bin/time -v $(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(samples)/set_cnv.fq -S /dev/null --threads 8 --sensitive --np 0 -k 1 ; } 2>&1 | tee $(HISAT_DIR)/mapcn.out
	{ /usr/bin/time -v $(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(DATA)/reads1.fq     -S /dev/null --threads 8 --sensitive --np 0 -k 1 ; } 2>&1 | tee $(HISAT_DIR)/mapr1.out
	{ /usr/bin/time -v $(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(DATA)/reads2.fq     -S /dev/null --threads 8 --sensitive --np 0 -k 1 ; } 2>&1 | tee $(HISAT_DIR)/mapr2.out

	$(call get_time_and_space,  $(HISAT_DIR))

$(HISAT_DIR)/map.eval: ../samples
	$(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(samples)/set_100a.fq -S $(HISAT_DIR)/tmps1.sam --threads 8  --sensitive  --np 0 -k 1 --reorder
	$(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(samples)/set_250a.fq -S $(HISAT_DIR)/tmps2.sam --threads 8  --sensitive  --np 0 -k 1 --reorder
	$(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(samples)/set_cnva.fq -S $(HISAT_DIR)/tmpcn.sam --threads 8  --sensitive  --np 0 -k 1 --reorder
	$(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(DATA)/reads1a.fq     -S $(HISAT_DIR)/tmpr1.sam --threads 8  --sensitive  --np 0 -k 1 --reorder
	$(hisat2_dir)/hisat2 -q -x $(hs_idx) -U $(DATA)/reads2a.fq     -S $(HISAT_DIR)/tmpr2.sam --threads 8  --sensitive  --np 0 -k 1 --reorder

	$(call all_my_eval, HISAT2, $(HISAT_DIR))



############################ VG ############################
#BUILD


vg_prefix = $(GIRAFFE_DIR)/human
vg_graph  = $(vg_prefix).vg
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


$(GIRAFFE_DIR)/index.out:
	mkdir -p $(GIRAFFE_DIR)
	{ /usr/bin/time -v $(vg) construct -r  $(FA) > $(vg_graph) -t 64 ; } 2>&1 | tee $(GIRAFFE_DIR)/index.out
	{ /usr/bin/time -v $(vg) index -x $(vg_prefix).xg $(vg_graph) -L -t 64 ; } 2>&1 | tee -a $(GIRAFFE_DIR)/index.out
	{ /usr/bin/time -v $(vg) view $(vg_graph) --threads 64 > $(vg_prefix).gfa  ; } 2>&1 | tee -a $(GIRAFFE_DIR)/index.out
	{ /usr/bin/time -v $(vg) autoindex --workflow giraffe -g $(vg_prefix).gfa -p $(vg_prefix)  -t 64 ; } 2>&1 | tee -a $(GIRAFFE_DIR)/index.out


$(GIRAFFE_DIR)/map.out: ../samples
	{ /usr/bin/time -v  $(vg) giraffe -Z $(gbz) -f $(samples)/set_100.fq -t 8 -o SAM > /dev/null ; } 2>&1 | tee $(GIRAFFE_DIR)/maps1.out
	{ /usr/bin/time -v  $(vg) giraffe -Z $(gbz) -f $(samples)/set_250.fq -t 8 -o SAM > /dev/null ; } 2>&1 | tee $(GIRAFFE_DIR)/maps2.out
	{ /usr/bin/time -v  $(vg) giraffe -Z $(gbz) -f $(samples)/set_cnv.fq -t 8 -o SAM > /dev/null ; } 2>&1 | tee $(GIRAFFE_DIR)/mapcn.out
	{ /usr/bin/time -v  $(vg) giraffe -Z $(gbz) -f $(DATA)/reads1.fq     -t 8 -o SAM > /dev/null ; } 2>&1 | tee $(GIRAFFE_DIR)/mapr1.out
	{ /usr/bin/time -v  $(vg) giraffe -Z $(gbz) -f $(DATA)/reads2.fq     -t 8 -o SAM > /dev/null ; } 2>&1 | tee $(GIRAFFE_DIR)/mapr2.out

	$(call get_time_and_space,  $(GIRAFFE_DIR))

$(GIRAFFE_DIR)/map.eval: ../samples
	$(vg) giraffe -Z $(gbz) -f $(samples)/set_100a.fq -t 8 -o SAM > $(GIRAFFE_DIR)/tmps1.sam
	$(vg) giraffe -Z $(gbz) -f $(samples)/set_250a.fq -t 8 -o SAM > $(GIRAFFE_DIR)/tmps2.sam
	$(vg) giraffe -Z $(gbz) -f $(samples)/set_cnva.fq -t 8 -o SAM > $(GIRAFFE_DIR)/tmpcn.sam
	$(vg) giraffe -Z $(gbz) -f $(DATA)/reads1a.fq     -t 8 -o SAM > $(GIRAFFE_DIR)/tmpr1.sam
	$(vg) giraffe -Z $(gbz) -f $(DATA)/reads2a.fq     -t 8 -o SAM > $(GIRAFFE_DIR)/tmpr2.sam

	$(call all_my_eval, GIRAFFE, $(GIRAFFE_DIR))

######################## Minimap ##########################
#BUILD


$(MINIMAP_DIR)/index.out:
	mkdir -p $(MINIMAP_DIR)
	{ /usr/bin/time -v $(minimap2) -x sr -d $(MINIMAP_DIR)/human.min $(FA) ; } 2>&1 | tee $(MINIMAP_DIR)/index.out

$(MINIMAP_DIR)/map.out: ../samples
	{ /usr/bin/time -v  $(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(samples)/set_100.fq > /dev/null ; } 2>&1 | tee minimap2/maps1.out
	{ /usr/bin/time -v  $(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(samples)/set_250.fq > /dev/null ; } 2>&1 | tee minimap2/maps2.out
	{ /usr/bin/time -v  $(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(samples)/set_cnv.fq > /dev/null ; } 2>&1 | tee minimap2/mapcn.out
	{ /usr/bin/time -v  $(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(DATA)/reads1.fq     > /dev/null ; } 2>&1 | tee minimap2/mapr1.out
	{ /usr/bin/time -v  $(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(DATA)/reads2.fq     > /dev/null ; } 2>&1 | tee minimap2/mapr2.out

	$(call get_time_and_space,  $(MINIMAP_DIR))

$(MINIMAP_DIR)/map.eval: ../samples
	$(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(samples)/set_100a.fq > $(MINIMAP_DIR)/tmps1.sam
	$(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(samples)/set_250a.fq > $(MINIMAP_DIR)/tmps2.sam
	$(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(samples)/set_cnva.fq > $(MINIMAP_DIR)/tmpcn.sam
	$(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(DATA)/reads1a.fq     > $(MINIMAP_DIR)/tmpr1.sam
	$(minimap2) -a $(MINIMAP_DIR)/human.min -t 8 $(DATA)/reads2a.fq     > $(MINIMAP_DIR)/tmpr2.sam

	$(call all_my_eval, MINIMAP2, $(MINIMAP_DIR))

###################################################################


VCF_GZ_V  = $(GZ_DIR)/chrom1.vcf.gz \
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

VCF_GZ_H  = $(GZ_DIR)/chrom1.vcf.gz,\
$(GZ_DIR)/chrom2.vcf.gz,\
$(GZ_DIR)/chrom3.vcf.gz,\
$(GZ_DIR)/chrom4.vcf.gz,\
$(GZ_DIR)/chrom5.vcf.gz,\
$(GZ_DIR)/chrom6.vcf.gz,\
$(GZ_DIR)/chrom7.vcf.gz,\
$(GZ_DIR)/chrom8.vcf.gz,\
$(GZ_DIR)/chrom10.vcf.gz,\
$(GZ_DIR)/chrom11.vcf.gz,\
$(GZ_DIR)/chrom12.vcf.gz,\
$(GZ_DIR)/chrom13.vcf.gz,\
$(GZ_DIR)/chrom14.vcf.gz,\
$(GZ_DIR)/chrom15.vcf.gz,\
$(GZ_DIR)/chrom16.vcf.gz,\
$(GZ_DIR)/chrom17.vcf.gz,\
$(GZ_DIR)/chrom18.vcf.gz,\
$(GZ_DIR)/chrom19.vcf.gz,\
$(GZ_DIR)/chrom20.vcf.gz,\
$(GZ_DIR)/chrom21.vcf.gz,\
$(GZ_DIR)/chrom22.vcf.gz,\
$(GZ_DIR)/chromX.vcf.gz,\
$(GZ_DIR)/chromY.vcf.gz


############################ FUNCTIONS ############################

define get_time
	grep Elapsed $(1) | awk '{print  "$(2)" , $$8}' >>  $(3)
endef

define get_space
	grep Maximum $(1) | awk '{print $$6}' >>  $(2)
endef

define get_time_and_space
	rm -f $(1)/map1.out $(1)/map2.out
	$(call get_time,  $(1)/maps1.out, s1, $(1)/map1.out)
	$(call get_time,  $(1)/maps2.out, s2, $(1)/map1.out)
	$(call get_time,  $(1)/mapcn.out, cnv, $(1)/map1.out)
	$(call get_time,  $(1)/mapr1.out, r1, $(1)/map1.out)
	$(call get_time,  $(1)/mapr2.out, r2, $(1)/map1.out)
	$(call get_space, $(1)/maps1.out,     $(1)/map2.out)
	$(call get_space, $(1)/maps2.out,     $(1)/map2.out)
	$(call get_space, $(1)/mapcn.out,     $(1)/map2.out)
	$(call get_space, $(1)/mapr1.out,     $(1)/map2.out)
	$(call get_space, $(1)/mapr2.out,     $(1)/map2.out)
	paste -d ' ' $(1)/map1.out $(1)/map2.out > $(1)/map.out
	rm -f $(1)/map1.out $(1)/map2.out
endef

define my_eval
	$(gedmap)/eval_sam $(1) | tail -n 1 | awk '{ printf("%.3f\n%.3f\n" ,  $$3/1000000 , $$5/1000000)}' >> $(2)
endef

define my_transpose
	awk 'BEGIN { print "ALG mapped100 correct100 mapped250 correct250 mappedcnv correctcnv  mappedr1 correctr1 mappedr2 correctr2" } { printf("%s ", $$1) } END { printf("\n") }' $(1) > $(2)
endef

define all_my_eval
	echo $(1) > $(2)/map.eval.tmp
	$(call my_eval, $(2)/tmps1.sam , $(2)/map.eval.tmp)
	$(call my_eval, $(2)/tmps2.sam , $(2)/map.eval.tmp)
	$(call my_eval, $(2)/tmpcn.sam , $(2)/map.eval.tmp)
	$(call my_eval, $(2)/tmpr1.sam , $(2)/map.eval.tmp)
	$(call my_eval, $(2)/tmpr2.sam , $(2)/map.eval.tmp)
	$(call my_transpose, $(2)/map.eval.tmp , $(2)/map.eval)
	rm  $(2)/map.eval.tmp
	cat $(2)/map.eval
endef
