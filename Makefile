#################################################
################ COMPLILE PROGRAMS  #############
#################################################
# use following command to compile              #
#                                               #
# make                                          #
#                                               #
#################################################
include PATHS
include $(sdsl_dir)/Make.helper
LIBS = -lsdsl -ldivsufsort -ldivsufsort64 -lstdc++fs
SRC = src
LIB  = $(SRC)/lib
PROG = $(SRC)/programs
OPTIONS = $(MY_CXX) $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -no-pie $(C_OPTIONS) -I$(INC_DIR) -L$(LIB_DIR) -std=c++2a

#################################################
##############   COMPILE PROGRAM   ##############
#################################################

help:
	# 'make bin' to compile
	# 'make data_download' to download data
	# 'make all_experiments' to make all experiments (set paths to vg and hisat2 before)
	# 'make gedmap_eperiments' to only run experiments with gedmap
	# 'make a_little_example' to run the little_example
	# 'make clean' to remove everything

bin: gedmap eval_sam comp_sams

gedmap: $(SRC)/* $(LIB)/* $(PROG)/*
	$(OPTIONS) -fopenmp \
	$(SRC)/gedmap.cpp -o gedmap \
	$(LIBS)

eval_sam: $(SRC)/evaluate_sam.cpp
	$(MY_CXX) $(SRC)/evaluate_sam.cpp -o eval_sam

comp_sams: $(SRC)/compare_sams.cpp
	$(MY_CXX) $(SRC)/compare_sams.cpp -o comp_sams

eval_idx: $(SRC)/evaluate_my_index.cpp
	$(OPTIONS) -fopenmp  $(SRC)/evaluate_my_index.cpp -o eval_idx $(LIBS)

#################################################

data_download:
	cd data && make

all_eperiments: bin data_download
	cd exp && make

gedmap_eperiments: bin data_download
	cd exp && make gedmap_all

a_little_example: bin
	cd a_little_example && make

#################################################

clean:
	rm -f gedmap eval_sam
	cd data && make clean
	cd exp && make clean
	cd a_little_example && make clean
