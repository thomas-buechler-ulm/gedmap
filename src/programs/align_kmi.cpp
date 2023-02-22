namespace gedmap_align_kmi{
//// STATE AS SUBMITTED TO WABI22. AN KMER INDEX WAS USED INSTEAD OF MINIMIZERS.

vector<sys_timer> tv(10);

const uint32_t tv_ALL  = 0;
const uint32_t tv_LOAD = 1;
const uint32_t tv_MAP  = 2;

//PARALELL TIME
const uint32_t tv_SEED_P = 3;
const uint32_t tv_EOCQ_P = 4;
const uint32_t tv_HSF_P  = 5;
const uint32_t tv_ALI_P  = 6;
const uint32_t tv_WRT_P  = 7;
const uint32_t tv_ALL_P  = 8;
const uint32_t tv_CHA_P  = 9;


template<class int_type>
void paralell_processor( linear_eoc_type & eoc, const std::string & EDS, const adjacency & ADJ, const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s);

template<class int_type>
void paralell_processor2( linear_eoc_type & eoc, const std::string & EDS, const adjacency & ADJ, const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s);

template<class int_type>
std::vector<fasta_read<int_type>> read_batch(std::ifstream & fastq_s);

template<class int_type>
fasta_read<int_type> get_read_from_stream(std::ifstream & fastq_s);

int main(int argc,  char** argv){
	using namespace gedmap_io;
	using namespace std;
	using namespace sdsl;
	print_prog_headline("GEDMAP ALIGN");

	tv[tv_ALL].start();
	tv[tv_LOAD].start();

	linear_eoc_type eoc;
	string EDS;
	adjacency ADJ = adjacency();
	pos_EDS_to_FA_type p2FA;

	ifstream fastq_in_stream;
	ofstream out_stream;

	handle_input(argc, argv, eoc, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	
	print_row("Search", argv[3]);
	print_row("in", argv[2]);
	print_row("using", argv[4]);
	dotline();

	tv[tv_LOAD].stop();
	tv[tv_MAP].start();

	const uint32_t int_width = bitsneeded(EDS.size());

	// DIFFERENTIATE BETWEEN INTEGER SIZES
	if		(int_width <=  8 ) paralell_processor2< uint8_t>(eoc, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	else if	(int_width <= 16 ) paralell_processor2<uint16_t>(eoc, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	else if	(int_width <= 32 ) paralell_processor2<uint32_t>(eoc, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	else if	(int_width <= 64 ) paralell_processor2<uint64_t>(eoc, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
	else throw runtime_error (" in align main: EDS too long");

	tv[tv_MAP].stop();

	print_row("Time for loading:",		tv[tv_LOAD]		.get(), " s");
	print_row("Time for mapping:",		tv[tv_MAP]		.get(), " s");
	print_row("-Time for seed determination:",tv[tv_SEED_P]	.get(), " s");
// 	print_row("-Time for eoc_queries:",		tv[tv_EOCQ_P]	.get(), " s");
	print_row("-Time for hotspots:",		tv[tv_HSF_P]	.get(), " s");
	print_row("-Time for aligning:",		tv[tv_ALI_P]	.get(), " s");
	print_row("-Time for chances:",		tv[tv_CHA_P]	.get(), " s");
	print_row("-Time for I/O (fastq+sam):",	tv[tv_WRT_P]	.get(), " s");
	print_row("Time for overall:", 		tv[tv_ALL].stop_and_get(), " s");
	dotline();
	return 0;
}

/** @brief 
 * read read from fastq file
 * map read in the index
 * write alignment to sam file
 */
template<class int_type>
void paralell_processor2( linear_eoc_type & eoc, const std::string & EDS, const adjacency & ADJ, const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s){

	environment<
		linear_eoc_type,
		std::string, // EDS
		adjacency,
		pos_EDS_to_FA_type> env(eoc, EDS, ADJ, p2FA);

// 	if(MAP_RC) BATCH_SIZE *= 2;
	
	//count number of reads in file
	string s;
	uint32_t read_count = 0;
	while( getline( fastq_s, s ) ) read_count++;
	read_count /= 4;
	fastq_s.clear();
	fastq_s.seekg(0);
	if(MAP_RC) read_count *= 2;
	

	gedmap_io::print_row("MAPPING:");
	uint32_t results	= 0;
	uint32_t mapped 	= 0;
	uint32_t count  	= 0;
	uint32_t batchcount	= 0;
	uint32_t chances_count	= 0;

	if(THREAD_COUNT) omp_set_num_threads(THREAD_COUNT);
	
	bool run = true;
	while(run){
		batchcount++;
		
		tv[tv_CHA_P].start();
		if(!MAP_RC){
			#pragma omp parallel for  schedule(dynamic,10)
			for(uint32_t r = 0; r < read_count; r++){
				fasta_read<int_type> read;
				bool ok = false;
				
				#pragma omp critical
				{
					if(fastq_s.good()){
						try{
							read = fasta_read<int_type>(fastq_s);
							ok = true;
						}catch(runtime_error & e){
							//NO READ LEFT
							run = false;
						}
					}
				}
				
				if(!ok) continue;
				
				uint32_t my_searches_chances = 0;
				
				std::vector<alignment<int_type>> alignments;
				for (size_t i = 0;
					i < FRAGMENT_COUNT.size() && (alignments.empty() || alignments[0].dist >= DOUBT_DIST);
					i++)
				{
					my_searches_chances++;

					uint32_t max_dist = MAX_DIST;
					if (!alignments.empty())
						max_dist = std::min(max_dist, alignments[0].dist - 1);

					auto hotspots = read_processor::find_hotspots(read, eoc,
										FRAGMENT_COUNT[i],
										SPOT_SIZE, SPOT_HITS[i], CHECK_COLLI);
					read_processor::append_alignments(
						alignments,
						read_processor::start_aligner(std::move(hotspots), read,
										env,
										max_dist,
										MAX_ALIGNS_C, MAX_ALIGNS_T, MAX_ALIGNS_O));
				}
				#pragma omp critical
				{
					const auto res = read_processor::write_alignment(
						std::move(alignments),
						read,
						o_s,
						WRITE_FAILURE,
						p2FA);
					results += res;
					if(res) mapped++;
					count++;
					if(!(count%1000))
						gedmap_io::flush_row("Searched reads", to_string(count));
					chances_count += my_searches_chances;
				}
			}
		}else{
			#pragma omp parallel for  schedule(dynamic,10)
			for(uint32_t r = 0; r < read_count; r+=2){
				
				fasta_read<int_type> read;
				bool ok = false;
				
				#pragma omp critical
				{
					if(fastq_s.good()){
						try{
							read = fasta_read<int_type>(fastq_s);
							ok = true;
						}catch(runtime_error & e){
							//NO READ LEFT
							run = false;
						}
					}
				}
				
				if(!ok) continue;
				uint32_t my_searches_chances = 0;
				
				vector<alignment<int_type>> alignments;
				for (size_t i = 0;
					 i < FRAGMENT_COUNT.size() && (alignments.empty() || alignments[0].dist >= DOUBT_DIST);
					 i++)
				{
					my_searches_chances++;

					uint32_t max_dist = MAX_DIST;
					if (!alignments.empty())
						max_dist = std::min(max_dist, alignments[0].dist - 1);
					
					const auto& read_rev = read.get_rev_compl();
					auto hotspots     = read_processor::find_hotspots(read, eoc,
											FRAGMENT_COUNT[i],
											SPOT_SIZE, SPOT_HITS[i], CHECK_COLLI);
					auto hotspots_r_c = read_processor::find_hotspots(read_rev, eoc,
											FRAGMENT_COUNT[i],
											SPOT_SIZE, SPOT_HITS[i], CHECK_COLLI);
					read_processor::append_alignments(
						alignments,
						read_processor::start_aligner(
							std::move(hotspots),     read,
							std::move(hotspots_r_c), read_rev,
							env, max_dist,
							MAX_ALIGNS_C, MAX_ALIGNS_T, MAX_ALIGNS_O));
				}
				#pragma omp critical
				{
					const auto res = read_processor::write_alignment(
						std::move(alignments),
						read,
						o_s,
						WRITE_FAILURE,
						p2FA);
					results += res;
					if(res) mapped++;
					count++;
					if(!(count%1000))
						gedmap_io::flush_row("Searched reads", to_string(count));
					chances_count += my_searches_chances;
				}
			}
		}
		tv[tv_CHA_P].stop();
	}
	fastq_s.close();
	o_s.close();

	gedmap_io::print_row("Searched reads", count);
	gedmap_io::print_row("Searches", chances_count);
	gedmap_io::print_row("Mapped reads", mapped);
	gedmap_io::print_row("Results", results);
	gedmap_io::dotline();
}

}

