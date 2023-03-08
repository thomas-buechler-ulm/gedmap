namespace gedmap_align_min {

vector<sys_timer> tv(3);

const uint32_t tv_ALL  = 0;
const uint32_t tv_LOAD = 1;
const uint32_t tv_MAP  = 2;


template<class int_type>
void paralell_processor( gedmap_mini::minimizer_index & eoc, const std::string & EDS, const adjacency & ADJ, const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s);

template<class int_type>
void paralell_processor_batch( gedmap_mini::minimizer_index & eoc, const std::string & EDS, const adjacency & ADJ, const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s);

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

	gedmap_mini::minimizer_index mini;
	string EDS;
	adjacency ADJ = adjacency();
	pos_EDS_to_FA_type p2FA;

	ifstream fastq_in_stream, fastq_in_stream2;
	ofstream out_stream;

	handle_input(argc, argv, mini, EDS, ADJ, p2FA, fastq_in_stream, fastq_in_stream2, out_stream);
	assert(ADJ.initialised);

	
	print_row("Search", argv[3]);
	print_row("in", argv[2]);
	print_row("using", argv[4]);
	if (fastq_in_stream2.is_open())
		print_row("Running in paired-end mode");
	dotline();

	ED_Graph<uint32_t> graph;
	if (fastq_in_stream2.is_open()) {
		// because sdsl::rank_support doesnt have move-assignment
		graph.~ED_Graph<uint32_t>();
		new(&graph) ED_Graph<uint32_t>(EDS);
	}

	tv[tv_LOAD].stop();
	tv[tv_MAP].start();

	const uint32_t int_width = std::max( bitsneeded(EDS.size()), static_cast<uint32_t>(2*mini.k) );
	edsm_align::init_lookuptables();

	// DIFFERENTIATE BETWEEN INTEGER SIZES
	if (fastq_in_stream2.is_open()) {
		// paired end mapping
		if		(int_width <=  8 ) map_pairs< uint8_t>(mini, graph, EDS, ADJ, p2FA, fastq_in_stream, fastq_in_stream2, out_stream);
		else if	(int_width <= 16 ) map_pairs<uint16_t>(mini, graph, EDS, ADJ, p2FA, fastq_in_stream, fastq_in_stream2, out_stream);
		else if	(int_width <= 32 ) map_pairs<uint32_t>(mini, graph, EDS, ADJ, p2FA, fastq_in_stream, fastq_in_stream2, out_stream);
		else if	(int_width <= 64 ) map_pairs<uint64_t>(mini, graph, EDS, ADJ, p2FA, fastq_in_stream, fastq_in_stream2, out_stream);
		else throw runtime_error (" in align main: EDS too long");
	} else {
		if		(int_width <=  8 ) paralell_processor< uint8_t>(mini, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
		else if	(int_width <= 16 ) paralell_processor<uint16_t>(mini, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
		else if	(int_width <= 32 ) paralell_processor<uint32_t>(mini, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
		else if	(int_width <= 64 ) paralell_processor<uint64_t>(mini, EDS, ADJ, p2FA, fastq_in_stream, out_stream);
		else throw runtime_error (" in align main: EDS too long");
	}

	tv[tv_MAP].stop();

	print_row("Time for loading:",		tv[tv_LOAD]		.get(), " s");
	print_row("Time for mapping:",		tv[tv_MAP]		.get(), " s");
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
void paralell_processor( gedmap_mini::minimizer_index & mini, const std::string & EDS, const adjacency & ADJ,const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s){

	if(IN_ORDER){
		paralell_processor_batch<int_type>(mini,EDS,ADJ,p2FA,fastq_s,o_s);
		return;
	}

	environment<
		gedmap_mini::minimizer_index,
		std::string, // EDS
		adjacency,
		pos_EDS_to_FA_type> env(mini, EDS, ADJ, p2FA);
	
	string s;
	size_t read_count = 0;
	while( getline( fastq_s, s ) ) read_count++;
	read_count /= 4;
	fastq_s.clear(); // reset error flags
	fastq_s.seekg(0, std::ios::beg);

	gedmap_io::print_row("MAPPING:");
	uint32_t results	= 0;
	uint32_t mapped 	= 0;
	uint32_t count  	= 0;
	uint32_t chances_count	= 0;
	

	if (THREAD_COUNT) omp_set_num_threads(THREAD_COUNT);
	//std::vector<std::ostringstream> outstreams(omp_get_max_threads());

	// set large buffers for o_s and fastq_s. Mostly useful for many threads
	std::vector<char> outbuffer(1<<23);
	o_s.rdbuf()->pubsetbuf(outbuffer.data(), outbuffer.size());
	//std::vector<char> inbuffer(1<<30);
	//fastq_s.rdbuf()->pubsetbuf(inbuffer.data(), inbuffer.size());

	
	// TODO: maybe use one thread only for reading/parsing data for larger THREAD_COUNT?
	if(!MAP_RC){
		#pragma omp parallel for  schedule(dynamic,10)
		for (size_t r = 0; r < read_count; r++){
			fasta_read<int_type> read;
			#pragma omp critical
			{
				read = fasta_read<int_type>(fastq_s);
			}
			
			uint32_t my_searches_chances = 0;
			
			vector<alignment<int_type>> alignments;
			for (size_t i = 0;
				i < FRAGMENT_COUNT.size() && (alignments.empty() || alignments[0].dist >= DOUBT_DIST);
				i++)
			{
				uint32_t max_dist = MAX_DIST[i];
				if (!alignments.empty())
					max_dist = std::min(max_dist, alignments[0].dist - 1);

				my_searches_chances++;
				auto hotspots = read_processor::find_hotspots(read, mini,
									FRAGMENT_COUNT[i],
									SPOT_SIZE, SPOT_HITS[i], CHECK_COLLI);
				read_processor::append_alignments(
					alignments,
					read_processor::start_aligner(std::move(hotspots), read,
									env,
									max_dist,
									MAX_ALIGNS_C[i], MAX_ALIGNS_T[i], MAX_ALIGNS_O));
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
				if (res > 0) mapped++;
				count++;
				if(!(count%1000))
					gedmap_io::flush_row("Searched reads", to_string(count));
				chances_count += my_searches_chances;
			}
		}
	}else{
		#pragma omp parallel for  schedule(dynamic,10)
		for (uint32_t r = 0; r < read_count; r++){
			
			fasta_read<int_type> read;
			
			#pragma omp critical
			{
				read = fasta_read<int_type>(fastq_s);
			}
			
			uint32_t my_searches_chances = 0;
			
			vector<alignment<int_type>> alignments;
			for (size_t i = 0;
				i < FRAGMENT_COUNT.size() && (alignments.empty() || alignments[0].dist >= DOUBT_DIST);
				i++)
			{
				my_searches_chances++;

				uint32_t max_dist = MAX_DIST[i];
				if (!alignments.empty())
					max_dist = std::min(max_dist, alignments[0].dist);

				const auto& read_rev = read.get_rev_compl();
				auto hotspots     = read_processor::find_hotspots(read, mini,
										FRAGMENT_COUNT[i],
										SPOT_SIZE, SPOT_HITS[i], CHECK_COLLI);
				auto hotspots_r_c = read_processor::find_hotspots(read_rev, mini,
										FRAGMENT_COUNT[i],
										SPOT_SIZE, SPOT_HITS[i], CHECK_COLLI);
				read_processor::append_alignments(
					alignments,
					read_processor::start_aligner(
						std::move(hotspots),     read,
						std::move(hotspots_r_c), read_rev,
						env, max_dist,
						MAX_ALIGNS_C[i], MAX_ALIGNS_T[i], MAX_ALIGNS_O));
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
				if (res) mapped++;
				count++;
				if (!(count%1000))
					gedmap_io::flush_row("Searched reads", to_string(count));
				chances_count += my_searches_chances;
			}
		}
	}
	
	fastq_s.close();
	o_s.close();

	gedmap_io::print_row("Searched reads", count);
	tv[tv_MAP].take();
	gedmap_io::print_row("Reads per second per thread", count / tv[tv_MAP].get() / omp_get_max_threads());
	gedmap_io::print_row("Searches", chances_count);
	gedmap_io::print_row("Mapped reads", mapped);
	gedmap_io::print_row("Results", results);
	gedmap_io::dotline();
}


/** @brief 
 * read read from fastq file
 * map read in the index
 * write alignment to sam file
 */
template<class int_type>
void paralell_processor_batch( gedmap_mini::minimizer_index & mini, const std::string & EDS, const adjacency & ADJ,const pos_EDS_to_FA_type & p2FA, std::ifstream & fastq_s, std::ofstream& o_s){
	fastq_s.clear();
	fastq_s.seekg(0);
	uint32_t batchsize = 10000;

	environment<
		gedmap_mini::minimizer_index,
		std::string, // EDS
		adjacency,
		pos_EDS_to_FA_type> env(mini, EDS, ADJ, p2FA);

	gedmap_io::print_row("MAPPING:");
	uint32_t results	= 0;
	uint32_t mapped 	= 0;
	uint32_t count  	= 0;
	uint32_t chances_count	= 0;
	

	if(THREAD_COUNT) omp_set_num_threads(THREAD_COUNT);
	
	bool run = true;
	while(run){
		vector<fasta_read<int_type>> reads;
		reads.reserve(batchsize);
		
		for(uint32_t i = 0; i < batchsize; i++){
			if(fastq_s.good()){
				try{
					fasta_read<int_type> read(fastq_s);
					reads.emplace_back(std::move(read));
				}catch(runtime_error & e){
					//NO READ LEFT
					run = false;
					break;
				}
			}
		}
		
		uint32_t my_searches_chances = 0;	

		std::vector< std::vector<alignment<int_type>> > alignments(reads.size());
		
		if(!MAP_RC){
			#pragma omp parallel for  schedule(dynamic,10)
			for(uint32_t r = 0; r < reads.size(); r++){				
				for (size_t i = 0;
					i < FRAGMENT_COUNT.size() && (alignments[r].empty() || alignments[r][0].dist >= DOUBT_DIST);
					i++)
				{
					my_searches_chances++;
					auto hotspots = read_processor::find_hotspots(reads[r], mini,
										FRAGMENT_COUNT[i],
										SPOT_SIZE, SPOT_HITS[i], CHECK_COLLI);

					uint32_t max_dist = MAX_DIST[i];
					if (!alignments[r].empty())
						max_dist = std::min(max_dist, alignments[r][0].dist - 1);

					read_processor::append_alignments(
						alignments[r],
						read_processor::start_aligner(
							std::move(hotspots),     reads[r],
							env, max_dist,
							MAX_ALIGNS_C[i], MAX_ALIGNS_T[i], MAX_ALIGNS_O));
				}
			}
		}else{
			#pragma omp parallel for  schedule(dynamic,10)
			for(uint32_t r = 0; r < reads.size(); r++){			
				for (size_t i = 0;
					i < FRAGMENT_COUNT.size() && (alignments[r].empty() || alignments[r][0].dist >= DOUBT_DIST);
					i++)
				{
					my_searches_chances++;

					const auto& read_rev = reads[r].get_rev_compl();
					auto hotspots     = read_processor::find_hotspots(reads[r], mini,
											FRAGMENT_COUNT[i],
											SPOT_SIZE, SPOT_HITS[i], CHECK_COLLI);
					auto hotspots_r_c = read_processor::find_hotspots(read_rev, mini,
											FRAGMENT_COUNT[i],
											SPOT_SIZE, SPOT_HITS[i], CHECK_COLLI);

					uint32_t max_dist = MAX_DIST[i];
					if (!alignments[r].empty())
						max_dist = std::min(max_dist, alignments[r][0].dist);

					read_processor::append_alignments(
						alignments[r],
						read_processor::start_aligner(
							std::move(hotspots),     reads[r],
							std::move(hotspots_r_c), read_rev,
							env, max_dist,
							MAX_ALIGNS_C[i], MAX_ALIGNS_T[i], MAX_ALIGNS_O));
				}
			}
		}
		
		for(uint32_t r = 0; r < reads.size(); r++){
			const auto res = read_processor::write_alignment(
				std::move(alignments[r]),
				reads[r],
				o_s,
				WRITE_FAILURE,
				p2FA);
			results += res;
			if (res > 0) mapped++;
			count++;
			if(!(count%1000))
				gedmap_io::flush_row("Searched reads", to_string(count));
			chances_count += my_searches_chances;
		}
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

