#pragma once
#include "vcf_line.cpp"

struct GEDS_builder{

	fa_stream fa;
	vcf_stream vcf;

	//ALL SEQENCES THAT ARE IN THE EDS
	vector<string> seq_names;
	vector<size_t> seq_start;

	//TEMPORARY FILES
	string fn_tmp_eds, fn_tmp_geds, fn_tmp_I, fn_tmp_II, fn_tmp_V;

	std::ofstream ofs_EDS;
	sdsl::int_vector_buffer<1> I_EDS; //indicater if EDS[i] stems from ref
	sdsl::int_vector_buffer<1> I_GEDS; //indicater if EDS[i] stems from ref
	sdsl::int_vector_buffer<1> V_EDS; //indicater if EDS[i] is inside alternative scope

	//TEMPORARY EDGE INFORMATION
	typedef pair <size_t,size_t> edge_type;
	typedef tuple<size_t,string,uint32_t> node_split;
	vector<std::set< node_split >>  insert_positions; // for each sequnce, the set of positions a string of the form: #( seq#  )*
	vector<vector<edge_type>> EDS_edges; // for each sequence, the set of edges to insert
	vector<edge_type> GEDS_edges; // the set of edges in the GEDS





	uint32_t lines;
	std::vector<uint32_t> variant_stats = std::vector<uint32_t>(9,0);
	static const uint32_t TOTAL		= 0;
	static const uint32_t SNP		= 1;
	static const uint32_t SMALL_ALT	= 2;
	static const uint32_t CNV		= 3;
	static const uint32_t SV_INDEL	= 4;
	static const uint32_t BIG_ALT		= 5;
	static const uint32_t IGNORED_OVERLAP	= 6;
	static const uint32_t UNSUPPORTED	= 7;
	static const uint32_t SMALL_INDEL	= 8;

	void parse2EDS();
	void parse2GEDS();

	GEDS_builder(string out_path){
		string fname = out_path.substr(out_path.find_last_of("/") + 1);
		fn_tmp_eds  = TMP_DIR + "/gedmap_tmp_EDS_"  + fname;
		fn_tmp_geds = out_path;
		fn_tmp_I    = TMP_DIR + "/gedmap_tmp_I_"    + fname;
		fn_tmp_II   = TMP_DIR + "/gedmap_tmp_II_"   + fname;
		fn_tmp_V    = TMP_DIR + "/gedmap_tmp_V_"    + fname;
	}

	void build(std::string fname_fa, std::string fname_vcf){

		sys_timer t;
		t.start();

		gedmap_io::print_row("Read VCF and FA and include eds-type variations");
		fa.open(fname_fa);
		vcf.open(fname_vcf);
		parse2EDS();
		gedmap_io::dotline();
		gedmap_io::print_row("Time for parsing short vars:", t.stop_and_get(), " s");
		gedmap_io::dotline();
		gedmap_io::print_row("Total number of included variants",	variant_stats[TOTAL]);
		gedmap_io::print_row("- SNPs",	variant_stats[SNP]);
		gedmap_io::print_row("- small indels",	variant_stats[SMALL_INDEL]);
		gedmap_io::print_row("- small alternatives",	variant_stats[SMALL_ALT]);
		gedmap_io::print_row("- copy numbers",	variant_stats[CNV]);
		gedmap_io::print_row("- large indels",	variant_stats[SV_INDEL]);
		gedmap_io::print_row("- large alternatives",	variant_stats[BIG_ALT]);
		gedmap_io::print_row("- ignored overlapping variants",	variant_stats[IGNORED_OVERLAP]);
		gedmap_io::print_row("- unsupported variants",	variant_stats[UNSUPPORTED]);
		gedmap_io::dotline();

		if(!gedmap_parse::INCLUDE_CNV){
			gedmap_io::file_move(fn_tmp_eds,fn_tmp_geds);
		}
		else{
			t = sys_timer();
			t.start();
			gedmap_io::print_row("Include sturctural variations");
			parse2GEDS();
			gedmap_io::dotline();
			gedmap_io::print_row("Time for including SV:", t.stop_and_get(), " s");;
		}


		sdsl::store_to_file(pos_EDS_to_FA_type(fn_tmp_II, seq_names, seq_start), fn_tmp_geds + "." + FEX_2FA);
		sdsl::store_to_file(adjacency(GEDS_edges)     , fn_tmp_geds + "." + FEX_ADJ);

		//CLEANUP
		gedmap_io::file_remove(fn_tmp_eds);
		//gedmap_io::file_remove(fn_tmp_geds);
		gedmap_io::file_remove(fn_tmp_I);
		gedmap_io::file_remove(fn_tmp_II);
		gedmap_io::file_remove(fn_tmp_V);
	}
};

/**
 * adds die bit @val to the bit_vector_buffer @bvb until its size is @p
 */
void fill_bvb_pos(sdsl::int_vector_buffer<1>& bvb, size_t p,bool val){
	while(bvb.size() < p) bvb.push_back(val);
}

void GEDS_builder::parse2EDS(){ //TODO LIMIT
	ofs_EDS.open(fn_tmp_eds,ios_base::out|ios_base::trunc);
	I_EDS = sdsl::int_vector_buffer<1> (fn_tmp_I, ios::out|ios_base::trunc);
	V_EDS = sdsl::int_vector_buffer<1> (fn_tmp_V, ios::out|ios_base::trunc);

	seq_names.push_back(fa.seq_name);
	seq_start.push_back(ofs_EDS.tellp());
	EDS_edges.push_back(vector<edge_type>());
	insert_positions.push_back( set<node_split>() );


	std::map<string,uint32_t> unsupported_var_type_count;

	vcf_line vcfl;

	while( vcf.get(vcfl) ){
		while(vcfl.chrom != fa.seq_name){
			fa.copy_to_ofs_until_next_seq(ofs_EDS); //IF NEEDED GO TO NEXT SEQUENCE IN FA
			seq_names.push_back(fa.seq_name);
			seq_start.push_back(ofs_EDS.tellp());
			EDS_edges.push_back(vector<edge_type>());
			insert_positions.push_back( set<node_split>() );
			fill_bvb_pos(I_EDS,(size_t) ofs_EDS.tellp()-1,1);
			fill_bvb_pos(I_EDS,ofs_EDS.tellp(),0);
			fill_bvb_pos(V_EDS,ofs_EDS.tellp(),0);
		}


		if( (vcfl.TYPE & vcf_line::BIG_ALT_TYPE || vcfl.TYPE & vcf_line::CNV_TYPE) && gedmap_parse::INCLUDE_CNV){
			std::pair<std::vector<node_split>,std::vector<edge_type>> nodes_and_edges = vcfl.get_pos_and_edges();
			for(auto ns : nodes_and_edges.first)  insert_positions.back().insert(ns);
			for(auto e  : nodes_and_edges.second) EDS_edges.back().push_back(e);

			if( vcfl.TYPE & vcf_line::CNV_TYPE)						variant_stats[GEDS_builder::CNV]++;
			else if (vcfl.TYPE & (vcf_line::IN_TYPE | vcf_line::DEL_TYPE) )	variant_stats[GEDS_builder::SV_INDEL]++;
			else  												variant_stats[GEDS_builder::BIG_ALT]++;
			variant_stats[TOTAL]++;

			continue;
		}


		if( fa.pos >= vcfl.pos){ //OVERLAP
			variant_stats[IGNORED_OVERLAP]++;
			continue;
		}

		fa.copy_to_ofs_until_t(vcfl.pos-1,ofs_EDS);
		fill_bvb_pos(I_EDS,ofs_EDS.tellp(),1);
		fill_bvb_pos(V_EDS,ofs_EDS.tellp(),0);

		if( vcfl.TYPE & vcf_line::SMALL_ALT_TYPE ){ // VARIATION IN EDS FORMAT
			string veds;
			vector<bool> I_VAR;
			std::tie(veds,I_VAR) = vcfl.get_eds();
			fa.skip(vcfl.ref.size());
			ofs_EDS << veds;
			for(bool b : I_VAR) I_EDS.push_back(b);
			V_EDS.push_back( ! (vcfl.TYPE | vcf_line::IN_TYPE)  );
			fill_bvb_pos(V_EDS,ofs_EDS.tellp(),1);

			if( vcfl.TYPE & vcf_line::SNP_TYPE)						variant_stats[SNP]++;
			else if(vcfl.TYPE & (vcf_line::IN_TYPE | vcf_line::DEL_TYPE))	variant_stats[GEDS_builder::SMALL_INDEL]++;
			else 												variant_stats[GEDS_builder::SMALL_ALT]++;
			variant_stats[TOTAL]++;

			continue;
		}

		variant_stats[ GEDS_builder::UNSUPPORTED]++;

		//COUNT UNSUPPOERTED TYPES
		auto it = unsupported_var_type_count.find(vcfl.alt);
		if (it == unsupported_var_type_count.end())	unsupported_var_type_count[vcfl.alt] = 1;
		else	it->second++;
	}

	for (auto var_count : unsupported_var_type_count)
		gedmap_io::print_error( "VCF ENTRY TYPE " + (var_count.first.size()>20?"":string(20-var_count.first.size(), ' ')) + var_count.first + " NOT SUPPORTED (ignored " + to_string(var_count.second) + " times)");

	while(true){ //COPY TILL END
		try{
			fa.copy_to_ofs_until_next_seq(ofs_EDS);
			seq_names.push_back(fa.seq_name);
			seq_start.push_back(ofs_EDS.tellp());
		}catch(invalid_argument &ia){break;} //NO NEXT SEQ
	}
	fill_bvb_pos(I_EDS,ofs_EDS.tellp(),1);
	fill_bvb_pos(V_EDS,ofs_EDS.tellp(),0);

	fa.close();
	vcf.close();
	ofs_EDS.close();
	I_EDS.close();
	V_EDS.close();
}

void GEDS_builder::parse2GEDS(){
	ifstream EDS_in;
	ofstream GEDS_out;
	EDS_in  .open(fn_tmp_eds ,ios_base::in);
	GEDS_out.open(fn_tmp_geds,ios_base::out|ios_base::trunc);
	I_GEDS = sdsl::int_vector_buffer<1> (fn_tmp_II, ios::out|ios_base::trunc);
	sdsl::bit_vector EDS_I;
	sdsl::bit_vector EDS_V;
	sdsl::load_from_file(EDS_I, fn_tmp_I);
	sdsl::load_from_file(EDS_V, fn_tmp_V);
	sdsl::select_support_mcl<1,1>  I_ss;
	sdsl::rank_support_v5<> I_rs;
	sdsl::util::init_support(I_ss,&EDS_I);
	sdsl::util::init_support(I_rs,&EDS_I);


	stream_to_stream_copy EDS_TO_GEDS(&EDS_in, &GEDS_out, 100*1024); //100 KB int_vector_buffer
	bv_to_bv_copy I_EDS_TO_I_GEDS(&EDS_I, &I_GEDS); //100 KB int_vector_buffer


	//WRITE GEDS
	for( uint32_t i = 0; i < seq_names.size(); i++){
		gedmap_io::flush_row("Add node boundaries and edge information to seq ", seq_names[i]);
		//COPY TILL SEQENCE STARTSls
		size_t EDS_seq_start = seq_start[i];
		size_t chromosom_offset = I_rs(EDS_seq_start) ;
		I_EDS_TO_I_GEDS.copy(EDS_seq_start - EDS_in.tellg());
		EDS_TO_GEDS	.copy(EDS_seq_start - EDS_in.tellg());
		seq_start[i] = GEDS_out.tellp();


		//MAP (seq,pos) -> (<GEDS_POS>)
		std::map< size_t, vector<size_t> > map_to_GEDS_pos;
		for( auto ip : insert_positions[i]){

			size_t split_pos = I_ss(   get<0>(ip) + chromosom_offset );
			// SHIFT SPLIT OUTSIDE THE BRACKETS
			while(EDS_V[split_pos]) split_pos++;

			//COPY TO SPLITPOS
			I_EDS_TO_I_GEDS.copy(split_pos - EDS_in.tellg());
			EDS_TO_GEDS	.copy(split_pos - EDS_in.tellg());

			//POSITIONS OF NEW BOUNDARY SYMBOLS
			vector<size_t> boundary_pos {get<2>(ip)};
			for(uint32_t i = 0; i < get<1>(ip).size(); i++ )
				if(get<1>(ip)[i] == '#')
					boundary_pos.push_back((size_t) GEDS_out.tellp() + i);
			map_to_GEDS_pos[ get<0>(ip) ] = boundary_pos;

			//WRITE NEW BOUNDARY SYMBOLS
			GEDS_out << get<1>(ip);
			for(uint32_t j = 0; j < get<1>(ip).size(); j++) I_GEDS.push_back(0);
		}

		//ADJUST EDGES
		for( auto e : EDS_edges[i]){
			vector<size_t> start = map_to_GEDS_pos[e.first];
			vector<size_t> end = map_to_GEDS_pos[e.second];
			if(start.size() < 2 || end.size() < 2)
				throw std::runtime_error("ERROR IN ADJUSTING EDGES, EDGE NOT LEADED CORRECTELY)");

			if( (start[0] == vcf_line::CNV_TYPE && end[0] == vcf_line::CNV_TYPE) || (start[0] == vcf_line::DEL_TYPE && end[0] == vcf_line::DEL_TYPE)){

				GEDS_edges.push_back( make_pair(start[1], end[1]));

			} else if (start[0] == vcf_line::IN_TYPE  && e.first == e.second){

				GEDS_edges.push_back( make_pair(start[1], start[1]));
				GEDS_edges.push_back( make_pair(start[2], start[2]));
				GEDS_edges.push_back( make_pair(start[1], start[2]));

			} else if (start[0] == vcf_line::BIG_ALT_TYPE && end[0] == vcf_line::BIG_ALT_TYPE){

				if (start.size()== 2 && end.size()==2) //READ REF
					GEDS_edges.push_back( make_pair(start[1], end[1]));
				else if (start.size() == 2) //SKIP REF
					for(uint32_t j = 1; j < end.size()-1 ; j++) GEDS_edges.push_back( make_pair(start[1], end[j]));
				else if (e.first == e.second) // JUMP TO END OF VARIANT
					for(uint32_t j = 1; j < start.size() ; j++) GEDS_edges.push_back( make_pair(start[j], start.back()));
				else{
					throw std::runtime_error("ERROR IN ADJUSTING EDGES, SUCH AN EDGE SHOULD NOT EXIST (BIG_ALT_TYPE)");
				}

			}  else {
				throw std::runtime_error("ERROR IN ADJUSTING EDGES, SUCH AN EDGE SHOULD NOT EXIST (TYPE " + to_string(start[0]) + "-"+ to_string(end[0]) +")");
			}
		}
		EDS_edges[i].resize(0);
	}
	EDS_TO_GEDS	.copy_all();
	I_EDS_TO_I_GEDS.copy_all();
	EDS_in.close();
	GEDS_out.close();
	I_GEDS.close();
}
