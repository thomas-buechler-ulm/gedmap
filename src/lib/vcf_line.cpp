#pragma once
#include <string>
#include <sstream>
#include <vector>

using namespace std;

struct vcf_line{
	//data
	std::string	chrom;
	size_t		pos;
	std::string	id;
	std::string	ref;
	std::string	alt;
	vector<string> alts;
	std::string	info;

	uint32_t TYPE;
	static const uint32_t SNP_TYPE = 1; //SNP
	static const uint32_t SMALL_ALT_TYPE = 2; //SMALL PLAIN VARIANTS
	static const uint32_t BIG_ALT_TYPE = 4; //BIG PLAIN VARIANTS
	static const uint32_t CNV_TYPE = 8; //CNV
	static const uint32_t IN_TYPE = 16; //PLAIN INSERTIONS
	static const uint32_t DEL_TYPE = 32; //PLAIN DELETIONS
	static const uint32_t STRUCT_TYPE = 64; // STRUCTURAL OR BIG





	/** \brief calculates the vector of alternatives */
	void init(){
		TYPE = 0;

		if (alt[0] == '<'){
			TYPE |= STRUCT_TYPE;
			if( alt.find("<CN") !=  string::npos ) TYPE |= CNV_TYPE;
			return;
		}

		//PARSE ALTERNATIVES
		alts = vector<string> (0);
		alts.push_back(ref);
		uint32_t alt_begin 	= 0;
		uint32_t alt_length	= 0;
		for(uint32_t i = 0; i < alt.size(); i++)
			if(alt[i] != ',')
				alt_length++;
		else{    //if "," add substring to alt
			alts.push_back(alt.substr(alt_begin,alt_length));
			alt_begin 	= i+1;
			alt_length 	= 0;
		}
		alts.push_back(alt.substr(alt_begin,alt_length));

		for( string a : alts){ //TODO SANITY CHECK
			if (a.size() < 1) throw std::runtime_error("ERROR IN PARSING ALT STRING >" + alt + "<");
		}

		//CHECK IF SNP
		bool is_snp = true;
		for( string a : alts) if(a.size() != 1) is_snp = false;
		if(is_snp){
			TYPE |= SNP_TYPE;
			TYPE |= SMALL_ALT_TYPE;
			return;
		}

		if(alts.size() == 2 && alts[0].size() == 1 && alts[1][0] == alts[0][0]  ) TYPE |= IN_TYPE;
		if(alts.size() == 2 && alts[1].size() == 1 && alts[1][0] == alts[0][0]  ) TYPE |= DEL_TYPE;


		for(auto allele : alts){
			if(allele.size() > gedmap_parse::PLAIN_ALT_LIMIT) {
				TYPE |= BIG_ALT_TYPE | STRUCT_TYPE;
				return;
			}
		}
		TYPE |= SMALL_ALT_TYPE;
	}

	vcf_line(){};
	vcf_line(const string & line){
		stringstream stream_line(line);
		string pos_s;
		stream_line >> chrom >> pos_s >> id >> ref >> alt >> info >> info >> info; //ignore qual and filter
		pos = stol(pos_s);
		init();
	}



	/** returns the eds */
	pair<string,vector<bool>> get_eds() const{
		if(TYPE & SNP_TYPE){
			string eds   (2*alts.size() +1, '(' );
			for(uint32_t i = 0; i < alts.size(); i++ )	{
				eds[2*i +1] = alts[i][0];
				eds[2*i +2] = '|';
			}
			eds.back()= ')';
			vector<bool> I(2*alts.size() +1,false);
			I[1] = true;
			return make_pair(eds,I);
		}

		if((TYPE & IN_TYPE) && (TYPE & SMALL_ALT_TYPE)){
			string eds = ref + "(|" + alt.substr(1) + ")";
			vector<bool> I(alt.size() +3, false);
			I[0] = true;
			return make_pair(eds,I);
		}

		if((TYPE & DEL_TYPE) && (TYPE & SMALL_ALT_TYPE)){
			string eds = ref.substr(0,1) + "(|" + ref.substr(1) + ")";
			vector<bool> I(alts[0].size() +3, true);
			I[1] = false;
			I[2] = false;
			I.back() = false;
			return make_pair(eds,I);

		}

		//OTHER
		if(TYPE & SMALL_ALT_TYPE){
			string eds = "(" + ref;
			for(uint32_t i = 1; i < alts.size(); i++)
				eds += "|" + alts[i];
			eds += ")";

			vector<bool> I(eds.size(),false);
			for(uint32_t i = 0; i < ref.size(); i++) I[i+1] = true;
			return make_pair(eds,I);
		}

		throw std::runtime_error("ERROR IN PARSING ALT TO EDS: INVALIDE TYPE " + to_string(TYPE));
	};


	typedef pair <size_t,size_t> edge_type;
	typedef tuple<size_t,string,uint32_t> node_split;

	size_t sv_end() const{
		size_t found1 = 0;
		while(true){
			found1 = info.find("END=",found1);
			if(found1 == 0 || info[found1-1] == ';') break;
			if(found1 == string::npos) throw std::runtime_error("ERROR IN PARSING VARIANT AT POS: " + to_string(pos));
			found1++;
		}
		found1 += 4;
		size_t found2 = info.find(";",found1);
		return stol(info.substr(found1,found2-found1));
	};

	const static uint32_t COPY_0 = 1;
	const static uint32_t COPY_2 = 2;
	pair< vector<node_split>, vector<edge_type> > get_pos_and_edges() {
		if( TYPE & CNV_TYPE) {
			std::size_t found1 = 0 , found2 = 0;
			uint32_t cn = 0;
			while(true){
				//alt has the form <CNx>,<CNy>,... with integers x,y.
				found1 = alt.find("<CN", found2);
				if(found1 == string::npos) break;
				found1 +=3;
				found2 = alt.find(">",found1);
				uint32_t x = stoi(alt.substr(found1,found2-found1));
				if(x == 0) cn |= COPY_0;
				if(x >= 2) cn |= COPY_2;
			}
			size_t end = sv_end();
			vector<edge_type> edges;
			vector<node_split> positions { make_tuple(pos,"#",CNV_TYPE), make_tuple(end, "#",CNV_TYPE)};
			edges.push_back( make_pair(pos,pos));
			edges.push_back( make_pair(end,end));
			if(cn & COPY_0) edges.push_back( make_pair(pos,end) );
			if(cn & COPY_2) edges.push_back( make_pair(end,pos) );
			return make_pair(positions,edges);
		}


		if( TYPE & DEL_TYPE){
			size_t begin = pos + 1;
			size_t end   = pos + ref.size();
			vector<node_split> positions { make_tuple(begin,"#",DEL_TYPE), make_tuple(end, "#",DEL_TYPE)};
			vector<edge_type> edges { make_pair(begin,begin), make_pair(begin,end), make_pair(end,end) };
			return make_pair(positions,edges);
		}

		if( TYPE & IN_TYPE){
			size_t begin = pos + 1;
			vector<node_split> positions { make_tuple(begin,"#"+alts[1].substr(1)+"#",IN_TYPE)};
			vector<edge_type> edges { make_pair(begin,begin)};
			return make_pair(positions,edges);
		}

		if( TYPE & BIG_ALT_TYPE){
			size_t end = pos + ref.size();
			string new_nodes = "#";
			for (uint32_t i = 1; i < alts.size(); i++) new_nodes += alts[i] + '#';
			vector<node_split> positions { make_tuple(pos,"#",BIG_ALT_TYPE), make_tuple(end, new_nodes,BIG_ALT_TYPE)};
			vector<edge_type> edges { make_pair(pos,pos), make_pair(pos,end), make_pair(end,end) };
			return make_pair(positions,edges);
		}
		throw std::runtime_error("ERROR IN PARSING ALT TO NODES AND EDGES: INVALIDE TYPE " + to_string(TYPE));
	};

};


