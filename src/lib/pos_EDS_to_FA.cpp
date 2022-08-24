#pragma once
using namespace std;
using namespace sdsl;



/**
 * @brief data structure that allowes to map a position of the EDS to original FA position
 * 
 * 
 * ref_ind = indicator: ref_ind[i]=1, iff EDS[i] was taken from FA
 * 
 * chrom_names: names of chromosomes in EDS
 * chrom_starts: chrom_names[i] start at EDS[chrom_start[i]]
 */
struct pos_EDS_to_FA_type{
	typedef sdsl::int_vector<>::size_type	size_type;
	
	
	sdsl::bit_vector 		ref_ind;
	sdsl::rank_support_v5<> rs_ref_ind;	
	
	std::vector<std::string> 	chrom_names;
	sdsl::int_vector<64>		chrom_starts;
	
	
	pos_EDS_to_FA_type();
	pos_EDS_to_FA_type(std::string ref_ind_fname, std::vector<std::string> chrom_names, std::vector<uint64_t> chrom_starts);
	pos_EDS_to_FA_type(sdsl::bit_vector ref_ind, std::vector<std::string> chrom_names, std::vector<uint64_t> chrom_starts);
	
	/** @brief 
	 * <c,p,o> = (pos),
	 * with pos is in chromosom c at position p (with variant offset o)
	 */
	std::tuple<std::string,uint64_t,uint32_t> operator()(uint64_t EDS_pos) const;
	
	bool empty() const;

	size_type 	serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const;
	void 		load(std::istream& in);
};



/** @brief default constructor */
pos_EDS_to_FA_type::pos_EDS_to_FA_type(){};


/** @brief parametrized constructor */
pos_EDS_to_FA_type::pos_EDS_to_FA_type(bit_vector ref_ind, vector<string> chrom_names, vector<uint64_t> chrom_starts)
	:ref_ind(ref_ind), chrom_names(chrom_names){
		
	this->chrom_starts = int_vector<64>(chrom_starts.size(),0,64);
	for(uint32_t i = 0; i < chrom_starts.size(); i++) this->chrom_starts[i] = chrom_starts[i];
	
	sdsl::util::init_support(rs_ref_ind,&ref_ind);
}


/** 
 * @brief mixed constructor, ref_ind from file residul parametrized
 * @param ref_ind_fname file name of ref_ind
 */
pos_EDS_to_FA_type::pos_EDS_to_FA_type(string ref_ind_fname, vector<string> chrom_names, vector<uint64_t> chrom_starts)
	:chrom_names(chrom_names){
		
	this->chrom_starts = int_vector<64>(chrom_starts.size(),0,64);
	for(uint32_t i = 0; i < chrom_starts.size(); i++) this->chrom_starts[i] = chrom_starts[i];

	load_from_file(ref_ind,	ref_ind_fname);
	sdsl::util::init_support(rs_ref_ind,&ref_ind);
};

/** 
 * @brief return position in FA
 * @return <c,p,o> , EDS[i] is in chromosom c at position p (with variant offset o)
 * if EDS[i] is in a variant scope, then EDS[i] does not occour in FA. 
 * The previos position that occours in FA is EDS[i-o]
 */
tuple<string,uint64_t,uint32_t> pos_EDS_to_FA_type::operator()(const uint64_t EDS_pos) const{
	
	
	uint32_t chrom_id = 0;
	while( chrom_id+1 < chrom_starts.size() && chrom_starts[chrom_id+1] <= EDS_pos )
		chrom_id++;
	
	uint32_t off = 0;
	for(;!ref_ind[EDS_pos-off];off++);
	
	
	uint64_t fa_symbols = rs_ref_ind(EDS_pos);
	if(chrom_id > 0) fa_symbols -= rs_ref_ind(chrom_starts[chrom_id]);

	
	return make_tuple( chrom_names[chrom_id], fa_symbols ,off);
}




/**
 * @brief converts a string to a int_vector<8>
 * for serialize
 */
int_vector<8> string_to_iv(const string & s){
	int_vector<8> iv = int_vector<8>(s.size(),0,8);
	for(uint32_t i = 0; i < s.size(); i++) iv[i] = s[i];
	return iv;
}

/**
 * @brief converts a int_vector<8> to a string
 * for load
 */
string iv_to_string(const int_vector<8> & iv){
	string s(iv.size(), '.');
	for(unsigned int i = 0; i < iv.size(); i++) s[i] = iv[i];
	return s;
}


/**
 * @brief converts a vector<int_vector<8>> to a vector<string>
 * for load
 */
vector<string> ivv_to_string_vector(const vector<int_vector<8>> & ivv){
	vector<string> sv = vector<string>(ivv.size());
	for(unsigned int i = 0; i < ivv.size(); i++)
		sv[i] = iv_to_string(ivv[i]);
	return sv;
}

/**
 * @brief converts a vector<string> to a vector<int_vector<8>>
 * for serialize
 */
vector<int_vector<8>> string_vector_to_ivv(const vector<string> & sv){
	vector<int_vector<8>> ivv = vector<int_vector<8>>(sv.size());
	for(unsigned int i = 0; i < sv.size(); i++)
		ivv[i] = string_to_iv(sv[i]);
	return ivv;
}

/**
 * @brief serialize pos_EDS_to_FA_type with sdsl methods
 */
pos_EDS_to_FA_type::size_type pos_EDS_to_FA_type::serialize(ostream& out, structure_tree_node* v,string name) const{
	
	
	structure_tree_node *child 	= structure_tree::add_child(v, name, sdsl::util::class_name(*this));
	
	pos_EDS_to_FA_type::size_type written_bytes 	= 0;
	
	
	written_bytes += ref_ind	.serialize(out, child, "ref_ind");
	written_bytes += rs_ref_ind	.serialize(out, child, "rank_support");	
	
	written_bytes += chrom_starts	.serialize(out, child, "chrom_starts");
	vector<int_vector<8>> c_names = string_vector_to_ivv(chrom_names);
	written_bytes += serialize_vector<int_vector<8>>(c_names,out,child,"chrom_names");
	

	structure_tree::add_size(child, written_bytes);
	return written_bytes;
}

bool pos_EDS_to_FA_type::empty() const{
	return (ref_ind.size() == 0);
}
/**
 * @brief loads pos_EDS_to_FA_type with sdsl methods
 */
void pos_EDS_to_FA_type::load(std::istream& in){
	ref_ind		.load(in);
	rs_ref_ind		.load(in,&ref_ind);
	chrom_starts	.load(in);
	
	vector<int_vector<8>> c_names = vector<int_vector<8>>(chrom_starts.size());		
	load_vector<int_vector<8>>(c_names,in);
	chrom_names = ivv_to_string_vector(c_names);	
}

/** \brief for debugging */
// std::ostream& operator<< (std::ostream& os, const pos_EDS_to_FA_type& t3){
// 	os << "pos_EDS_to_FA_type:" << endl;
// 	os << t3.ref_ind << endl;
// 	os << t3.chrom_starts << endl;
// 	os << t3.chrom_names << endl;
// 	return os;	
// }
