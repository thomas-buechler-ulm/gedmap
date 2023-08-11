#pragma once
#include "pos_EDS_to_FA.cpp"
using namespace std;
using namespace sdsl;




/**
 * @brief data structure that allowes to map a position of the EDS to the original GFA node and position
 * 
 * 
 * node_start_ind = indicator: node_start_ind[i]=1, iff EDS[i] is the start of a node in the original GFA file
 * 

 */
struct pos_EDS_to_GFA_type: public Transform{
	typedef sdsl::int_vector<>::size_type	size_type;
	

	//ENTRY PER NODE
	sdsl::int_vector<0> node_number; //NUMBER IN GFA FILE
	sdsl::int_vector<1> orientation;
	sdsl::int_vector<0> seq_name_idx; //SEQ IDX IN SEQ_NAME
	sdsl::int_vector<0> offset;       //SEQ OFFSET OF NODE

	//ONCE
	vector<string>    seq_name;

	//ENTRY PER CHAR IN EDSG
	sdsl::bit_vector 		node_start_ind;
	sdsl::rank_support_v5<> 	rs_node_start_ind;
	//sdsl::... ss_node_start_ind; //TODO
	
	
	
	/** @brief default constructor */
	pos_EDS_to_GFA_type(){};

	/** @brief parametrized constructor */
	pos_EDS_to_GFA_type(
		sdsl::int_vector<0> node_number,
		sdsl::int_vector<1> orientation,
		sdsl::int_vector<0> seq_name_idx,
		sdsl::int_vector<0> offset,
		vector<string>  seq_name,
		sdsl::bit_vector node_start_ind
	):node_number(node_number),orientation(orientation),seq_name_idx(seq_name_idx),offset(offset),seq_name(seq_name),node_start_ind(node_start_ind){
		sdsl::util::init_support(rs_node_start_ind,&node_start_ind);
	};
	
	/** @brief 
	 * <u,p,o> = (pos),
	 * with pos is in node u at position p (with orientation o)
	 */
	std::tuple<std::string,uint64_t,std::string> operator()(uint64_t EDS_pos) const;
	
	bool empty() const;

	size_type 	serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const;
	void 		load(std::istream& in);
};



/** 
 * @brief return position in FFA
 * @return <c,p,n> ,
 * c sequence / node;
 * p = position
 * n = offset???
 */
tuple<std::string,uint64_t,std::string> pos_EDS_to_GFA_type::operator()(uint64_t EDS_pos) const{
	EDS_pos--;

	//cout << "operator(" << EDS_pos << ")" << endl;
	
	if(EDS_pos >= node_start_ind.size()) throw runtime_error("TRANSFORM: EDS_POS TO BIG");

	uint64_t p = 0;
	while( !node_start_ind[EDS_pos - p] ){
		p++;
		if(EDS_pos < p) throw runtime_error("TRANSFORM: WEIRD");
	};


	//cout << "p=" << p << endl;

	size_t node = rs_node_start_ind(EDS_pos-p);
	//cout << "node=" << node << endl;
	if(seq_name_idx.size() <= node) throw runtime_error("TRANSFORM: WEIRD2");
	if(seq_name.size() <= seq_name_idx[node]) throw runtime_error("TRANSFORM: WEIRD3");
	string seq = seq_name[ seq_name_idx[node] ];


	//cout << "seq=" << seq << endl;
	//cout << "offset=" << offset[node] << endl;

	if(seq != "" && !orientation[node])	return make_tuple( seq, offset[node]+p, "" );
	if(seq_name_idx.size() <= node) throw runtime_error("TRANSFORM: WEIRD4");
	if(seq == "" && !orientation[node])                  return make_tuple("Node:"+to_string(node_number[node]), p, "");
	return make_tuple("Node_reverse:"+to_string(node_number[node]), p, "");
}


sdsl::int_vector<32> flatten(const std::vector<std::string> & v){
	size_t c = 1 + v.size() + 1;
	for(auto s : v) c+= s.size();
	sdsl::int_vector<32> out(c,0,32);


	out[0] = v.size();
	out[1] = 1 + v.size() + 1;

	for(size_t i = 0; i < v.size(); i++){
		out[i+2] = out[i+1];
		for(auto c : v[i]) out[ out[i+2]++ ] = c;
	}
	return out;
}

std::vector<std::string> deFlatten(const sdsl::int_vector<32> & iv ){
	vector<string> out(iv[0]);
	for(size_t i = 0; i < out.size(); i++){
		out[i] = "";
		for(size_t j = iv[i+1]; j < iv[i+2]; j++)
			out[i] += (char) iv[j];
	}
	return out;
}



/**
 * @brief serialize pos_EDS_to_FA_type with sdsl methods
 */
pos_EDS_to_GFA_type::size_type pos_EDS_to_GFA_type::serialize(ostream& out, structure_tree_node* v,string name) const{
	
	
	structure_tree_node *child 	= structure_tree::add_child(v, name, sdsl::util::class_name(*this));
	
	pos_EDS_to_GFA_type::size_type written_bytes 	= 0;

	written_bytes += node_number	.serialize(out, child, "node_number");
	written_bytes += orientation	.serialize(out, child, "orientation");
	written_bytes += seq_name_idx	.serialize(out, child, "seq_name_idx");
	written_bytes += offset	.serialize(out, child, "offset");
	written_bytes += node_start_ind   .serialize(out, child, "node_start_ind");
	written_bytes += rs_node_start_ind.serialize(out, child, "rs_node_start_ind");

	
	sdsl::int_vector<32> seqs = flatten(seq_name);
	written_bytes += seqs.serialize(out, child, "seq_names");

	structure_tree::add_size(child, written_bytes);
	return written_bytes;
}


bool pos_EDS_to_GFA_type::empty() const{
	return (node_number.size() == 0);
}
/**
 * @brief loads pos_EDS_to_FA_type with sdsl methods
 */
void pos_EDS_to_GFA_type::load(std::istream& in){
	node_number		.load(in);
	orientation		.load(in);
	seq_name_idx		.load(in);
	offset			.load(in);
	node_start_ind		.load(in);
	rs_node_start_ind	.load(in,&node_start_ind);
	sdsl::int_vector<32> seqs;
	seqs.load(in);
	seq_name = deFlatten(seqs);
}



std::ostream& operator<< (std::ostream& os, const pos_EDS_to_GFA_type& T){
	os << "pos_EDS_to_GFA_type:" << endl;
	os << "node_number " << T.node_number << endl;
	os << "orientation " << T.orientation << endl;
	os << "seq_name_idx " << T.seq_name_idx << endl;
	os << "offset " << T.offset << endl;
	os << "node_start_ind " << T.node_start_ind << endl;
	os << "seqs: " << flush;
	for(auto s : T.seq_name) cout << "\""<< s << "\""<< endl;
	os <<  endl;
	return os;
}
