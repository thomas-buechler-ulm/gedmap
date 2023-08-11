#pragma once
#include <map>

using namespace std;
using namespace sdsl;

/** @brief Implements an adjaceny list on the positions of '#'-Symbols */
struct adjacency{
	typedef sdsl::int_vector<>::size_type	size_type;
	typedef pair <uint64_t,uint64_t> edge_type;
	const static bool FORWARD = false;
	const static bool BACKWARD = true;
	bool initialised;
	
	// let A[i] be the set of adjacent positions of i
	// targets is the concatenation A[0]A[1]A[2]...
	// forward_edges is a map that maps i to the beginning of A[i] in the concatenation targets
	// so we need only to 1D vectors instead of a 2D vector
	// the problem of a 2D vector is that the i-values are positions in the EDS and therefore cannot be used as indices

	std::map<uint64_t,uint64_t> forward_edges;
	std::vector<uint64_t> forward_targets;


	std::map<uint64_t,uint64_t> backward_edges;
	std::vector<uint64_t> backward_targets;
	
	adjacency(){initialised=false;};
	adjacency(std::vector<edge_type> edges);


	using result_type = it_range<std::vector<uint64_t>::const_iterator>;
	result_type operator()(uint64_t node, bool direction = FORWARD) const {
		//if(!initialised) throw runtime_error("adjacency: not initialised."); RETURN EMPTY VECTOR WHEN ADJ NOT GIVEN
		if(direction == FORWARD)
			return get_tagets(node,forward_edges,forward_targets);
		else
			return get_tagets(node,backward_edges,backward_targets);
	};
	
	
	size_type 	serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const;
	void 		load(std::istream& in);
	
	
private:
	/** @brief returns the set of adjacent positions of '#' */
	result_type get_tagets(uint64_t node, const std::map<uint64_t,uint64_t> & nodes, const std::vector<uint64_t> & targets ) const {
		std::vector<uint64_t>::const_iterator first;
		std::vector<uint64_t>::const_iterator last;
		
		//find first element of A[node] in targets
		std::map<uint64_t,uint64_t>::const_iterator n = nodes.find(node);
		if(n == nodes.cend()) return result_type();//throw runtime_error("adjacency: node is not in adjacency list");
		first = targets.cbegin() + n->second;
		
		//find last element of A[node] in targets
		if(++n != nodes.cend())
			last = targets.cbegin() + n->second;
		else 
			last = targets.cend();


		return result_type(first, last);
	};
};

/** 
 */
adjacency::adjacency(std::vector<edge_type> edges){
	std::map<uint64_t,uint64_t>::iterator hint;
	
	//FORWARD DIRECTION
	std::sort (edges.begin(), edges.end(), [](edge_type a, edge_type b) {return (a.first < b.first || (a.first == b.first && a.second < b.second)); }  );
	forward_targets = std::vector<uint64_t>(edges.size());
	forward_edges = std::map<uint64_t,uint64_t>();
	hint = forward_edges.begin();


	for(size_t i = 0; i < edges.size(); i++){
		forward_targets[i] = edges[i].second;
		if(i == 0 || edges[i].first != edges[i-1].first )
			hint = forward_edges.insert( hint, pair<uint64_t,uint64_t> (edges[i].first, i) );
	}
	
	//BACKWARD DIRECTION
	std::sort (edges.begin(), edges.end(), [](edge_type a, edge_type b) {return (a.second < b.second || (a.second == b.second && a.first < b.first)); }  );
	
	backward_targets = std::vector<uint64_t>(edges.size());
	backward_edges = std::map<uint64_t,uint64_t>();
	hint = backward_edges.begin();
	
	for(size_t i = 0; i < edges.size(); i++){
		backward_targets[i] = edges[i].first;
		if(i == 0 || edges[i].second != edges[i-1].second )
			hint = backward_edges.insert( hint, pair<uint64_t,uint64_t> (edges[i].second, i) );
	}
	initialised = true;
};

	
std::vector<uint64_t> map2vec(std::map<uint64_t,uint64_t> m){
	std::vector<uint64_t> v  = std::vector<uint64_t>(2*m.size());
	uint64_t i = 0;
	for (std::map<uint64_t,uint64_t>::iterator it=m.begin(); it!=m.end(); ++it){
		 v[i++] = it->first;
		 v[i++] = it->second;
	}
	return v;
};


std::map<uint64_t,uint64_t> vec2map(std::vector<uint64_t> v){
	std::map<uint64_t,uint64_t> m  = std::map<uint64_t,uint64_t>();
	std::map<uint64_t,uint64_t>::iterator hint = m.begin();
	for(size_t i = 0; i < v.size(); i+=2)
		hint = m.insert( hint, pair<uint64_t,uint64_t> (v[i], v[i+1]) );
	return m;
};

/**
 * @brief serialize adjacency with sdsl methods
 */
adjacency::size_type adjacency::serialize(ostream& out, structure_tree_node* v,string name) const{
	
	
	structure_tree_node *child 	= structure_tree::add_child(v, name, sdsl::util::class_name(*this));
	
	adjacency::size_type written_bytes 	= 0;

	if(forward_targets.size() != backward_targets.size()) throw runtime_error("edge sizes not equal");


	sdsl::int_vector<64> sizes = sdsl::int_vector<64>(3,0,64);

    std::vector<uint64_t> fn = map2vec(forward_edges);
	std::vector<uint64_t> bn = map2vec(backward_edges);

	sizes[0] = forward_targets.size();
	sizes[1] = fn.size();
	sizes[2] = bn.size();

	written_bytes += sizes.serialize(out, child, "sizes");
	written_bytes += serialize_vector<uint64_t>(forward_targets,out,child,"forward_targets");
	written_bytes += serialize_vector<uint64_t>(backward_targets,out,child,"backward_targets");
	written_bytes += serialize_vector<uint64_t>(fn,out,child,"forward_edges");
	written_bytes += serialize_vector<uint64_t>(bn,out,child,"backward_edges");

	structure_tree::add_size(child, written_bytes);
	return written_bytes;
}

/**
 * @brief loads pos_EDS_to_FA_type with sdsl methods
 */
void adjacency::load(std::istream& in){
	sdsl::int_vector<64> sizes;
	sizes.load(in);
	
	forward_targets  = vector<uint64_t>(sizes[0]);
	backward_targets = vector<uint64_t>(sizes[0]);

	load_vector<uint64_t>(forward_targets,in);
	load_vector<uint64_t>(backward_targets,in);
	
	vector<uint64_t> tmp = vector<uint64_t>(sizes[1]);
	load_vector<uint64_t>(tmp,in);
	forward_edges = vec2map(tmp);

    tmp = vector<uint64_t>(sizes[2]);
	load_vector<uint64_t>(tmp,in);
	backward_edges = vec2map(tmp);
	
	initialised =true;
}

/** \brief for debugging */
std::ostream& operator<< (std::ostream& os, const adjacency& a){
	os << "adjacency-list:" << endl;
	os << "FORWARD:" << endl;
	for(auto i = a.forward_edges.cbegin(); i != a.forward_edges.cend(); i++ ){
		cout << i->first << " -> " << a(i->first) << endl;
	}	
	os << "BACKWARD:" << endl;
	for(auto i = a.backward_edges.cbegin(); i != a.backward_edges.cend(); i++ ){
		cout << i->first << " -> " << a(i->first,adjacency::BACKWARD) << endl;
	}
	return os;	
}
