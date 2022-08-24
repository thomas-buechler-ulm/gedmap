#include <algorithm>    // std::max //min
#include <utility>      // std::swap
#include <tuple>      // std::tie
#include "encodings.cpp"
#include "alignment_table.cpp"

using namespace std;
using namespace sdsl;




namespace edsm_levinstein{


// template<typename pos_type>
// tuple<string,uint32_t,uint32_t>	align_forward (const string & EDS, const string & P, const uint32_t D, pos_type i, const int j, const adjacency & adj); //TODO WHY INT j?
// template<typename pos_type>
// tuple<string,uint32_t,pos_type,uint32_t> 	align_backward	(const string & EDS, const string & P, const uint32_t D, pos_type i, const uint32_t j, const adjacency & adj);


template<typename pos_type>
tuple<string,uint32_t,pos_type> align(const string & EDS, const adjacency & adj, const string & P, uint32_t D, pos_type i, uint32_t j, bool direction);


/**
 * @brief Map EDS[i] to P[j] and aling in both directions with max D errors
 * 
 * @template pos_type uint32_t iff long enough to address EDS uint64_t otherwise
 * 
 * @param EDS	chrom
 * @param P	pattern
 * @param D	max edit dist
 * @param i	pos in EDS
 * @param j	pos in pattern
 * @param adj	adjacency list
 * 
 * @return triple with cigar string, distance, start position in EDS
 */
template<typename pos_type>
tuple<string,uint32_t,pos_type> levinstein(const string & EDS, const adjacency & adj, const string & P, uint32_t D, pos_type i, uint32_t j){
	typedef tuple<string,uint32_t,pos_type> 	alignment; //CIGAR, D, current i
	
	if( j >= P.size())	throw runtime_error("in levinstein algorithm: j >= size of pattern");
	if( i >= EDS.size())	throw runtime_error("in levinstein algorithm: i >= size of eds");
	
	alignment a = align(EDS,adj,P,D,i,j,FORWARD);
	if(get<1>(a) > D) return make_tuple("",0,EDS.size());//alignment not possible
	
	if(j == 0)	return make_tuple( gedmap_encode::RL_encode(get<0>(a)) , get<1>(a), i);
	
	alignment a2 = align(EDS,adj,P,D-get<1>(a),i,j,BACKWARD);
	if(get<1>(a2) + get<1>(a) > D) return make_tuple("",0,EDS.size());//alignment not possible
	
	return make_tuple(
		gedmap_encode::RL_encode(get<0>(a2) + get<0>(a)),
		get<1>(a2) + get<1>(a),
		get<2>(a2)
	);
}


template<typename pos_type>
tuple<string,uint32_t,uint32_t> 		align_forward (const string & EDS, const string & P, const adjacency & adj, pos_type i, ali_column & init_col);
template<typename pos_type>
tuple<string,uint32_t,pos_type,uint32_t>	align_backward(const string & EDS, const string & P, const adjacency & adj, pos_type i, ali_column & init_col);



/**
 * @brief calls align_forward() or align_backward() with initial column
 * @param EDS 	EDS
 * @param P		pattern
 * @param adj	adjacency list for '#'-Symbols
 * @param max_D	number of maximally allowed errors
 * @param i		position in EDS
 * @param j		position in pattern
 * 
 * @return <cig,D,pos> with:
 * 	cig = cigar string of alignment
 * 	D   = number of errors in alignment
 * 	pos = starting position in alignment
 */
template<typename pos_type>
tuple<string,uint32_t,pos_type> align(const string & EDS, const adjacency & adj, const string & P, uint32_t max_D, pos_type i, uint32_t j, bool direction){
	string	cig;
	uint32_t	D,row,max_e;
	pos_type	pos;
	
	if(direction == FORWARD){
		
		//def init_col = column befor alignment starts
		uint32_t	j_init	= j-1;
		uint32_t	init_size 	= min(max_D+1,(uint32_t)( P.size()-j_init));
		ali_column	init_col 	= ali_column(init_size,j_init);
		for(uint32_t k = 0; k < init_size; k++) init_col[k] = max_D-k+1;
		
		
		//call align_forward();
		tie(cig,row,max_e) = align_forward(EDS,P,adj,i,init_col);
 
		
		
		if(row > 0) cig = string(row,'I') + cig;
		pos = i;
	}else{
		//def init_col = column befor alignment starts
			
	
		uint32_t init_size = min( max_D+1, j+1 );	
		ali_column init_col = ali_column(init_size, j);
	
		for(uint32_t k = 0; k < init_size; k++) init_col[k] = max_D-k+1;
		
		
		tie(cig,row,pos,max_e) = align_backward(EDS,P,adj,i,init_col);
		if(row > 0) cig = cig + string(row,'I');
		
	
	}
	D = max_D - (max_e-1);
	return make_tuple(cig,D,pos);
	
}





/**
 * @brief calculates a new column of A in base case
 */
ali_column new_column( ali_column & last_col, char c, const string & P, bool direction, uint32_t min_value = 0){
	if(last_col.empty()) 
		return ali_column(0,0,c);
	
		
	uint32_t new_j = (direction==BACKWARD)?
		last_col.get_end(direction) :
		last_col.get_start();
	
		
	uint32_t new_size = (direction==BACKWARD)? 
		min(last_col.get_size() + 1, new_j + 1) :
		min(last_col.get_size() + 1, (uint32_t) P.size() - last_col.get_start());
	
	
	ali_column new_col = ali_column(new_size,new_j,c);
	
	uint32_t off = last_col.first_non_zero;//offseali_tablet between last ali and new ali  start;
	
	
	
	for( uint32_t k = 0; k < new_size; k++){
		
		uint32_t p_pos = (direction==BACKWARD)?
			new_col.j - k:
			new_col.j + k;
		
		if(k > 0){
			//(MIS)MATCH
			uint32_t penalty = (c == P[p_pos] || c == 'N' || P[p_pos] == 'N' )?0:1;
			if(new_col[k] + penalty < last_col[k-1+off])
				new_col[k] = last_col[k-1+off] - penalty;
			
			//INSERTION
			if(new_col[k] + 1 < new_col[k-1])
				new_col[k] = new_col[k-1]-1;
		}

		if(k < new_size - 1){
			//DELETION
			if(new_col[k] + 1 < last_col[k+off] )
				new_col[k] = last_col[k+off] - 1;
		}
		
		if(new_col[k] < min_value) new_col[k] = 0;
	}
	new_col.set_first_last_zero();
	return new_col;
}



/**
 * @brief aligns pattern to EDS in forward direction.
 * @param EDS 	EDS
 * @param P		pattern
 * @param adj	adjacency list for '#'-Symbols
 * @param i		position in EDS
 * @param j		position in pattern
 * @param init_col first column in alignment table
 * 
 * @return <cig,row,max_e> with:
 * 	cig = cigar string of alignment
 * 	row = row in init_col where traceback ended
 * 	max_e = maximum number of left errors + 1
 */
template<typename pos_type>
tuple<string,uint32_t,uint32_t> align_forward(const string & EDS, const string & P, const adjacency & adj, pos_type i, ali_column & init_col){
	ali_table A = ali_table(init_col);

	uint32_t max_e = 0; //maxium of left erros
	uint32_t col = 0;
	uint32_t row = 0;
	string   cig = "";
	
	// WORK ON A SMALL COPY OF EDS
	//###### 
	//size_t EDS_BUFFER_START = i;
	//string EDS_BUFFER = EDS.substr(EDS_BUFFER_START,STRING_BUFFER_SIZE);
	//###### 
	
	
	while(!A.ended()){
		ali_column* last_col = &(A.table.back());
		
		
		if(!last_col->empty() && last_col->get_end() != uint32_t (-1) && last_col->get_end() >= P.size() -1 ){ //END OF PATTERN COVERED
			assert(last_col->D.size() <= P.size());			
			uint32_t val = (*last_col)[last_col->last_non_zero];
			if( val >  max_e){	//set new min
				max_e = val;
				col = A.table.size()-1;
				row = last_col->last_non_zero;
				
				//eliminate worse values TODO BACK
// 				last_col->enforce_min_value(val);
// 				last_col->set_first_last_zero();
				if(last_col->first_non_zero == last_col->D.size()-1
					&& (A.enter_index == (uint32_t) -1 )
					&& A.poss_index .empty()){ //no better values possible
					break;
				}
// 				if(A.enter_index != (uint32_t) -1 )A.table[A.enter_index].enforce_min_value(val);
// 				for(uint32_t pred : A.poss_index){
// 					A.table[pred].enforce_min_value(val);
// 				}
			}	
		}	
		//######   
		char c=EDS[i++] ;
		/*
		if( i >=  EDS_BUFFER_START + EDS_BUFFER.size()){
			if (i >= EDS.size() ) break; //TODO WHY IS THIS FASTER THAN EVERYTHING ELSE
			EDS_BUFFER_START = i;
			EDS_BUFFER = EDS.substr(EDS_BUFFER_START,STRING_BUFFER_SIZE);
		}
		char c = EDS_BUFFER[i-EDS_BUFFER_START];i++;
		*/
		//######
		
		if(c == EDS_NODE_BOUNDARY){
			ali_column bound_col(*last_col);
			bound_col.c = '#';
			
			std::vector<uint64_t> targets = adj(i-1).to_vector();
			for(uint64_t pos : targets){
				string   rec_cig;
				uint32_t rec_row;
				uint32_t rec_max_e;
				try{
					tie(rec_cig,rec_row,rec_max_e) = align_forward(EDS,P,adj,pos+1,bound_col);
					if(rec_max_e > max_e){
						max_e = rec_max_e;
						col = A.table.size()-1;//Col befor #
						row = rec_row;
						cig = cig + rec_cig;
					}
				}catch(runtime_error & e){
					//ali not possible
				}
			}break;
		}
		
		if( c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' )
			A.add_col(new_column(*last_col,c,P,FORWARD,max_e));
		
		else if( c ==  '(') //Enter ALT	
			A.enter_index = A.table.size()-1;
				
		else if( c =='|'){ // NEXT OPTION
				if(A.enter_index == (u_int32_t) -1){ //IF ENTER POS UNKNOWN -> JUMP TO ALT END 
					while(c != ')'){
						
						//######   
						c=EDS[i++];
						/*
						if( i >=  EDS_BUFFER_START + EDS_BUFFER.size() ){
							EDS_BUFFER = EDS.substr(i,STRING_BUFFER_SIZE);
							EDS_BUFFER_START = i;
						}
						c = EDS_BUFFER[i-EDS_BUFFER_START];
						i++;*/
						//###### 
					}
					continue;
				}		
				if(!A.table.back().empty())
					A.poss_index.push_back(A.table.size()-1); //ADD COLUMN AS POSSIBLE
				A.add_column_copy(A.enter_index,'|',vector<uint32_t>(1,A.enter_index));
		}
			
		else if( c == ')'){
			if(!A.table.back().empty())
				A.poss_index.push_back(A.table.size()-1);
			A.union_poss(')',FORWARD);
			A.enter_index = (uint32_t) -1; //REMOVE ENTER VALUE
		}
		
		else throw runtime_error("af in levinstein: UNEXPECTED SYMBOL: ASCII=(" + to_string((uint32_t) c) + ") at EDS[" + to_string(i-1) +"]" );
			
		
	}		
	
	if(max_e == 0)	return make_tuple("",0,0); //alignment not possible


	string cig_tb;
	tie(cig_tb,row) = A.trace_back(col,row,FORWARD);
	return make_tuple(cig_tb+cig,row,max_e);
}


/**
 * @brief aligns pattern to EDS in backward direction.
 * @param EDS 	EDS
 * @param P		pattern
 * @param adj	adjacency list for '#'-Symbols
 * @param i		position in EDS
 * @param j		position in pattern
 * @param init_col first column in alignment table
 * 
 * @return <cig,row,pos,max_e> with:
 * 	cig = cigar string of alignment
 * 	row = row in init_col where traceback ended
 * 	pos = position in EDS where alignment starts
 * 	max_e = maximum number of left errors + 1
 */
template<typename pos_type>
tuple<string,uint32_t,pos_type,uint32_t>	align_backward (const string & EDS, const string & P,  const adjacency & adj, pos_type i, ali_column & init_col){
	ali_table A = ali_table(init_col);

	
	uint32_t max_e = 0;
	uint32_t col = 0;
	uint32_t row = 0;
	string   cig = "";
	pos_type start_i = 0;
	
	
	//###### 
	// WORK ON A SMALL COPY OF EDS
	//size_t EDS_BUFFER_START = (i<STRING_BUFFER_SIZE)? 0 : (i-STRING_BUFFER_SIZE);
	//string EDS_BUFFER = EDS.substr(EDS_BUFFER_START,STRING_BUFFER_SIZE);
	//###### 
	
	while(!A.ended()){		
		ali_column* last_col = &(A.table.back());
		
		
		if(last_col->get_start(BACKWARD) == 0){ //END OF PATTERN COVERED
			
			uint32_t val = (*last_col)[last_col->last_non_zero];
			
			//eliminate worse values
			for(uint32_t j = last_col->first_non_zero; j < last_col->last_non_zero; j++)
				if(last_col->D[j] <= val)
					last_col->D[j] = 0;
			
			if( val >  max_e){	//set new min				
				max_e = val;
				col = A.table.size()-1;
				row = last_col->last_non_zero;
				start_i = i;
				
				
				last_col->set_first_last_zero();
				if(last_col->first_non_zero == last_col->D.size()-1
					&& (A.enter_index == (uint32_t) -1 )
					&& A.poss_index .empty()){ //no better values possible
					break;
				}
			}	
		}	
 		if( i == 0) break;	
		
		
		//###### 
		char c = EDS[--i]; 
		/*
		--i;
		if(i < EDS_BUFFER_START){			
			EDS_BUFFER_START = (i<STRING_BUFFER_SIZE)? 0 : (i-STRING_BUFFER_SIZE+1);
			EDS_BUFFER = EDS.substr(EDS_BUFFER_START,STRING_BUFFER_SIZE);
		}
		char c = EDS_BUFFER[i-EDS_BUFFER_START];
		*/
		//###### 
		
		if(c == EDS_NODE_BOUNDARY){
			ali_column bound_col(*last_col);
			bound_col.c = '#';
			
			std::vector<uint64_t> targets = adj(i, adjacency::BACKWARD).to_vector();
			for(uint64_t pos : targets){
				string   rec_cig;
				uint32_t rec_row;
				uint32_t rec_max_e;
				uint32_t rec_start_i;
				try{
					tie(rec_cig,rec_row,rec_start_i,rec_max_e) = align_backward(EDS,P,adj,pos,bound_col);
					if(rec_max_e > max_e){
						max_e = rec_max_e;
						col = A.table.size()-1;//Col befor #
						row = rec_row;
						cig = rec_cig + cig;
						start_i = rec_start_i;
					}
				}catch(runtime_error & e){
					//ali not posssible
				}
			}break;
		}
		
		if( c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' ) 
			A.add_col(new_column(*last_col,c,P,BACKWARD,max_e));
			
				
		else if (c ==  '|'){ // NEXT OPTION
			if(A.enter_index == (uint32_t) -1){
				while(c != '('){
					
					//###### 
					c = EDS[--i];
					/*
					i--;
					if(i < EDS_BUFFER_START){			
						EDS_BUFFER_START = (i<STRING_BUFFER_SIZE)? 0 : (i-STRING_BUFFER_SIZE+1);
						EDS_BUFFER = EDS.substr(EDS_BUFFER_START,STRING_BUFFER_SIZE);
					}
					c = EDS_BUFFER[i-EDS_BUFFER_START];
// 					*/
					//###### 
					
					
				}
				continue;
			}	
			if(!A.table.back().empty())
				A.poss_index.push_back(A.table.size()-1);
			A.add_column_copy(A.enter_index,'|',vector<uint32_t>(1,A.enter_index));
		}
				
		else if (c ==  ')' ) //Enter ALT	
			A.enter_index = A.table.size()-1;
		
		else if (c == '('){
			
			if(A.enter_index == (uint32_t) -1)
				continue; //ignore syntax symbol
			
			if(!A.table.back().empty())
				A.poss_index.push_back(A.table.size()-1);
			
			A.union_poss('(',BACKWARD);
			A.enter_index = (uint32_t) -1;
		}
		
		else	throw runtime_error("ab in levinstein: UNEXPECTED SYMBOL: ASCII=(" + to_string((uint32_t) c) + ")= "+c+" at EDS[" + to_string(i) +"]" );
		
	}		
	
	if(max_e == 0)	return make_tuple("",0,0,0); //alignment not possible

	//rebuild CIGAR
	string cig_tb;
	tie(cig_tb,row) = A.trace_back(col,row,BACKWARD);	
	return make_tuple(cig+cig_tb,row,start_i,max_e);
}

}
