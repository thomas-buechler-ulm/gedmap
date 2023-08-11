#pragma once
/**
 * \author Thomas Buechler
 * @brief defintion aligingment column and alignment table
 */


#include<vector>
#include<string>
namespace edsm_levinstein{
/**
 * @brief a column of the alignment table
 */
struct ali_column{
	sdsl::int_vector<8>  D; //contains number of errors still allowed + 1 (0 == alignment failed)
	uint32_t j; //corresponding postion in pattern
	uint32_t first_non_zero; // smallest i with D[i] != 0
	uint32_t last_non_zero; // greatest i with D[i] != 0
	
	char c;
	bool resolved;
	
	std::vector<uint32_t> predecessors; //for columns that were calculated by resolve_poss the poss vector;
	
	/** @brief empty constructor*/
	ali_column(){};
	
	/** @brief copy constructor */ //USE IMPICIT CONDTRUCTOR
// 	ali_column(const ali_column & original); 
	
	/** @brief constructor that creates vector D*/
	ali_column(uint32_t l, uint32_t j, char c);
	
	/** @brief sets first_non_zero and last_non_zero*/
	void set_first_last_zero();

	void swap(ali_column & x);
	
	/** @brief true if no value >0 contained*/
	bool empty(){return D.size() == 0 || first_non_zero == D.size();}
	
	uint32_t get_start(bool direction = FORWARD){
		if(direction==BACKWARD) return j - last_non_zero;
		else return j + first_non_zero;
	}
	
	uint32_t get_end(bool direction = FORWARD){
		if(direction==BACKWARD) return j - first_non_zero;
		else 	return j + last_non_zero;
	}
	
	uint32_t get_size(){
		if(empty())	return 0;
		return last_non_zero - first_non_zero + 1;
	}
	
	sdsl::int_vector<8>::reference  operator[] (uint32_t n){return D[n];}
	
	uint32_t max(){
		uint32_t max = 0;
		for(uint32_t i = 0; i < D.size(); i++)
			if( ((uint32_t) D[i]) > max )
				max = D[i];
		return max;
	}
};


/**
 * @brief the alignment table as vector of columns
 */
struct ali_table{
	std::vector<ali_column> table;		// alignment table	
	std::vector<uint32_t>	poss_index;		// columns that are possabilitys (will be unioned at next ')' or '(' )
	uint32_t		enter_index;	// column where alternative part was entered

	/** @brief initialisise table with first collumn*/
	ali_table(ali_column init);
	
	/** @brief unions all possible columns (specified by poss_index vector) and adds new col to A*/
	void union_poss(char c, bool direction);
	
	
	bool ended(){ return table.back().empty() && (enter_index == (uint32_t) -1 ) && poss_index .empty();}
	
	/** @brief copies column A[id] and add it to A*/
	void add_column_copy(uint32_t id, char c, vector<uint32_t> pred = vector<uint32_t>());
	
	/** @brief adds col it to A*/
	void add_col(const  ali_column & ac){table.push_back(ac);} 
	
	/**
	 * @brief trace back alignment in table
	 * traces back table[col][row] to table[0][i]. 
	 * Genereates cigar cig string (encoding of trace back)
	 * @return pair of cig and i.
	 */
	pair<string,uint32_t>  trace_back( uint32_t col, uint32_t row, bool direction);
};

// DEBUG PRINT
/** @brief prints vector for debugging*/
std::ostream& operator <<(std::ostream& os, const std::vector<uint32_t> & vec){
	os << '<';
	for (auto it = vec.begin() ; it != vec.end(); ++it)
		os << (it==vec.begin()?"":",") << *it;
	os << '>';
	return os;
};
std::ostream& operator <<(std::ostream& os, const ali_column& ac){	
	char c = ac.c;
	if(ac.c == 0)
		c = '0';
	os << c << "("<< ac.j<<"):";
	
	for(uint32_t i = 0; i < ac.D.size();i++){
		os << (uint32_t) ac.D[i] << " ";
	}
	
	os << " (" << ac.first_non_zero << "-" << ac.last_non_zero << ")";
	os << "  p:" << ac.predecessors;
		
	return os;
}



std::ostream& operator <<(std::ostream& os, const ali_table& at){
	using namespace std;
	std::cout  << "*************TABLE***************" << endl;
	for(uint32_t i = 0; i < at.table.size(); i++){
		os << setfill(' ') << std::setw(4) << i << ":" << setfill(' ') << std::setw(10) << at.table[i].j << ": ";
		os << at.table[i];
		os << endl;
	}
	std::cout << "enter: " << at.enter_index << endl;
	std::cout << "pos  : " << at.poss_index << endl;
	std::cout << "*********************************" << endl;
	return os;
}

ali_column::ali_column(uint32_t l, uint32_t j, char c  = 0):j(j),c(c){
	D = int_vector<8>(l,0,0);
	first_non_zero = 0;
	last_non_zero = l-1;
	resolved = false;
}

// ali_column::ali_column(const ali_column & original){
// 	D = int_vector<8>(original.D.size(),0,0);
// 	for(uint32_t i = 0; i < D.size(); i++ )
// 		D[i] = original.D[i];
// 	j = original.j;
// 	c = original.c;
// 	first_non_zero = original.first_non_zero;
// 	last_non_zero = original.last_non_zero;
// 	resolved = original.resolved;
// 	predecessors = vector<uint32_t>(original.predecessors);
// }

void ali_column::set_first_last_zero(){
	if(D.size() == 0) return;  //trivial
	
	
	first_non_zero = D.size();
	last_non_zero = 0;			
	
	for(uint32_t i = 0; i < D.size(); i++){//find first non zero
		if(D[i]){
			first_non_zero = i;
			break;
		}
	}
	
	if(first_non_zero == D.size()){//Only Zeros
		last_non_zero = 0;
		return;
	}
	
	for(uint32_t i = D.size(); i > 0; i--){//find last non zero
		if(D[i-1]){
			last_non_zero = i-1;
			break;
		}
	}
}

void ali_column::swap(ali_column & x){
	D.swap(x.D);		
	std::swap(j,x.j);	
	std::swap(first_non_zero,x.first_non_zero);	
	std::swap(last_non_zero,x.last_non_zero);
}

ali_table::ali_table(ali_column init){
	table 	= vector<ali_column>(1,init);
	poss_index 	= vector<uint32_t>(0);
	enter_index = -1;
}

void ali_table::union_poss(char c, bool direction){
	if(poss_index.empty()){
		table.push_back(ali_column());
		return;
	}		

	int start = 2147483647; //INT MAX
	uint32_t end = 0;
	
	for(uint32_t pred : poss_index){
		int i_start		= table[pred].get_start(direction);
		uint32_t i_end 	= table[pred].get_end(direction);
		
		if(i_start < start)
			start = i_start;
		if(i_end > end)
			end = i_end;			
	}
	
	
	uint32_t new_j = (direction==BACKWARD)? end : start;

			
	ali_column new_col = ali_column(end-start + 1,new_j,c);
	
	for(uint32_t pred : poss_index){
		int off;
		if(direction==BACKWARD)
			off = new_j - table[pred].j;
		else
			off = table[pred].j - new_j;
		
				
		for(uint32_t k = table[pred].first_non_zero; k <= table[pred].last_non_zero; k++)
			if( new_col[ k+off ] < table[pred][k])
				new_col[ k+off ] = table[pred][k];
		
		
	}
	
	new_col.set_first_last_zero();
	new_col.resolved = true;
 	new_col.predecessors.swap(poss_index);
	table.push_back(new_col);
}

void ali_table::add_column_copy(uint32_t id, char c, vector<uint32_t> pred){
	ali_column ac = ali_column(table[id]);
	ac.predecessors = pred;
	ac.c = c;
	ac.resolved = true;
	table.push_back (ac);
}

std::pair<std::string,std::uint32_t> ali_table::trace_back( uint32_t col, uint32_t row, bool direction ){
	using namespace edsm_levinstein; //FORWARD/BACKWARD
	uint32_t D_val 		= table[col][row];
	uint32_t max_Length 	= table.size() + table[0].max() - D_val;
	std::string cig 			= std::string(max_Length,'?');
	uint32_t j;
	int one;
	
	if(direction == edsm_levinstein::FORWARD){
		j = max_Length - 1;
		one = 1;
	}else{
		j = 0;
		one = -1;
	}
	
	while( col != (uint32_t) 0){


		if(table[col].resolved == false) { //ALL LETTERS			
		
			uint32_t pred = col-1;
			

			
			if(table[col].predecessors.size() == 1)
				pred = table[col].predecessors[0];
			
			uint32_t off = 0;
			if(pred != (uint32_t) -1) 
				off = table[pred].first_non_zero;
			
/*
			//INSETION BASE CASE
			if( pred == (uint32_t) -1){ //TODO DELETE THIS CASE
				cig[j] = 'I';
				j =  j - one;
				row--,
				D_val++; 
			}*/

			//DELETION BASE CASE
			if(row==0 && row < table[col].D.size() - 1 
				&& D_val+1== (uint32_t)  table[pred][row + off]){
				cig[j] = 'D';
				j =  j - one;
				col = pred;
				row += off;
				D_val++; 
			}
			//INSERTION
			else if(row > 0 && D_val+1 == table[col][row-1]){  
				cig[j] = 'I';
				j =  j - one;
				row--,
				D_val++; 
			}
			//DELETION
			else if(row < table[col].D.size() - 1 
				&& D_val+1== (uint32_t)  table[pred][row + off]){
				cig[j] = 'D';
				j =  j - one;
				col = pred;
				row += off;
				D_val++; 
			}
			//MATCH (we do not check if letters are equal. But if the value cant be get by insertion or deletion it has to be a (mis)match)
			else if( D_val == (uint32_t) table[pred][row+off-1] ){
				cig[j] = 'M';
				col = pred;
				row = row + off -1;
				j =  j - one;
			}
			//MISMATCH
			else if ( D_val +1  == (uint32_t)  table[pred][row+off-1] ){ 
				cig[j] = 'M';					
				j =  j - one;
				col = pred;
				row = row + off -1; 
				D_val++; 
			}
			else				
				throw runtime_error ("trace_back failed");	
		}else{
		
			bool ok = false;
			
			for(uint32_t pred : table[col].predecessors){					
				int pred_row = row;
				if(direction==FORWARD)
					pred_row += table[col].j - table[pred].j;
				else						
					pred_row -= table[col].j - table[pred].j;
				
				if (pred_row < 0 || pred_row >= (int)(table[pred].D.size()))
					continue;					
				if( D_val == (uint32_t) table[pred][pred_row]){
					col = pred;
					row = pred_row;
					ok = true;
					break;
				}
			}
			if(!ok) throw runtime_error ("trace_back failed2");
		}
	}
	
	
// 	if (!backw && j != (uint32_t) -1 )
// 		return cig.substr(j+1);			
// 	if( backw && j !=  max_Length)
// 		return cig.substr(0,j);
		
	
	if (direction==FORWARD && j != (uint32_t) -1 ){
		cig = cig.substr(j+1);
	}
	
	if(direction==BACKWARD && j !=  max_Length){
		cig = cig.substr(0,j);
	}

	return make_pair(cig,row);
}
}
