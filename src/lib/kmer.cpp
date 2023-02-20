#pragma once
/**
 * @brief 64 bit number that represents a string with alphabet {A,C,G,T}
 * 
 * this will not be checked here, but is essential:
 * 
 * length <  2* (bit width of interger_type)
 */
template<typename interger_type>
struct KMER{
	static constexpr uint8_t A_CODE = 0;
	static constexpr uint8_t C_CODE = 1;
	static constexpr uint8_t G_CODE = 2;
	static constexpr uint8_t T_CODE = 3;
	
	interger_type content;
	uint32_t length;
	
	/** @brief default constructor */
	KMER(){
		content = 0;
		length = 0;
	}
	
	/** @brief parameterized constructor */
	KMER(interger_type content, uint32_t length):content(content),length(length){};	
	
	/** @brief copy constructor */
	KMER(const KMER & copy):content(copy.content),length(copy.length){};
	
	/** @brief constructor from string*/
	KMER(const std::string_view & kmer){
		content 	= 0;
		length	= kmer.size();
		for (const auto& c : kmer) {
			content <<= 2;
			switch (c) {
				case 'C': content += C_CODE; break;
				case 'G': content += G_CODE; break;
				case 'T': content += T_CODE; break;
				// TODO: what if 'N'?
			}
		}
	};
	
	/** @brief return string representation*/
	std::string get_kmer(){		
		std::string kmer = std::string(length,'A');
		interger_type cont = content;		
		for(uint32_t i = length; i > 0; i--){		
			switch(cont % 4){
				case C_CODE: kmer[i-1] = 'C';break;
				case G_CODE: kmer[i-1] = 'G';break;
				case T_CODE: kmer[i-1] = 'T';break;
				case A_CODE: break;
			}
			cont >>= 2;
		}		
		return kmer;
	};
	
	/** @brief adds a letter at front*/
	void add_front(char c){
		interger_type code;
		switch(c){
			case 'A': code = A_CODE; break;
			case 'C': code = C_CODE; break;
			case 'G': code = G_CODE; break;
			case 'T': code = T_CODE; break;
			default: [[unlikely]] throw std::invalid_argument("ADD_FRONT: illegal letter ASCII=(" + std::to_string((uint32_t) c) + ")");
		}		
		code <<= (2*length);
		content += code;
		length++;		
	}
	
	/** @brief adds a letter at back*/
	void add_back(char c){
		content = KMER::add_back(content, c);
		length++;		
	}
	
	/** @brief adds a letter at back
	 * add_back_force does not throw error if c is no base letter, c will instead be converted to 'A'
	 */
	void add_back_f(char c){
		if( c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'A';
		content = KMER::add_back(content, c);
		length++;		
	}
	
	/** @brief drops the front letter*/
	void rm_front(){
		assert(length > 0);
		//interger_type front = content;
		//front >>= (2*length - 2);
		//front <<= (2*length - 2);
		//content -= front;
		length--;
		content &= (static_cast<interger_type>(1) << (2 * length)) - 1;
	}
	
	
	static interger_type add_back(interger_type content, char c){
		interger_type code;
		switch(c){
			case 'A': code = A_CODE; break;
			case 'C': code = C_CODE; break;
			case 'G': code = G_CODE; break;
			case 'T': code = T_CODE; break;
			default: [[unlikely]] throw std::invalid_argument("KMER ADD_BACK: illegal letter ASCII=(" + std::to_string((uint32_t) c) + ")");
		}		
		content <<= 2;
		content += code;
		return content;
	}
	

	static interger_type number_of_kmers(uint32_t k){
		interger_type kmer_count = 1;
		for(uint32_t i = 0; i < k; i++)
			kmer_count <<= 2;
		return kmer_count;
	}
	
	static interger_type max_id(uint32_t k){
		return  number_of_kmers(k)-1;
	}
	
	static interger_type to_id(std::string_view pattern){
		KMER kmer(pattern);
		return kmer.content;
	}	
};
