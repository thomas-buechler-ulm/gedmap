#pragma once
namespace gedmap_encode{

/** \brief returns the complement to a base */
char complement_char(char & c){
	switch(c){
		case 'A': return 'T';
		case 'T': return 'A';
		case 'G': return 'C';
		case 'C': return 'G';
		case 'N': return 'N';
		case 'B': return 'V';
		case 'D': return 'H';
		case 'H': return 'D';
		case 'K': return 'M';
		case 'M': return 'K';
		case 'R': return 'Y';
		case 'S': return 'S';
		case 'V': return 'B';
		case 'W': return 'W';
		case 'Y': return 'R';
	}
	throw std::invalid_argument("INVALID_UPAC_LETTER in complement_char(): >" + std::string(1,c) + "<");
}

/** \brief calculates the 5'-3' reverse of a DNA pattern
 *  if complement = false, this just calculates the reverse pattern
 */
std::string rev_complement(std::string & pattern, bool complement = true){
  std::string inv = std::string(pattern);
  for(unsigned int i = 0; i < pattern.length() ; i++)
    inv[i] = complement?complement_char(pattern[pattern.length() - i - 1]):pattern[pattern.length() - i - 1];
  return inv;
}

std::string RL_encode(const std::string & s){
	std::stringstream  out;	
	std::string::const_iterator it = s.begin();	
	while(it != s.end()){
		char c = *it;
		unsigned int run_length = 1;
		it++;		
		while(it != s.end() && *it == c){
			it++;
			run_length++;
		}
		out << run_length << c;		
	}
	return out.str();
}

std::string RL_decode(const std::string & s){
	std::stringstream  out;	
	std::istringstream in = std::istringstream(s); 	
	while(true){					
		unsigned int run_length;			
		in >> run_length;		
		if(in.eof())
			return out.str();		
		if(in.fail())
			throw std::invalid_argument("invalid RL-encoding");		
		char c;		
		in.get(c);		
		for(;run_length > 0; run_length--)
			out << c;	
	}
}

}
