#pragma once
#include <fstream>
#include <string>
#include <memory>
#include "vcf_line.cpp"

struct fa_stream{
	std::ifstream ifs;
	std::string seq_name;
	size_t pos;

	void open(std::string fname){
		ifs.open(fname ,ios_base::in);
		getline(ifs, seq_name);
		seq_name = seq_name.substr(1,seq_name.find('|'));
		pos = 0;
		gedmap_io::print_row("READ SEQUENCE: " + seq_name);
	}

	void close(){
		ifs.close();
	}
	void skip(size_t x){
		//ifs.seekg(c,ios_base::cur); WRONG because of linebreaks
		//pos += c;
		size_t target_pos = pos + x;
		char c;
		while(pos < target_pos){
			if(!ifs.get(c))	throw std::invalid_argument( "FA STREAM ENDED BEVOR POSITION " + std::to_string(target_pos) +" OF VARIATION END WAS REACHED");
			if(isspace(c))		continue;
			if(c == '>') 		throw std::invalid_argument( "FA: POSITION " + std::to_string(target_pos) +" OF VARIATION END NOT FOUND IN SEQ "  + seq_name);
			if(!(c == 'A' || c == 'T' || c == 'G'|| c == 'C'|| c == 'N')){
				gedmap_io::print_error("UNEXPECTED SYMBOL: ASCII=(" + std::to_string((uint32_t) c) + ")='" + c + "' IN FA" + " at pos " + std::to_string(ifs.tellg()) + " replaced by ref of vcf");
			}
			pos++;
		}
	}

	void copy_to_ofs_until_t(size_t target_pos, std::ofstream& of){
		char c;
		while(pos < target_pos){
			if(!ifs.get(c)) throw std::invalid_argument( "FA STREAM ENDED BEVOR POSITION " + std::to_string(target_pos) +" OF VARIATION WAS REACHED");
			if(isspace(c)) continue;
			if(c == '>') throw std::invalid_argument( "FA: POSITION " + std::to_string(target_pos) +" OF VARIATION NOT FOUND IN SEQ "  + seq_name);
			if(!(c == 'A' || c == 'T' || c == 'G'|| c == 'C'|| c == 'N')){
				gedmap_io::print_error("UNEXPECTED SYMBOL: ASCII=(" + std::to_string((uint32_t) c) + ")='" + c + "' IN FA" + " at pos " + std::to_string(ifs.tellg()) + " replaced by 'N'");
				c = 'N';
			}
			of << c;
			pos++;
		}
	}

	void copy_to_ofs_until_next_seq(ofstream& ofs){
		char c;
		while(true){
			if(!ifs.get(c)) throw std::invalid_argument( "FA STREAM ENDED BEVOR SEQUENCE NEXT SEQUENCE WAS REACHED");
			if(c == '\n') continue;
			if(c == '>'){
				ofs << EDS_NODE_BOUNDARY;
				getline(ifs, seq_name);
				seq_name.substr(0,seq_name.find('|')); //cut additional infos
				pos = 0;
				gedmap_io::print_row("READ SEQUENCE: " + seq_name);
				return;
			}
			if ( !(c == 'A' || c == 'T' || c == 'G'|| c == 'C'|| c == 'N') ){
				gedmap_io::print_error("UNEXPECTED SYMBOL: ASCII=(" + std::to_string((uint32_t) c) + ")='" + c + "' IN FA" + " at pos " + to_string(ifs.tellg()) + "replaced by 'N'");
				c = 'N';
			}
			ofs << c;
			pos++;
		}
	}
};

struct vcf_stream{
	std::ifstream ifs;
	uint32_t lines_read;


	void open(std::string fname){
		ifs.open(fname ,ios_base::in);
		lines_read = 0;
	}

	void close(){
		ifs.close();
	}

	bool get(vcf_line& vcfl){
		string line;
		while(true){
	 		if(! getline(ifs, line) ) return false;
			if(line.size() == 0 || line[0] == '#') continue; //skip comments
			vcfl = vcf_line(line);
			if(++lines_read % 200000 == 0) gedmap_io::flush_row("Read vcf lines", to_string(lines_read));
			return true;
		}
	}

};

struct stream_to_stream_copy{
	ifstream* in;
	ofstream* out;
	std::unique_ptr<char[]> buffer;
	uint32_t buffer_size;

	stream_to_stream_copy(ifstream* in, ofstream* out, uint32_t buffer_size):in(in),out(out),buffer_size(buffer_size){
		buffer = std::make_unique<char[]>(buffer_size);
	};

	/**
	 * @brief writes specified amount of data from one stream to another
	 * @param in input stream
	 * @param out output stream
	 * @param buffer buffer for data
	 * @param buffer_size
	 * @param count number of chars to copy
	 */
	void copy(uint64_t count){
		while(count > buffer_size){
			in ->read (buffer.get(),buffer_size);
			out->write(buffer.get(),buffer_size);
			count -= buffer_size;
		}
		in->read (buffer.get(),count);
		out->write(buffer.get(),count);
	};

	/**
	 * @brief writes all data left in one stream to another
	 * @param in input stream
	 * @param out output stream
	 * @param buffer buffer for data
	 * @param buffer_size
	 */
	void copy_all(){
		uint64_t pos =  in->tellg();
		in->seekg (0, in->end);
		uint64_t pos_end =  in->tellg();
		in->seekg(pos);
		copy(pos_end - pos);
	};
};

struct bv_to_bv_copy{
	sdsl::bit_vector* in; //input bitvector
	sdsl::int_vector_buffer<1>* out;//output bitvector stream
	uint64_t	next_pos;

	bv_to_bv_copy(sdsl::bit_vector*  in, sdsl::int_vector_buffer<1>* out):in(in),out(out){next_pos = 0;};

	void copy(uint64_t count){
		uint64_t last_pos = next_pos + count;
		while(next_pos < last_pos) out->push_back( (*in)[next_pos++] );
	};

	void copy_all(){
		copy( in->size() - next_pos );
	};
};

