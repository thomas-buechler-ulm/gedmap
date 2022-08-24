#pragma once
/**
 * This File includes code from the Master Thesis of Thomas Buechler.
 * 
 * Later published as "An improved encoding of genetic variation in a Burrows - Wheeler transform" 
 * in  "Bioinformatics, Volume 36, Issue 5, March 2020, Pages 1413 - 1419"
 * 
 * 
 * Purpuse is to read variations from a VCF  file and combine variations if  needed.
 * 
 */
using namespace std;



/************************************************************************/
/************************      VCF - Line        ************************/
/************************************************************************/
/**
 * (In Master Thesis: 2.3.3 VCF)
 * \brief this struct contains the data of a vcf line
 */
struct vcfl{
	//data
	std::string	chrom;
	uint64_t	pos;
	std::string	id;
	std::string	ref;
	std::string	alt;
	std::string	info;
	//uint32_t line; //line in the vcf file

	//constructors:
	
	/** \brief parameterized constructor*/ 
	vcfl(const std::string &chrom, const int pos, const std::string &id, const std::string & ref, const std::string & alt, const std::string & info);
	/** \brief input parsing constructor*/
	vcfl(const string& line);
	/** \brief default constructor*/   
	vcfl();
	
	std::deque<std::string> genaltQ() const;
	/** \brief true if the line encodes a SNP */
	bool isSNP() const;
	
	/** \brief checks if variant is an Alt (if alternative is specified)*/
	bool isALT() const;
	
	/** \brief checks if variant is an INDEL*/
	bool isINDEL() const;
	
	/** \brief checks if variant is a Copy Number variation*/
	bool isCNV() const;
	
	/** \brief checks wheater a vcf line is equal to this line */
	bool equal(vcfl & buf) const;
	
	/** \brief checks wheater a vcf line is the empty struct */
	bool empty() const;

	/** \brief returns number of commas in alt string +1 */
	int altCount() const;
	
	/** \brief returns max copy number*/
	uint32_t getCN() const;
	
	/** \brief searches for ";END=x;" in info string and returns x */
	uint64_t sv_end() const;
	
	
	/** \brief VAL FOR COUNT, BEGIN, END */
	typedef  std::tuple<uint32_t,uint32_t,uint32_t> CN_INFO;
	

	const static uint32_t COPY_0 = 1;
	const static uint32_t COPY_2 = 2; // COPY >= 2;
};



std::ostream& operator <<(std::ostream& outputStream, const vcfl& v);

/** 
 * (In Master Thesis: 6.6 Konflikt bei ueberschneidungen) 
 *\brief struct that contains several vcf lines that are in conflict
 */





/** \brief calculates queue of alternatives */
deque<string> vcfl::genaltQ() const{
	deque<string> altQ;
	int alt_begin 	= 0;
	int alt_length	= 0;
	for(unsigned int i = 0; i < alt.size(); i++)
		if(alt[i] != ',')	
			alt_length++;
		else{    //if "," add substring to altQ
			altQ.push_back(alt.substr(alt_begin,alt_length));
			alt_begin 	= i+1;
			alt_length 	= 0;
		}
	altQ.push_back(alt.substr(alt_begin,alt_length));
	
	return altQ;
}


vcfl::vcfl(const string &chrom, const int pos, const string &id, const string & ref, const string & alt, const string & info)
:chrom(chrom),pos(pos),id(id),ref(ref),alt(alt),info(info){}

/** \brief parses VCF line entry read from file*/
vcfl::vcfl(const string & line){
	unsigned int i 		= 0;
	unsigned int chromend	= 0;
	unsigned int idbegin	= 0;
	unsigned int idend 	= 0;
	unsigned int posbegin	= 0;	
	unsigned int posend 	= 0;
	unsigned int refbegin	= 0;
	unsigned int refend 	= 0;
	unsigned int altbegin	= 0;
	unsigned int altend 	= 0;
	unsigned int infobegin	= 0;
	unsigned int infoend 	= 0;
	//parse the line
	while(!isspace(line[i])		&& i < line.length()) i++; 	//i = chromend
	chromend	= i;  
	
	while( isspace(line[i])		&& i < line.length()) i++;	//i = posbegin
	posbegin	= i;
	while(!isspace(line[i])		&& i < line.length()) i++; 	//i = posend  
	posend	= i;
	
	while( isspace(line[i])		&& i < line.length()) i++;	//i = idbegin
	idbegin	= i;
	while(!isspace(line[i]) 	&& i < line.length()) i++; 	//i = idend  
	idend = i;

	while( isspace(line[i]) 	&& i < line.length()) i++;	//i = refbegin
	refbegin	= i;
	while(!isspace(line[i]) 	&& i < line.length()) i++; 	//i = refend 
	refend	= i;  
	
	while( isspace(line[i]) 	&& i < line.length()) i++;	//i = altbegin
	altbegin	= i;
	while(!isspace(line[i]) 	&& i < line.length()) i++; 	//i = altend
	altend	= i;    
	
	while( isspace(line[i])		&& i < line.length()) i++;	//i = QUALbegin
	while(!isspace(line[i]) 	&& i < line.length()) i++; 	//i = QUALend  

	while( isspace(line[i])		&& i < line.length()) i++;	//i = FILTERbegin
	while(!isspace(line[i]) 	&& i < line.length()) i++; 	//i = FILTERend  

	while( isspace(line[i])		&& i < line.length()) i++;	//i = INFObegin
	infobegin = i;
	while(!isspace(line[i]) 	&& i < line.length()) i++; 	//i = INFOend 
	infoend = i;
	//	initalise data
	chrom = line.substr(0,chromend);
	id	= line.substr(idbegin,idend-idbegin);
	ref	= line.substr(refbegin,refend-refbegin);
	alt 	= line.substr(altbegin,altend-altbegin);  
	info	= line.substr(infobegin,infoend-infobegin); 
	pos	= stoi(line.substr(posbegin,posend-posbegin)); 
}


vcfl::vcfl(){}


/**
 * \brief checks if a vcf line is a SNP 
 */
bool vcfl::isSNP() const{
	//ref must be one base
	if (ref.size()>1)
		return false;

	int altlength = 0;
	for(unsigned int i = 0; i < alt.size(); i++){
		altlength++;

		// with ',' the next alternative starts
		if(alt[i] == ',')
			altlength = 0;

		//length of each alternative must be one
		if (altlength>1)
			return false;    
    }
    return true;  
}

/**
* \brief checks if alternative is specified
*/
bool vcfl::isALT() const{
	return alt[0]!='<';
}

bool vcfl::isINDEL() const{
	if(alt.find(',') !=  string::npos||ref.size()==alt.size())
		return false;

	if(ref.size()==1 && alt[0]==ref[0])
		return true;
	
	if(alt.size()==1 && alt[0]==ref[0])
		return true;
	
	return false;
}

bool vcfl::isCNV() const{
	return (alt.find("<CN") !=  string::npos);
}

	
	
	/**
	 * \brief checks wheater a vcf line is equal to this line
	 */
	bool vcfl::equal(vcfl & buf) const{
		return(pos == buf.pos && ref == buf.ref && alt == buf.alt && chrom == buf.chrom);
	}
	
	bool vcfl::empty() const{
		return (pos == 0 && chrom == "");		
	}

	/** \brief returns number of commas in alt string +1 */
  int vcfl::altCount() const{
		unsigned int out = 1;    
			for(unsigned int i = 0; i < alt.size(); i++)
				if(alt[i] == ',') out++;
		return out;
	}
	
	/** \brief returns if CN0 and/or CNX X>1 exist */
  uint32_t vcfl::getCN() const{
		//alt has the form <CNx1>,<CNx2>,... with integers x1,x2..
		//reads numbers between 'N' and '>' and 
		size_t x_pos2 = 0;
		uint32_t ret_val = 0; 
		int ac = altCount();
		for(int i = 0; i < ac; i++){
			size_t x_pos1 = alt.find('N',x_pos2)+1;
			x_pos2 = alt.find('>',x_pos1);
			int x = stoi(alt.substr(x_pos1,x_pos2-x_pos1));
			if(x == 0) ret_val |= COPY_0;
			if(x >= 2) ret_val |= COPY_2;
		}	
		return ret_val;
	}
	
	
/**
 * \brief searches for ";END=x;" in info string and returns x
 */
uint64_t vcfl::sv_end() const{
	size_t end_pos1 = info.find(";END=")+5;
	size_t end_pos2 = info.find(';',end_pos1);
	return stol(info.substr(end_pos1,end_pos2-end_pos1));
}

/** \brief defiening an string for a variation for streams, for debugging*/
ostream& operator <<(ostream& outputStream, const vcfl& v){    
	if(v.isINDEL())
		outputStream << v.chrom << '\t' << v.pos << '\t' << v.id << '\t' << v.ref << '\t' << v.alt << "\t*\t*\t*";
	else
		outputStream << v.chrom << '\t' << v.pos << '\t' << v.id << '\t' << v.ref << '\t' << v.alt << "\t*\t*\t" << v.info;
	return outputStream;
}

/*
 * 
 * 
 * 
 * 
 *		CONFLICT GROUP 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 */
struct conflict_group{
	//data
	/** lines that are in conflict */
	std::deque<vcfl>	vcf_lines;
	/** last position  that is in conflict*/
	uint32_t maxend;
	/** \brief defalut constructor */
	conflict_group(){
		vcf_lines = std::deque<vcfl>();
		maxend 	= 0; 
	};
	
	//methods  
	/** \brief check if line is in conflict with the group */
	bool conflict(const vcfl & line) const;  
	/** \brief solves conflict, remains clean struct*/
	vcfl solve();  
	/** \brief adds a line to the group */
	void add(const vcfl & line);   
	/** \brief clears object */
	void clear();   
	/** \brief if group is empty*/
	bool empty() const;   
	/** \brief generates the reference string of the conflict group */
	std::string generate_ref() const;  
	/** \brief generates alternatives and adds them to a list */
	void combine(const std::string & ref, const uint32_t & p, std::string & new_alts, uint32_t i = 0,uint32_t off = 0,  std::string tmpalt = "") const;  
	/** \brief applies all SNPS in the list to ref */
	void apply_snps(std::string & ref);  
	/** \brief returns true if all lines are SNPs */
	bool only_snps() const;  
};

std::ostream& operator <<(std::ostream& outputStream, const conflict_group& c);




/** \brief check if (the next) line is in conflict with the group */
bool conflict_group::conflict(const vcfl & line) const{  
	if(empty())
		return false;  
	if(vcf_lines.front().chrom != line.chrom)
		return false; 
	return   maxend >= line.pos; 
};

/** \brief adds a line to the group */
void conflict_group::add(const vcfl & line){  
	unsigned int end 	= line.pos + line.ref.size() -1;	//update maxend
	if(end > maxend)
		maxend 	= end;  
	vcf_lines		.push_back(move(line)); 		//add line  
};

/** \brief if group is empty*/
bool conflict_group::empty() const{
	return vcf_lines.empty();
};

/**
 * \brief generates the reference string of the conflict group 
 * 
 * reads the ref of the first vcf_line in the queue
 * if chars for whole reference are missing, check next vcf_line
 * 
 * ends if all chars of the reference are read out of the vcf_lines
 */
string conflict_group::generate_ref() const{

	unsigned int startpos	= vcf_lines.front().pos;		//reference starts at the pos of the first ref in the queue
	unsigned int refsize 	= maxend - startpos +1;			//reference ends at maxend


	char ref [refsize];						//buffer for refstring
	int j 	= 0;						// current vcf line

	for(unsigned int i = 0; i <  refsize; i++)  			//loop over chars of reference
		if(vcf_lines[j].pos + vcf_lines[j].ref.size() -1 >= i + startpos) 	// what is possible >= what is needed
			ref[i] 	= vcf_lines[j].ref[i + startpos - vcf_lines[j].pos];	// define char
		else{ // info lies in next vcf line
			j++; // skip to next line
			i--; // redo loop for this char
		}
	return move(string(ref,refsize));
};


/**
 * \brief generates list of all possible alternative combinations
 * 
 * uses rekursion over vcf lines,
 * 
 * base case is when no vcf line is left.
 * 
 * in each step 2 rekursions can be made
 * 
 * the first is ignoring the line
 * the second is applying the line (only called when possible)
 * 
 * 
 * \param ref		reference string
 * \param p		reference-string-position 	(vcf_lines[0].pos)
 * \param alts		queue where result is stored
 * \param i		line counter 			(rekursion control, default = 0)
 * \param rso		reference-string-offset 	(default = 0)
 * \param tmpalt	temporary alt 			(default = "")
 */
void conflict_group::combine(const string& ref, const unsigned int& p, string& new_alts, unsigned int i, unsigned int off, string tmpalt) const{

	if(i == vcf_lines.size()){ 		//base case, no more vcf lines
		tmpalt	+= ref.substr(off);	//copy end of ref to alt
		new_alts	+= "," +(tmpalt);	//add alt to list
		return;				
	}
	//case 1: ignore vcf_line i
	combine(ref, p, new_alts, i+1, off, tmpalt);			// go to rekursion case 1 
	//case 2: if vcf_lines[i] is possible apply it
	unsigned int off_l = vcf_lines[i].pos - p;
	if(off <= off_l){ 					// if line i can be applied, apply it
		tmpalt 	+= ref.substr( off, off_l - off);	// copy refstring till begin of ref of this line
		// update reference-string-offset 
		off 	= off_l + vcf_lines[i].ref.size();				// ignore ref of line i
		//iterate over alts of this vcfl
		deque<string> altQ = vcf_lines[i].genaltQ();
		for(auto alt = altQ.begin(); alt != altQ.end(); alt++)
		     combine(ref, p, new_alts, i+1, off, tmpalt + *alt);//add alternative and go to rekursion case 2
	}
};

/** \brief this method fueses the alternatives of the snps */
vcfl fuse_snp(deque<vcfl> & vcf_lines){
	auto vl = vcf_lines.begin();
	string 	alt	= vl->alt;
	vl++;
	while(vl != vcf_lines.end()){
		alt += "," + vl->alt;
		vl++;
	}
	vcf_lines[0].alt = alt;
	return vcf_lines[0];  
};

/** \brief solves conflict, cleans struct this */
vcfl conflict_group::solve(){  
	switch(vcf_lines.size() ){
		case(0):	throw invalid_argument("ERROR: TRY TO SOLVE EMPTY CONFLICT_GROUP");	//empty group cannot be solved			
		case(1):	return(vcf_lines.front());					//just use front
		default:	

			if(only_snps())	
				return(move(fuse_snp(vcf_lines)));
				
				//generate alt
				//generate altQs
				//for(auto vl = vcf_lines.begin(); vl != vcf_lines.end(); vl++){
				//	vl->genaltQ();
				//}
				
			string 		ref	= generate_ref();
// 			apply_snps(ref);
			string 	new_alts	= "";
			combine(ref ,vcf_lines[0].pos,new_alts);
			new_alts = new_alts.substr(new_alts.find (',',1)+1);//ignore first, this is ref  
			vcf_lines[0].id += "_extended";
			vcf_lines[0].ref = ref;
			vcf_lines[0].alt = new_alts;
			//return
			return vcf_lines[0]; //make vcfl struct
	}
};

/** \brief solves conflict, cleans struct this */
void conflict_group::clear(){
	vcf_lines	.clear();
	maxend 	= 0;  
};

/** \brief extracts vcf_lines that are SNPS and integrates the SNPs in ref.  
 * changes ref and vcf_lines 
 */
// void conflict_group::apply_snps(string & ref){
// 	if(empty()) return;
// 
// 	for( deque<vcfl>::iterator vl = vcf_lines.begin() ; vl != vcf_lines.end(); vl++)
// 		if(vl->isSNP()){
// 			ref[vl->pos - vcf_lines.front().pos] = getSNPLetter(vl->alt,vl->ref[0]);
// 			vl = vcf_lines.erase(vl);
// 			vl--;
// 		}
// }

/** \brief returns true if all lines are SNPs */
bool conflict_group::only_snps() const{  
	for( auto vl = vcf_lines.begin() ; vl != vcf_lines.end(); vl++)
		if(!vl->isSNP()) return false;
	return true;
};

/** \brief defiening an string for a variation for streams, for debugging*/
ostream& operator <<(ostream& outputStream, const conflict_group& c){
	outputStream << "/*COGRU:" << endl;
	for(auto l = c.vcf_lines.begin(); l != c.vcf_lines.end(); l++)
		outputStream << "-" << *l;
	outputStream << "mxe: " << c.maxend << endl;
	outputStream << "ref: " << c.generate_ref() << "*/" << endl;
	return outputStream;
};
