#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>      // std::setw

using namespace std;

void front(string & s){s = s.substr(0, s.find(' '));} //sets s to its first word

vector<string> string_split(string & s, char delim){
	vector<string> out(0);
	size_t l = 0;
	size_t r;
	while(l < s.size()){
		r = s.find(delim,l);
		if(r == string::npos){
			out.push_back(s.substr(l));
			return out;
		}else{
			out.push_back(s.substr(l,r-l));
			l = r+1;
		}
	}
	throw runtime_error ("Strange things in my string_split implementation\n");
}


int main(int argc, char *argv[]){

	if(argc != 4){	
		cout << "3 Argument  expected fastq-file  sam-file-1 sam-file-2" << endl; 
		return 0;
	}

	ifstream fq_in, sam1_in, sam2_in;
	fq_in.open(argv[1],ios_base::in);
	sam1_in.open(argv[2],ios_base::in);
	sam2_in.open(argv[3],ios_base::in);

	string last_read_id,read_id, ali1, ali2;
	
	//get alignment line
	getline(sam1_in, ali1);
	getline(sam2_in, ali2);
	
	//skip comments '@'
	while(ali1[0] == '@' && getline(sam1_in, ali1));
	while(ali2[0] == '@' && getline(sam2_in, ali2));
	
	vector<uint32_t> errors = vector<uint32_t>(40,0);

	uint32_t only_1 = 0;
	uint32_t only_2 = 0;
	uint32_t both = 0;
	uint32_t none = 0;
	uint32_t eq = 0;
	uint32_t better_1 = 0;
	uint32_t better_2 = 0;

	while(getline(fq_in, read_id)){
		read_id = read_id.substr(1);
		front(read_id);
		
		
		vector<string> f1 = string_split(ali1,'\t');front(f1[0]);
		vector<string> f2 = string_split(ali2,'\t');front(f2[0]);
		
		
		//skip doubled lines
		while(f1[0] == last_read_id && getline(sam1_in, ali1) ){
			f1 = string_split(ali1,'\t'); front(f1[0]);
		}
		while(f2[0] == last_read_id && getline(sam2_in, ali2) ){
			f2 = string_split(ali2,'\t'); front(f2[0]);
		}
		
		front(read_id);
		
		if(f1[0] == f2[0] && f2[0] == read_id && !( (stoi(f1[1]) & 4) || (stoi(f2[1]) & 4 ))){
			both++;
			uint32_t e1,e2 = 0;
			
			//find dist
			for(uint32_t i = 0; i < f1.size(); i++)
				if(f1[i].compare(0,5,"NM:i:") == 0)
					e1 = stoi (f1[i].substr(5));
			//find dist
			for(uint32_t i = 0; i < f2.size(); i++)
				if(f2[i].compare(0,5,"NM:i:") == 0)
					e2 = stoi (f2[i].substr(5));
			
			if(e1 == e2) eq++; else if (e1 > e2) better_2++; else better_1++;
		} else if (f2[0] == read_id &&  !(stoi(f2[1]) & 4)){
			only_2++;
		} else if (f1[0] == read_id && !(stoi(f1[1]) & 4)){
			only_1++;
		} else {
			none++;
		}
		last_read_id = read_id;
		getline(fq_in, read_id); //read
		getline(fq_in, read_id); //+
		getline(fq_in, read_id); //qual
	}
	
	cout << "sam 1: " << argv[2] << endl;
	cout << "sam 2: " << argv[3] << endl;
	cout << "only in sam 1   "  << only_1 << endl;
	cout << "only in sam 2   "  << only_2 << endl;
	cout << "in both         "  << both << endl;
	cout << "in none         "  << none << endl;
	cout << "equal in both   "  << eq << endl;
	cout << "better in sam 1 " << better_1 << endl;
	cout << "better in sam 2 " << better_2 << endl;
}
