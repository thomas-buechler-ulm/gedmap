#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>      // std::setw

using namespace std;



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

	if(argc != 2)	cout << "1 Argument (fq-file) expected" << endl; 

	ifstream fq_in;
	fq_in.open(argv[1],ios_base::in);


	vector<uint32_t> errors = vector<uint32_t>(40,0);

	string entry;

	while(getline(fq_in, entry)){
		
		vector<string> fields  = string_split(entry, '_');
		uint32_t d = stoi(fields[5]);
		
		d<40? errors[d]++ : errors[39]++;
		
		getline(fq_in, entry); //query
		getline(fq_in, entry); //+
		getline(fq_in, entry); //qual
	}
	
	cout << "d:\n";
	for(auto i : errors)	cout << i << '\t';
	cout << endl;
	
	uint32_t sum = 0;
	uint32_t count = 0;
	
	for(uint32_t i = 0; i < 20; i++){
		count += errors[i];
		sum += i*errors[i];
	}
	
	double avg = sum/ (double)count;
	uint32_t med = 0;
	uint32_t com = 0;
	for(; com < (count/2); med++ ){
		com += errors[med];
	}
	

	cout	<< "median      " << setw(6)<< med-1 << endl
		<< "avg         " << setw(6)<< avg << endl;
}
