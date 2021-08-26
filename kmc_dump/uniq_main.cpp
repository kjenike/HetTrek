#include <iostream>
#include <string>
#include <fstream>
#include <map> 
#include <vector>
#include <tuple>
#include <algorithm>
#include <sstream>
#include <set>
#include <iterator>
#include <numeric>
#include "timerSK.hpp"
#include "uniq.hpp"

using namespace std;

int main()
{
	vector<string> pos_names;
	vector<vector <uint32_t> > pos_mins;
	vector<int> pos_overlaps;
	vector<int> read_starts;
	vector<int> read_stops;
	int k;
	int m;
	int cov_t;

	//TimerSK t;
	//cout << "Hello World!\n";
      	k = 31;
	m = 25; 
	cov_t = 35;
	hash<string> h;
	std::vector <int> read_number(3000);
	std::iota (std::begin(read_number), std::end(read_number), 0);

	//std::cout << t.lap() << " \n";
	string fasta_file_que; fasta_file_que="simulatedreads_template1_g1000000_cov30_rl10000_err1.0_het1.0_indel1.fasta";
	vector<string> seqs;
	vector<string> que_names;
	
	//Read in the fasta file 
	int seq_cntr; seq_cntr=0;
	int name_cntr; name_cntr=0;
	std::ifstream file(fasta_file_que);
	if (file.is_open()) {
		//cout << "File open\n";
		std::string line;
		while (std::getline(file, line)) {
			//cout << line << endl;
			if (name_cntr > seq_cntr)
			{
				seqs.push_back(line); //seq_cntr] = line;
				seq_cntr++;
			}
			else
			{
				que_names.push_back(line);
				name_cntr++;
			}
		}
		file.close();
	}
	//cout << "Fasta file read\n" ;
	//std::cout << t.lap() << " Read in Fasta file\n";
	
	map<string, vector<uint32_t>> kmer_db;
	
	
	
	
	
	
	
	/*
	make_kmer_db(seqs, kmer_db, k);
	
	std::cout << t.lap() << " Made the Kmer DB\n";
	//cout << "Kmer DB made\n";
	
	vector <set <int>> min_db;
	//Now it's time to find the minimizers 
	for (int x=0; x < seq_cntr; x++){
		//Get the counts 
		vector <uint32_t> cnts_q;
		get_counts(seqs[x], kmer_db, cnts_q, k);
		//get the minimizers for that seq 
		set <int> line_mins;
		find_mins(seqs[x], cnts_q, line_mins, h, k, m, cov_t);
		//Save these minimizers 
		min_db.push_back(line_mins);
	}
	std::cout << t.lap() << " Found minimizers \n";
	//cout << "All minimizers found\n";
	//cout << "***********************\n";
	//Now we are going to go through and identify the reads that we think fit 
	for (int x=1; x < read_number.size(); x++){
		compare_q(read_number[x], kmer_db, seqs, que_names, min_db, h, pos_names, k, m, cov_t, pos_mins, pos_overlaps, read_starts, read_stops);	
	}
	std::cout << t.lap() << " \n";
	return 0;
	*/
}

