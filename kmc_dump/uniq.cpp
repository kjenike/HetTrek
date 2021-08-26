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
#include "uniq.hpp"
#include <unordered_map>
#include <unordered_set>

using namespace std; 

int vals[256];

void always_call_this()   
{
    for (int i=0; i<256; i++) vals[i]=0;
    vals['A']=0;
    vals['C']=1;
    vals['G']=2;
    vals['T']=3;
}

void reverse_c(string& r)
{
	//string r_string;
	reverse(r.begin(), r.end());
	
	//https://stackoverflow.com/questions/33074574/creating-complement-of-dna-sequence-and-reversing-it-c
	
	
	for (std::size_t i = 0; i < r.length(); ++i){
        	switch (r[i]){
        	case 'A':
            		r[i] = 'T';
            		break;    
        	case 'C':
            		r[i] = 'G';
            		break;
        	case 'G':
            		r[i] = 'C';
            		break;
        	case 'T':
            		r[i] = 'A';
           		break;
        }
    }
	//return s; 
}

string to_numbers(string& s)
{	
	//int the_number = 0;
	string tmp_s = s;
	for (int i = 0; i < s.size(); i++) {
		switch (s[i]){
                	case 'A':
                        	tmp_s[i] = '0';
                        	break;
                	case 'C':
                        	tmp_s[i] = '1';
                        	break;
                	case 'G':
                        	tmp_s[i] = '2';
                        	break;
                	case 'T':
                        	tmp_s[i] = '3';
                        	break;
			default :
				tmp_s[i] = '4';
				break;
		}
	}
	return tmp_s;
}

void find_kmers(string& read, vector<string>& read_kmers, int k)
{
	//vector<string> read_kmers;
	//std::hash<std::string> h;	
	for (int x=0; x < (read.length()-k+1); x++)
	{
		string kmer = read.substr(x, k);
		//string kmer; kmer = to_numbers(kmer_tmp);
		//kmer_tmp.clear();
		//vector<string> tmp_arr = 
		read_kmers.push_back(kmer);
		//Come back here and add the RC of the string. Right now that seems too complicated 
		//string r_kmer = reverse_c(kmer);
		//read_kmers.push_back(r_kmer);
		
	}

	//return read_kmers;
}

void find_kmers_local(string& read, vector<string>& read_kmers, int k, std::unordered_map<string, uint32_t>& kmer_db)
{
        //vector<string> read_kmers;

        for (int x=0; x < (read.length()-k+1); x++)
        {
                string kmer; kmer = read.substr(x, k);
                read_kmers.push_back(kmer);
                //Come back here and add the RC of the string. Right now that seems too complicated 
                string r_kmer = kmer;
	       	reverse_c(r_kmer);
                read_kmers.push_back(r_kmer);
                kmer_db[r_kmer] = 0;
                kmer_db[kmer] = 0;
        }

        //return read_kmers;
}

void find_mins(string& read, vector<uint32_t>& cnts, vector<int>& line_mins, int& k, int& m, int& cov_t)
{
	std::hash<std::string> h;
	//cout << "Looking for minimizers\n";
	//vector<string> line_mins;
	//line_mins.push_back(0);
	for (size_t x=0; x <= (read.length() - k) ; x+=15)
	{
		string forward; 
		forward = read.substr(x, k);
		//static_cast<std::uint32_t>(min);
		int min = 1000000000; //h("Vvvvvoyager");
		//static_cast<std::uint32_t>(min);
		for (int j=0; j < k-m ; j++)
		{
			int sub; sub=static_cast<int>(h(forward.substr(j, m))) ;
			if (sub < min)
			{
				//std::cout << "At least the hashed value is smaller :) \n" ;
				//std::cout << cnts[x] << "\n";
				//std::cout.flush();
				if (cnts[x] < cov_t)
				{
					if (cnts[x] > 3) {
						//std::cout << "At least the hashed value is smaller :) \n" ;
                                		//std::cout << cnts[x] << "\n";
                                		//std::cout.flush();
						min = sub;
					}
			
				}
			}
			
		}
		if (min != 1000000000 ) { //static_cast<std::uint32_t>(h("Vvvvvoyager"))){
			//std::cout << "We changed the min value! \n";
			//std::cout.flush();
			
			//if (std::find(line_mins.begin(), line_mins.end(), min)==line_mins.end()){
				
				//std::cout << "And we addedthe min value to the line mins\n" ;
				//std::cout.flush();
				line_mins.push_back(min);
			//}
			//line_mins.push_back(min); //Changed insert to pushback
		}
	}
	//std::cout << line_mins.size() << "\n" ;
	//std::cout.flush();

	//sort(line_mins.begin(), line_mins.end());
	/*if (line_mins.size() < 2) {
		std::cout << "THERE ARE NO MINIMIZERS :( \n";
		std::cout.flush();
	}*/

	//return line_mins;
}

void update_kmerdb(string& kmer, unordered_map<string, vector<uint32_t>>& kmer_db, int read_number)
{
	if (kmer_db.count(kmer))
	{
		vector<uint32_t> tmp = kmer_db[kmer];
		tmp.push_back(read_number);
		kmer_db[kmer] = tmp;
	}
	else 
	{
		vector<uint32_t> tmp = {read_number};
		kmer_db[kmer] = tmp;
	}
	
}

void make_kmer_db(vector<string>& seqs, unordered_map<string, vector<uint32_t>>& kmer_db, int k)
{
	//map<string, int> kmer_db;
	//What about finding all of the kmers, then sorting, then countint up all of  the kmers? 	
	/*for (int x=0; x < seqs.size(); x++){
		vector<int> read_kmers;
		find_kmers(seqs[x], read_kmers, k);
		//cout << "Found the kmers!\n";
	
		for (int j=0; j < read_kmers.size() ; j++){
			//cout << j << "\n";
			update_kmerdb(read_kmers[j], kmer_db, x);
		}
	}*/
	/*vector<string> read_kmers;
	for  (int x=0; x < seqs.size(); x++){
		for (int y=0; y<(seqs[x].size()-k+1), y++){
			//
			string kmer; kmer = seqs[x].substr(y, k);
                	read_kmers.push_back(kmer);
                	//Come back here and add the RC of the string. Right now that seems too complicated
                	string r_kmer = reverse_c(kmer);
               		read_kmers.push_back(r_kmer);
		}
	}
	//Sort
	sort(read_kmers.begin(), read_kmers.end())	
	//Go through and make the counts 
	string current_kmer;
	string last_kmer = read_kmers[0];
	int cntr=0;
	for (int j=0; j<read_kmers.size(); j++){
		if (read_kmers[j] == last_kmer)
	}*/
	//return kmer_db;
}

void get_counts(string& read, unordered_map<string, vector<uint32_t>>& kmer_db, vector <uint32_t>& cnts_q, int k)
{
	//vector <int> cnts_q;
	for (int x=0; x < (read.length() - k + 1  ) ; x++){
		string kmer;
		kmer = read.substr(x, k);

		int cnt = kmer_db[kmer].size();
		cnts_q.push_back(cnt);
	}
	/*rev_read = reverse_c(read);
	for (int x=0; x < (rev_read.length() - k + 1) ; x++){
                string kmer;
                kmer = rev_read.substr(x, k);

                int cnt = kmer_db[kmer].size();
                cnts_q.push_back(cnt);
        }*/
	//return cnts_q;
}

void find_intersection(const vector<uint32_t>& kmer, vector <int>& shared_reads, int& cnt){

	cnt = 0;
	for (uint32_t ele: shared_reads){
		if (std::binary_search (kmer.begin(), kmer.end(), ele)){
                        ++cnt;
                }
	}
	/*for (auto ele: shared_reads){
		std::cout << ele << "\n";
	}*/
	



}

uint64_t ReverseComp(const uint64_t mer, uint8_t kmerSize)
{
    uint64_t res = ~mer;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);

    return (res >> (2 * (32 - kmerSize)));
}

void characterT0Bits(char character)
{
    int vals[256];
    for (int i=0; i<256; i++) vals[i]=0;
    vals['A']=0;
    vals['C']=1;
    vals['G']=2;
    vals['T']=3;
}

int kmerize(const string& s, int k)
{
        int kmer = 0;
        for(int i = 0; i<k; i++) kmer = (kmer << 2) | vals[s[i]];
        return kmer;
}


void turn_kmer_string_to_int(string& kmer, int& kmer_int, int k)
{
    //int kmer_int = 0;
    for (int y=0; y< kmer.length(); ++y) {
        char character = kmer[y];
        int bits = 0;//kmerize(character, k);

        kmer_int = kmer_int << 2;
        kmer_int = kmer_int|bits;
    }
}

void get_local_counts(string& read, const std::unordered_map<int, vector<int>>& kmer_db, vector<uint32_t>& v, int k)
{
	for (int x=0; x < (read.length() - k + 1  ) ; x++){
		//string kmer_tmp = read.substr(x, k);
		//string kmer = to_numbers(kmer_tmp);
		//int kmer_suf = kmerize(kmer_tmp, k) ; //0;
        	int kmer_pre = kmerize((read.substr(x,23)), 23) ;
		int kmer_suf = kmerize((read.substr((x+23),8)), 8) ;
		
        	int kmer_pre_r = ReverseComp(kmer_pre, 23) ; //0;
        	int kmer_suf_r = ReverseComp(kmer_suf, 8) ;
		
		uint32_t cnt1 = 0;
        	uint32_t cnt2 = 0;
		//int this_len = kmer.size();
		auto it = kmer_db.find(kmer_pre);
		if (it != kmer_db.end()){
			int up = std::upper_bound(it->second.begin(), it->second.end(), kmer_suf) - it->second.begin();
			int low = std::lower_bound(it->second.begin(), it->second.end(), kmer_suf) - it->second.begin();
			cnt1 = up-low ;
		} 
		it = kmer_db.find(kmer_pre_r);
		if (it != kmer_db.end()) {
			int up = std::upper_bound(it->second.begin(), it->second.end(), kmer_suf_r) - it->second.begin();
			int low =  std::lower_bound(it->second.begin(), it->second.end(), kmer_suf_r) - it->second.begin();
			cnt2 = up-low ;
		}
		
		/*if (kmer_db.count(kmer_pre)) {
			int up = std::upper_bound(kmer_db.at(kmer_pre).begin(), kmer_db.at(kmer_pre).end(), kmer_suf) - kmer_db.at(kmer_pre).begin();
			int low = std::lower_bound(kmer_db.at(kmer_pre).begin(), kmer_db.at(kmer_pre).end(), kmer_suf) - kmer_db.at(kmer_pre).begin();
			cnt1 = up-low ;//count(kmer_db.at(kmer_pre).begin(), kmer_db.at(kmer_pre).end(), kmer_suf);
        	}
       		if (kmer_db.count(kmer_pre_r)) {
            		int up = std::upper_bound(kmer_db.at(kmer_pre_r).begin(), kmer_db.at(kmer_pre_r).end(), kmer_suf_r) - kmer_db.at(kmer_pre_r).begin();
			int low =  std::lower_bound(kmer_db.at(kmer_pre_r).begin(), kmer_db.at(kmer_pre_r).end(), kmer_suf_r) - kmer_db.at(kmer_pre_r).begin();
			cnt2 = up-low ;//count(kmer_db.at(kmer_pre_r).begin(), kmer_db.at(kmer_pre_r).end(), kmer_suf_r);
        	}*/

		v.push_back(cnt1+cnt2);
	}
}

void make_local_kmer_db( std::vector<int>& shared_reads, vector<string>& seqs, std::unordered_map<int, vector<int>>& kmer_db, int k)
{
    
	//https://stackoverflow.com/questions/1204313/counting-occurrences
	for ( int x: shared_reads){ //(int x=0; x < seqs.size(); x++){
		
        	//vector<string> read_kmers;
        	//find_kmers(seqs[x], read_kmers, k);
		
		for (int y=0; y < (seqs[x].length()-k+1); y++) {
			
            		int kmer_pre_int = kmerize((seqs[x].substr(y,23)), 8) ; //0;
			int kmer_suf_int = kmerize((seqs[x].substr((y+23),8)), 8) ;

            		if (kmer_db.count(kmer_pre_int)) {
                		//int tmp_cntr = kmer_db.at(kmer_int);
                		//++tmp_cntr;
                		kmer_db[kmer_pre_int].push_back(kmer_suf_int); // = tmp_cntr;
            		} else {
                		kmer_db[kmer_pre_int] = {kmer_suf_int};   
            		} 
		}
	}

	for (auto it : kmer_db) {
    		sort(it.second.begin(), it.second.end());
	}
}



vector<string> split(string line, char delim)
{
	//https://stackoverflow.com/questions/20755140/split-string-by-a-character
	//cout << line << "\n";
	vector<string> arr_split;
	//replace(line.begin(), line.end(), delim, ' '); 
	//cout << line << "\n";
	stringstream ss(line);
	string item;

	while (getline (ss, item, delim)){
		arr_split.push_back(item);
		//cout << item << "\n";
	}
	return arr_split;
}

/*int find_start(vector <string> s)
{
	int start;

	//vector <string> s; s = split(header, '|');
        string strand; strand = s[2];
        if (strand == "forward") {
		string tmp; tmp = s[3];
                start = std::stoi(tmp);
        } else {
		string tmp2; tmp2 = s[5];
                int len; len = std::stoi(tmp2);
		string tmp3; tmp3 = s[3];
                int tmp4; tmp4 = std::stoi(tmp3);
                int stop; stop = 1000000 - tmp4;
		start = stop - len;
        }

	return start;
}

int find_stop(vector <string> s)
{
        int stop;
	//vector <string> s; s = split(header, '|'); 
	string strand; strand = s[2];
	if (strand == "forward") {
		string tmp1; tmp1 = s[3];
		int start; start = std::stoi(tmp1);
		string tmp2; tmp2 = s[5];
		int len; len = std::stoi(tmp2); //TODO 
		stop = start + len;
	} else { 
		string tmp3; tmp3 = s[5];
		int len; len = std::stoi(tmp3);
		string tmp4; tmp4 = s[3] ;
		int start; start = std::stoi(tmp4);
		stop = 1000000 - start;
	}

        return stop;
}

void find_pos_controls(int read_number, vector<string>& que_names, vector <string>& seqs, int ref_start, int ref_stop, unordered_map<string, vector<uint32_t>>& kmer_db, hash<string>& h, vector<string>& pos_names, int& k, int& m, int& cov_t, vector<set <int> >& pos_mins, vector<int>& pos_overlaps, vector<int>& read_starts, vector<int>& read_stops) 
{
	//cout << que_names[0] << "\n";
	//cout << "Find Pos Controls \n";
	pos_overlaps.clear();
	pos_names.clear();
	pos_mins.clear();
	for (int i=0; i < que_names.size() ; i++){
		//
		string r; r = que_names[i];
		vector<string> rd; rd=split(r, '|');
		//cout << "Split the read name \n";
		string this_strand=rd[2]; //TODO
		//cout << "Split the strand \n";
		//cout << r;
		//cout << rd[5] << "\n";
		vector<string> this_l = split(rd[5], ' '); //TODO
		
		//cout << r;
		string this_len = this_l[0] ;
		//cout << "Found the length \n";
		int this_start;
		int this_stop;
		//cout << "About to check the strand \n";
		if (this_strand == "forward") { 
			//
			this_start = std::stoi(rd[3]);
			this_stop = this_start + std::stoi(this_len);
		} else {
			//
			this_stop = 1000000 - std::stoi(rd[3]);
			this_start = this_stop - std::stoi(this_len);
		}
		//cout << "Found the strand \n";
		//Calculate the overlap 
		if (ref_start <= this_start) {
			if (this_start <= ref_stop) {
				//
				int overlap; overlap = ref_stop - this_start;
				pos_overlaps.push_back(overlap);
				pos_names.push_back(r);
				vector<uint32_t> cnts_q;
				get_counts(seqs[i], kmer_db, cnts_q, k);
				vector<uint32_t> this_min; 
				find_mins(seqs[i], cnts_q, this_min, h, k, m, cov_t);
				//pos_mins.push_back(this_min);
				//read_stops.push_back(this_stop);
				//read_starts.push_back(this_start);
			}
		} else if (ref_start <= this_stop) {
			if (this_stop <= ref_stop) {
				//
				int overlap; overlap = this_stop - ref_start;
				pos_overlaps.push_back(overlap);
				pos_names.push_back(r);
				vector<uint32_t> cnts_q;
				get_counts(seqs[i], kmer_db, cnts_q, k);
				vector<uint32_t> this_min;
				find_mins(seqs[i], cnts_q, this_min, h, k, m, cov_t);
				//pos_mins.push_back(this_min);
				//read_stops.push_back(this_stop);
				//read_starts.push_back(this_start);
			}
		}
		//cout << "Overlap calculated \n";
		
	}
	//cout << "Size of Pos's: " << pos_names.size() <<"\n";
	return void();
}
*/
int find_sims(vector <uint32_t>& s1, vector <uint32_t>& s2){
	
	//int total_sim; 
	//set<string> s1(v1.begin(), v1.end());
	//set<string> s2(v2.begin(), v2.end());
	//std::multiset<int> s3;
	//std::vector<int> s3;
	//set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(), std::inserter(s3, s3.begin())); //std::back_inserter(s3));
	//set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(), std::back_inserter(s3));
	int cnt = 0;
        for (uint32_t ele: s2){
		//if (std::find(s2.begin(), s2.end(), ele) != s2.end()) { //(shared_reads.count(ele)){
                //        ++cnt;
                //}
		//Binary search, this should be the fastest? 
        	if (std::binary_search (s1.begin(), s1.end(), ele)){	
			++cnt;
		}
	}
	
	
	//total_sim = s3.size();
	
	return cnt;
}

/*void find_shared_reads(int t, std::vector<int>& shared_reads, std::vector<std::string>& seqs, std::vector<uint32_t>& ref_counts,std::vector<uint32_t>& rev_counts , int k, std::hash<std::string>& h, int m, int cov_t, int read_number, std::vector <vector <int>>& min_db, std::string& read, std::string& read_rev)
{
	//cout << "Looking for shared reads\n";
        //vector <uint32_t> ref_cnts;
        vector<int> forw = min_db[read_number];
        string r;
        //vector <uint32_t> rev_cnts;
        vector <int> ref_mins;
	//find_mins(read, ref_counts, forw, h, k, m, cov_t);
        
	
	find_mins(read_rev, rev_counts, ref_mins, k, m, cov_t);
	//cout << "Found minimizers R \n";
        
	//What if we just found the shared reads independently?
	for (int x=0; x < seqs.size()  ; x++){
		uint32_t similar_f = find_sims(min_db[x], forw);
		uint32_t similar_r = find_sims(min_db[x], ref_mins);
		if ((similar_f+similar_r) > t) {
			shared_reads.push_back(x);
			//std::cout << x << "\n" ;
		}
		

	}
	
	
	ref_mins.insert(ref_mins.end(), forw.begin(), forw.end()); 
	//Sort the ref_mins 
	sort(ref_mins.begin(), ref_mins.end());

        //vector<string> above_t_names;
	for (int x=0; x < seqs.size()  ; x++){
                int similar = 0;
                similar = find_sims(min_db[x], ref_mins);
                //que_sim.push_back(similar);
                if (similar > t) {
			shared_reads.push_back(x);
                }
        }

	return void();
}*/

/*void compare_q(int read_number, map<string, vector<uint32_t>>& kmer_db, vector <string>& seqs, vector <string>& que_names, vector <set <int>>& min_db, hash<string>& h, vector<string>& pos_names, int& k, int& m, int& cov_t, vector<set <int> >& pos_mins, vector<int>& pos_overlaps, vector<int>& read_starts, vector<int>& read_stops)
{
	//Make the reference min DB 
	//Make the references constant 
	//Also, make these references ^^^^^^
	map<int, int> minDB_ref;
	int i; i=0;
	vector <string> all_reads;
	vector <uint32_t> ref_cnts;
	vector<uint32_t> forw;
	string r;
	vector <uint32_t> rev_cnts;
	vector <uint32_t> ref_mins;
	int ref_start;
	int ref_stop;
	vector <string> ref_header;
	string ref_header_tmp;
	
	//Get the ref start and stop (from the que names) 
	ref_header_tmp = que_names[read_number];
	ref_header = split(ref_header_tmp, '|');
	//cout << ref_header_tmp << "\n";
	ref_start = find_start(ref_header);
	ref_stop = find_stop(ref_header);

	//Now we need to get the counts for the reference
	get_counts(seqs[read_number], kmer_db, ref_cnts, k);
	find_mins(seqs[read_number], ref_cnts, forw, h, k, m, cov_t);
	
	r = reverse_c(seqs[read_number]);
	get_counts(r, kmer_db, rev_cnts, k);
	find_mins(r, rev_cnts, ref_mins, h, k, m, cov_t);

	ref_mins.insert(ref_mins.end(), forw.begin(), forw.end()); //ref_mins.insert( ref_mins.end(), forw.begin(), forw.end() ); //ref_mins.insert(forw.begin(), forw.end());
	

	find_pos_controls(read_number, que_names, seqs, ref_start, ref_stop, kmer_db, h, pos_names, k, m, cov_t, pos_mins, pos_overlaps, read_starts, read_stops);

	//cout << "Positive file processed\n";
	
	//Find the number of similarities btwn the reads 
	vector<int> que_sim;
	vector<int> que_sim_true;
	int tp=0;
	int fp=0;
	int fn=0;
	int t=75;

	vector<string> above_t_names;
	for (int x=0; x < que_names.size() ; x++){
		int similar;
		similar = find_sims(min_db[x], ref_mins);
		//que_sim.push_back(similar);
		if (similar > t) {
			if (std::find(pos_names.begin(), pos_names.end(), que_names[x]) != pos_names.end()){
				//This is a TP!
				tp = tp+1;
				//cout << "TP similars: " << similar << "\n";
                                vector <string> tmp; tmp = split(que_names[x], '|');
                                above_t_names.push_back(tmp[0]);
			}else{
				//This is an FP :( 
				fp = fp+1;
				//cout << "FP similars: " << similar << "\n";
                                vector <string> tmp; tmp = split(que_names[x], '|');
                                above_t_names.push_back(tmp[0]);
			}
			
		}
	}
	//cout << "Found similarities \n";
	//cout << "Found TPs and FPs  \n";
	//cout << pos_names[0] << "\n" ;	
	for (int x=0; x < pos_names.size(); x++) {
		//Now we'll look for FNs 
		vector <string> tmp; tmp = split(pos_names[x], '|');
		if (std::find(above_t_names.begin(), above_t_names.end(), tmp[0]) == above_t_names.end()){
		//if (pos_names[x] not in above_t_names) {
			//
			if (pos_overlaps[x] > 1000) {
				//Want to make sure we are considering only stuff that is over x bp overlap
				fn = fn + 1;

			}
		}	

	}

//	cout << read_number << "\n" ;
//	cout << "TPs: " << tp << "\n";
//	cout << "FPs: " << fp << "\n";
//	cout << "FNs: " << fn << "\n";
//	cout << "****************************************\n" ;
	return void();
}
*/
/*int main()
{
	vector<string> pos_names;
	vector<set <int> > pos_mins;
	vector<int> pos_overlaps;
	vector<int> read_starts;
	vector<int> read_stops;
	int k;
	int m;
	int cov_t;

	Timer t;
	//cout << "Hello World!\n";
      	k = 31;
	m = 25; 
	cov_t = 35;
	hash<string> h;
	std::vector <int> read_number(3000);
	std::iota (std::begin(read_number), std::end(read_number), 0);

	std::cout << t.lap() << " \n";
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
	std::cout << t.lap() << " Read in Fasta file\n";
	
	map<string, vector<int>> kmer_db;
	make_kmer_db(seqs, kmer_db);
	
	std::cout << t.lap() << " Made the Kmer DB\n";
	//cout << "Kmer DB made\n";
	
	vector <set <int>> min_db;
	//Now it's time to find the minimizers 
	for (int x=0; x < seq_cntr; x++){
		//Get the counts 
		vector <int> cnts_q;
		get_counts(seqs[x], kmer_db, cnts_q);
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
	for (int x=0; x < read_number.size(); x++){
		compare_q(read_number[x], kmer_db, seqs, que_names, min_db, h, pos_names, k, m, cov_t, pos_mins, pos_overlaps, read_starts, read_stops);	
	}
	std::cout << t.lap() << " \n";
	return 0;
}*/

