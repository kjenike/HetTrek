#ifndef UNIQ_HPP
#define UNIQ_HPP
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
#include <unordered_map>
#include <unordered_set>


using namespace std;

void always_call_this();

void reverse_c(string& r);

void find_kmers(string& read, vector<string>& read_kmers, int k);

void find_mins(string& read, vector<uint32_t>& cnts, vector<int>& line_mins,  int& k, int& m, int& cov_t);

void update_kmerdb(string& kmer, unordered_map<string, vector<uint32_t>>& kmer_db, int read_number);

void make_kmer_db(vector<string>& seqs, unordered_map<string, vector<uint32_t>>& kmer_db, int k);

void get_counts(string& read, unordered_map<string, vector<uint32_t>>& kmer_db, vector <uint32_t>& cnts_q, int k);

vector<string> split(string line, char delim);

int kmerize(const string& s, int k);

uint64_t ReverseComp(const uint64_t mer, uint8_t kmerSize);

void turn_kmer_string_to_int(string& kmer, int& kmer_int, int k);

void make_local_kmer_db(std::vector<int>& shared_reads, vector<string>& seqs, std::unordered_map< int, vector<int>>& kmer_db, int k);

void get_local_counts(string& read, const std::unordered_map< int, vector<int>>& kmer_db, vector <uint32_t>& cnts_q, int k);

void find_kmers_local(string& read, vector<string>& read_kmers, int k, std::unordered_map<string, uint32_t>& kmer_db);

int find_start(vector <string> s);

int find_stop(vector <string> s);

void find_pos_controls(int read_number, vector<string>& que_names, vector <string>& seqs, int ref_start, int ref_stop, unordered_map<string, vector<uint32_t>>& kmer_db, hash<string>& h, vector<string>& pos_names, int& k, int& m, int& cov_t, vector<set <int> >& pos_mins, vector<int>& pos_overlaps, vector<int>& read_starts, vector<int>& read_stops);

int find_sims(vector <int>& s1, vector <int>& s2);

void find_intersection(const vector<uint32_t>& kmer, vector <int>& shared_reads, int& cnt);

//void find_shared_reads(int t, std::vector<int>& shared_reads, std::vector<std::string>& seqs, std::vector<uint32_t>& ref_counts,std::vector<uint32_t>& rev_counts , int k, std::hash<std::string>& h, int m, int cov_t, int read_number, std::vector <vector <int>>& min_db, std::string& read, std::string& read_rev);

#endif
