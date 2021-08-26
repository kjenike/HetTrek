#include "./../kmc_api/kmc_file.h"
#include "./../kmc_api/kmer_api.h"
#include "nc_utils.h"
#include "stdafx.h"
#include <algorithm>
#include <bits/stdc++.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <mutex>
#include <numeric>
#include <stdlib.h>
#include <thread>
#include <unistd.h>
#include <vector>
#include "uniq.hpp"
#include <unordered_map>
#include <chrono>

//smooth kmcdb reads.fastq errremoved.fasta erredits.fasta errpaths.fasta hetremoved.fasta hetedits.fasta hetpath1.fasta hetpath2.fasta hetpath3.fasta hetpath4.fasta hetpath5.fasta hetpath6.fasta error_threshold het_threshold unique_threshold anchor_threshold allowed_err_fraction allowed_rep_fraction max_nodes_to_search distance_multiplier strict > counts.txt

int get_type(uint32_t& coverage, int& error_threshold, int& het_threshold, int& unique_threshold)
{
	if (coverage <= error_threshold)
	{
		return 0;
	}
	else if (coverage <= het_threshold)
	{
		return 1;
	}
	else if (coverage <= unique_threshold)
	{
		return 2;
	}
	else
	{
		return 3;
	}
}

std::vector<std::string> get_adjacent(CKMCFile& file, std::vector<int>& shared_reads, const unordered_map<auto, vector<string>>& local_kmer_db, std::string& kmer, int& error_threshold, int& het_threshold, int& unique_threshold, bool& going_right, int k)
{
	//std::cout << "Getting adjacent\n";
	//std::cout.flush();

	std::vector<uint32_t> v;
	std::vector<std::string> adjacent_kmers;
	//for all possible nucleotide extensions from kmer
	for (char const &c: "ACGT")
	{
		
		std::string adjacent_kmer;
		if (going_right)
		{
			
			adjacent_kmer = kmer.substr(1)+c;
		}
		else
		{
			//cout << "Going left\n";
			//cout << kmer << "\n";
			//cout << kmer.length() << "\n";
			//cout << c << "\n";
			adjacent_kmer = c+kmer.substr(0, kmer.length()-1);
		}
		//std::cout << kmer << "\n" ; 
		//std::cout.flush();
		v.clear();
		if (shared_reads.size() > 30) {		
			get_local_counts(adjacent_kmer, local_kmer_db, v, k);//TODO
		} else {
			file.GetCountersForRead(adjacent_kmer, v);
		}
		if (v.size() == 0){
			v.push_back(0);
		}
		//std::cout << v.size() << "\n" ;
		//std::cout.flush();
		//std::cout << kmer << "\t" << v[0] << "\n" ; 
		/*for (uint32_t something: v) {
			std::cout << something << "\n" ;
		}*/
		int current_type = get_type(v[0], error_threshold, het_threshold, unique_threshold);
		
		//if the adjacent kmer is not an error
		if ((current_type > 0)) //&& (current_type < 3)) //Added the 2nd condition 
		{
		  adjacent_kmers.push_back(adjacent_kmer);
		}
	}
	return adjacent_kmers;
}

bool is_left_anchor(CKMCFile& file, std::vector<int>& shared_reads, std::string& previous_kmer, int& previous_count, int& k, const unordered_map<auto, vector<string>>& local_kmer_db, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold)
{
	if ((previous_kmer.length() != k) || (previous_count > unique_threshold)) // > unique_threshold or anchor_threshold
	//if (previous_kmer.length() != k)
	{
		return false;
	}
	bool going_right = true;
	std::vector<std::string> adjacent_kmers = get_adjacent( file, shared_reads, local_kmer_db, previous_kmer, error_threshold, het_threshold, unique_threshold, going_right, k);
	if (adjacent_kmers.size() == 2)
	{
		std::vector<uint32_t> v1;
		//file.GetCountersForRead(adjacent_kmers[0], v1);
		if (shared_reads.size() > 30) {
			get_local_counts(adjacent_kmers[0], local_kmer_db, v1, k);
		} else {
                        file.GetCountersForRead(adjacent_kmers[0], v1);
                }
		int count1 = v1[0];
		std::vector<uint32_t> v2;
		//file.GetCountersForRead(adjacent_kmers[1], v2);
		if (shared_reads.size() > 30) {
			get_local_counts(adjacent_kmers[1], local_kmer_db, v2, k);
		} else {
                        file.GetCountersForRead(adjacent_kmers[1], v2);
                }
		int count2 = v2[0];
		//If the sum of the coverages of the two branches is within 3 of the homozygous portion
		if ((previous_count - 3 <= count1 + count2) && (count1 + count2 <= previous_count + 3))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

bool is_right_anchor(CKMCFile& file, std::vector<int>& shared_reads, std::string& current_kmer, int& current_count, int& k, const unordered_map<auto, vector<string>>& local_kmer_db, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold)
{
	if ((current_kmer.length() != k) || (current_count > unique_threshold)) // > unique_threshold or anchor_threshold?
	//if (current_kmer.length() != k)
	{
		return false;
	}
	bool going_right = false;
	std::vector<std::string> adjacent_kmers = get_adjacent(file, shared_reads, local_kmer_db, current_kmer, error_threshold, het_threshold, unique_threshold, going_right, k);
	if (adjacent_kmers.size() == 2)
	{
		std::vector<uint32_t> v1;
		//file.GetCountersForRead(adjacent_kmers[0], v1);
		if (shared_reads.size() > 30) {
			get_local_counts(adjacent_kmers[0], local_kmer_db, v1, k);
		} else {
                        file.GetCountersForRead(adjacent_kmers[0], v1);
                }
		int count1 = v1[0];
		std::vector<uint32_t> v2;
		//file.GetCountersForRead(adjacent_kmers[1], v2);
		if (shared_reads.size() > 30) {
			get_local_counts(adjacent_kmers[1], local_kmer_db, v2, k);
		} else {
                        file.GetCountersForRead(adjacent_kmers[1], v2);
                }
		int count2 = v2[0];
		//If the sum of the coverages of the two branches is within 3 of the homozygous portion
		if ((current_count - 3 <= count1 + count2) && (count1 + count2 <= current_count + 3))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

int get_type_het(CKMCFile& file, std::vector<int>& shared_reads, int& previous_type, std::string& previous_kmer, std::string& current_kmer, int& previous_count, int& current_count, int& k, const unordered_map<auto, vector<string>>& local_kmer_db, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold, std::string& anchor_found)
{
	//If "previously" at the beginning of the read
	if (previous_type == -1)
	{
		//if current kmer is a right anchor
		if ((current_count > (previous_count + error_threshold)) && is_right_anchor(file,shared_reads, current_kmer, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		//((current_count > (previous_count + error_threshold)) && is_right_anchor(shared_reads, current_kmer, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			
			anchor_found = "right";
			return 2;
		}
		//if current kmer is a left anchor
		else if (is_left_anchor(file,shared_reads, current_kmer, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			
			anchor_found = "left";
			return 2;
		}
		//To keep it consistent with previous get_type, we need type to be -1 only before we start the read
		//So, in this case we will do what was done before, using the coverage and coverage thresholds to
		//determine whether the kmer is homozygous or not
		else if ((het_threshold < current_count) && (current_count <= unique_threshold))
		{
			anchor_found = "none";
			return 2;
		}
		else
		{
			anchor_found = "none";
			return 1;
		}
	}
	//If previously in a nonhomozygous region
	else if (previous_type == 1)
	{
		//If current kmer is a right anchor
		if ((current_count > (previous_count + error_threshold )) && is_right_anchor(file,shared_reads, current_kmer, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		//((current_count > (previous_count + error_threshold )) && is_right_anchor(shared_reads, current_kmer, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			//we have left the nonhom region
			
			anchor_found = "right";
			return 2;
		}
		//If current kmer is a left anchor
		else if (is_left_anchor(file,shared_reads, current_kmer, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			//this is a weird case where we find two left anchors in a row before a right anchor
			
			anchor_found = "left";
			return 2;
		}
		//else current kmer is not an anchor
		else
		{
			//we continue the nonhom region
			anchor_found = "none";
			return 1;
		}
	}
	//If previously in a homozygous region
	//else if (previous_type == 2)
	else
	{
		//If current kmer is a right anchor
		//if ((current_count > (previous_count + error_threshold)) && is_right_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		//{
			//this is a weird case where we find two right anchors in a row before a left anchor
			//we don't need to set anchor_found or return 2 because we really care
			//if current kmer is a left anchor, ending the homozygous region
			
		//}
		//If current kmer is a left anchor
		if (is_left_anchor(file,shared_reads, current_kmer, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			
			anchor_found = "left";
			return 2;
		}
		//If previous kmer is a left anchor
		else if ((current_count < (previous_count - error_threshold )) && is_left_anchor(file,shared_reads, previous_kmer, previous_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		//((current_count < (previous_count - error_threshold )) && is_left_anchor(shared_reads, previous_kmer, previous_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
		{
			//we have left the hom region
			
			anchor_found = "left";
			return 1;
		}
		//else previous kmer is not a left anchor (and current kmer is not a left anchor)
		else
		{
			//we continue the hom region
			anchor_found = "none";
			return 2;
		}
	}
}

std::vector<std::string> get_paths(CKMCFile& file, std::vector<int>& shared_reads, const unordered_map<auto, vector<string>>& local_kmer_db, int& error_threshold, int& het_threshold, int& unique_threshold, std::string& left_anchor_kmer, std::string& right_anchor_kmer, int& min_distance_of_path, int& max_distance_of_path, int& max_nodes_to_search, int& k, bool& queue_broken)
{
	//std::cout << "Getting paths\n";
	//std::cout.flush();

	std::string starting_anchor_kmer;
	bool going_right;
	
	if (left_anchor_kmer.empty())
	{
		//We are at the beginning of the read.
		//This function will find paths starting from right_anchor_kmer continuing left to the
		//beginning of the read where each kmer in the path is a nonerror kmer.
		starting_anchor_kmer = right_anchor_kmer;
		going_right = false;
	}
	else if (right_anchor_kmer.empty())
	{
		//We are at the end of the read.
		//This function will find paths starting from left_anchor_kmer continuing right to the
		//end of the read where each kmer in the path is a nonerror kmer.
		starting_anchor_kmer = left_anchor_kmer;
		going_right = true;
	}
	else
	{
		//We are in the middle of the read.
		//This function will find paths starting from left_anchor_kmer continuing right until
		//right_anchor_kmer where each kmer in the path is a nonerror kmer.
		starting_anchor_kmer = left_anchor_kmer;
		going_right = true;
	}
	//std::cout << "We are either going left or right\n";
	//std::cout.flush();
	
	//In every case, we follow all paths, but cut the depth of any path to max_distance_of_path.
	//We also ensure a minimum depth equal to min_distance_of_path.
	std::list<std::string> queue;
	queue.push_back(starting_anchor_kmer);
	
	//Initialize paths to store all the paths that are found.
	std::vector<std::string> paths;
	//We use i as a counter for how many nodes have been visited in the search.
	//If we haven't finished the search within max_nodes_to_search nodes, 
	//we break the search.
	//This drastically speeds up the run time for some regions.
	//Thankfully, it doesn't seem to impact effectiveness, since most searches
	//complete before this threshold.
	int i = 0;
	//This flag keeps track of whether we had to stop the search early.
	queue_broken = false;
	while(!queue.empty())
	{
		//std::cout << "Queue not empty\n";
		//std::cout.flush();	
		i++;
		std::string current_path = queue.front();
		std::string current_kmer;
		if (going_right)
		{
			
			current_kmer = current_path.substr(current_path.length()-k);
			
		}
		else
		{
			
			current_kmer = current_path.substr(0, k);
		}
		queue.pop_front();
		int current_depth = current_path.length()-k;
		//If we have to terminate search early
		if (i > max_nodes_to_search)
		{
			
			queue_broken = true;
			break;
		}
		//cout << "Checking the node depth\n";
		//If the depth of this node hasn't exceeded the max distance of the path
		if (current_depth <= max_distance_of_path)
		{
			//std::cout << "About to get adjacent\n";
			//std::cout.flush();
			//Extend the path by one nucleotide, keep the ones that are not error kmers
			std::vector<std::string> adjacent_kmers = get_adjacent(file,shared_reads, local_kmer_db, current_kmer, error_threshold, het_threshold, unique_threshold, going_right, k);
			//std::cout << "Extending the path by one nuc\n";
			//std::cout.flush();
			for (auto adjacent_kmer : adjacent_kmers)
			{
				std::string path;
				if (going_right)
				{
					path = current_path + adjacent_kmer.back();
				}
				else
				{
					path = adjacent_kmer.front() + current_path;
				}
				
				bool end_condition;
				//If we are in the middle of the read, we end when we have found a path
				//of nonerror kmers which bridges the anchor kmers
				//and doesn't terminate too early (i.e. before min_distance_of_path)
				if (!left_anchor_kmer.empty() && !right_anchor_kmer.empty())
				{
					end_condition = ((adjacent_kmer == right_anchor_kmer) && (current_depth + 1 >= min_distance_of_path));
				}
				//If we are at either end of the read, we end when we have found a path
				//of nonerror kmers which continues until the end of the read
				//and doesn't terminate too early (i.e. before min_distance_of_path)
				else
				{
					end_condition = ((current_depth + 1 == max_distance_of_path) && (current_depth + 1 >= min_distance_of_path));
				}
				if (end_condition)
				{
					//added this case for if the right and left anchors overlap
					if (!left_anchor_kmer.empty() && !right_anchor_kmer.empty() && (path.size() < 2*k))
					{
						int total_overlaps = 2*k - path.size();
						path.clear();
						for (int number_overlaps = 0; number_overlaps < total_overlaps; number_overlaps++)
						{
							path += "-";
						}
						paths.push_back(path);
					}
					else
					{
						if (!left_anchor_kmer.empty())
						{
							//remove left_anchor_kmer from path
							path.erase(path.begin(), path.begin()+k);
						}
						if (!right_anchor_kmer.empty())
						{
							//remove right_anchor_kmer from path
							path.erase(path.end()-k, path.end());
						}
						paths.push_back(path);
					}
				}
				// "No path yet\n";
				//Else we haven't found a path yet
				else
				{
					// "No path yet\n";
					queue.push_back(path);
				}
			}
		}
	}
	return paths;
}

int minDis(std::string s1, std::string s2, int n, int m, std::vector<std::vector<int>> &dp)
{
	// If any string is empty,
	// return the remaining characters of other string
	if (n==0)
	{
		return m;
	}
	if (m==0)
	{
		return n;
	}
	// To check if the recursive tree
	// for given n & m has already been executed	
	if (dp[n][m]!=-1)
	{
		return dp[n][m];
	}
	// If characters are equal, execute
	// recursive function for n-1, m-1
	
	if (s1[n-1] == s2[m-1])
	{
		if (dp[n-1][m-1] == -1)
		{
			return dp[n][m] = minDis(s1, s2, n-1, m-1, dp);
		}
		else
		{
			return dp[n][m] = dp[n-1][m-1];
		}
	}
	// If characters are nt equal, we need to
	// find the minimum cost out of all 3 operations.
	else
	{
		int m1, m2, m3; // temp variables
		
		if (dp[n-1][m]!=-1)
		{
			m1 = dp[n-1][m];
		}
		else
		{
			m1 = minDis(s1, s2, n-1, m, dp);
		}
		
		if (dp[n][m-1]!=-1)
		{
			m2 = dp[n][m-1];
		}
		else
		{
			m2 = minDis(s1, s2, n, m-1, dp);
		}
		
		if (dp[n-1][m-1]!=-1)
		{
			m3 = dp[n-1][m-1];
		}
		else
		{
			m3 = minDis(s1, s2, n-1, m-1, dp);
		}
		return dp[n][m] = 1 + std::min(m1, std::min(m2, m3));
	}
}

void write_error_paths(CKMCFile& file, std::vector<int>& shared_reads, int k, bool& queue_broken, std::vector<std::string>& edited_error_portions, std::ofstream& errpaths_queuecomplete0_numpaths0_output_file, std::ofstream& errpaths_queuecomplete0_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& errpaths_queuecomplete1_numpaths0_output_file, std::ofstream& errpaths_queuecomplete1_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& erredits_output_file, std::string& edited_read, int& read_number, int& first_error_idx, int& last_error_idx, std::string& before_first_error_kmer, std::string& original_error_portion, std::string& after_last_error_kmer, const unordered_map<auto, vector<string>>& local_kmer_db, std::mutex& outputfileMutex)
{
	//std::cout << "Writing error paths\n";
	//std::cout.flush();

	std::string left_part;
	if (!before_first_error_kmer.empty())
	{
		left_part = before_first_error_kmer.substr(1);
	}
	else
	{
		left_part = before_first_error_kmer;
	}
	std::string right_part;
	if (!after_last_error_kmer.empty())
	{
		right_part = after_last_error_kmer.substr(0, after_last_error_kmer.length()-1);
	}
	else
	{
		right_part = after_last_error_kmer;
	}
	std::size_t found = original_error_portion.find("-");
	std::string original;
	if (found!=std::string::npos)
	{
		original = left_part + right_part.substr(original_error_portion.size());
	}
	else
	{
		original = left_part + original_error_portion + right_part;
	}
	auto is_less_than = [&](std::string edited_error_portion1, std::string edited_error_portion2)
	{
		int n = original.length();
		std::size_t found1 = edited_error_portion1.find("-");
		std::string portion1;
		if (found1!=std::string::npos)
		{
			portion1 = left_part + right_part.substr(edited_error_portion1.size());
		}
		else
		{
			portion1 = left_part + edited_error_portion1 + right_part;
		}
		int m = portion1.length();
		std::vector<std::vector<int>> dp(n+1, std::vector<int>(m+1, -1));
		int dist1 = minDis(original, portion1, n, m, dp);
		std::size_t found2 = edited_error_portion2.find("-");
		std::string portion2;
		if (found2!=std::string::npos)
		{
			portion2 = left_part + right_part.substr(edited_error_portion2.size());
		}
		else
		{
			portion2 = left_part + edited_error_portion2 + right_part;
		}
		m = portion2.length();
		std::vector<std::vector<int>> dp2(n+1, std::vector<int>(m+1, -1));
		int dist2 = minDis(original, portion2, n, m, dp2);
		return dist1 < dist2;
	};
	
	std::sort(edited_error_portions.begin(), edited_error_portions.end(), is_less_than);
	
	std::ofstream* errwrite_output_file;
	bool was_edited = false;
	//we finished the search and presumably we have found one homozygous path or two heterozygous paths
	//if ((!queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
	if (edited_error_portions.size() >= 1)
	{
		//errwrite_output_file = &erredits_output_file;
		//std::cout << "Finished the search\n";
        	//std::cout.flush();
		found = edited_error_portions[0].find("-");
		if (found!=std::string::npos)
		{
			edited_read += left_part.substr(0, left_part.size() - edited_error_portions[0].size());
		}
		else
		{
			edited_read += left_part + edited_error_portions[0];
		}
		if (edited_error_portions[0] != original_error_portion)
		{
			was_edited = true;
		}
		//if (!before_first_error_kmer.empty())
		//{
		//	edited_read += before_first_error_kmer.substr(1) + edited_error_portions[0];
		//}
		//else
		//{
		//	edited_read += edited_error_portions[0];
		//}
	}
	//there are no paths, or there are more than two paths, or the search wasn't finished.
	//We are currently not editing.
	else
	{
		//errwrite_output_file = &errpaths_output_file;
		found = original_error_portion.find("-");
		if (found!=std::string::npos)
		{
			edited_read += left_part.substr(0, left_part.size() - original_error_portion.size());
		}
		else
		{
			edited_read += left_part + original_error_portion;
		}
		//if (!before_first_error_kmer.empty())
		//{
		//	edited_read += before_first_error_kmer.substr(1) + original_error_portion;
		//}
		//else
		//{
		//	edited_read += original_error_portion;
		//}
	}
	if ((!queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
	{
		errwrite_output_file = &errpaths_queuecomplete1_numpaths1to2_output_file;
	}
	if ((!queue_broken) && (edited_error_portions.size() < 1))
	{
		errwrite_output_file = &errpaths_queuecomplete1_numpaths0_output_file;
	}
	if ((!queue_broken) && (edited_error_portions.size() > 2))
	{
		errwrite_output_file = &errpaths_queuecomplete1_numpaths3plus_output_file;
	}
	if ((queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
	{
		errwrite_output_file = &errpaths_queuecomplete0_numpaths1to2_output_file;
	}
	if ((queue_broken) && (edited_error_portions.size() < 1))
	{
		errwrite_output_file = &errpaths_queuecomplete0_numpaths0_output_file;
	}
	if ((queue_broken) && (edited_error_portions.size() > 2))
	{
		errwrite_output_file = &errpaths_queuecomplete0_numpaths3plus_output_file;
	}
	outputfileMutex.lock();
	std::string original_error_block;
	found = original_error_portion.find("-");
	if (found!=std::string::npos)
	{
		original_error_block = before_first_error_kmer + after_last_error_kmer.substr(original_error_portion.size());
	}
	else
	{
		original_error_block = before_first_error_kmer + original_error_portion + after_last_error_kmer;
	}
	*errwrite_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
	*errwrite_output_file << before_first_error_kmer << " " << original_error_portion << " " << after_last_error_kmer << '\n';
	if (was_edited)
	{
		erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
		erredits_output_file << before_first_error_kmer << " " << original_error_portion << " " << after_last_error_kmer << '\n';
	}
	std::vector<uint32_t> w;
	//file.GetCountersForRead(original_error_block, w);
	if (shared_reads.size() > 30) {
		get_local_counts(original_error_block, local_kmer_db, w, k);
	} else {
		file.GetCountersForRead(original_error_block, w);
        }

	for (int j=0; j < w.size(); j++)
	{
		*errwrite_output_file << w.at(j) << " ";
		if (was_edited)
		{
			erredits_output_file << w.at(j) << " ";
		}
	}
	*errwrite_output_file << '\n';
	if (was_edited)
	{
		erredits_output_file << '\n';
	}
	for (int l = 0; l < edited_error_portions.size(); l++)
	{
		std::string edited_error_portion = edited_error_portions[l];
		std::string edited_error_block;
		found = edited_error_portion.find("-");
		if (found!=std::string::npos)
		{
			edited_error_block = before_first_error_kmer + after_last_error_kmer.substr(edited_error_portion.size());
		}
		else
		{
			edited_error_block = before_first_error_kmer + edited_error_portion + after_last_error_kmer;
		}
		*errwrite_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_path" << l << '\n';
		*errwrite_output_file << before_first_error_kmer << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
		if ((l==0) && (was_edited))
		{
			erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << '\n';
			erredits_output_file << before_first_error_kmer << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
		}
		//file.GetCountersForRead(edited_error_block, w);
		w.clear();

		if (shared_reads.size() > 30) {
			get_local_counts(edited_error_block, local_kmer_db, w, k);
		} else {
        	        file.GetCountersForRead(edited_error_block, w);
	        }


		for (int j=0; j < w.size(); j++)
		{
			*errwrite_output_file << w.at(j) << " ";
			if ((l==0) && (was_edited))
			{
				erredits_output_file << w.at(j) << " ";
			}
		}
		*errwrite_output_file << '\n';
		if ((l==0) && (was_edited))
		{
			erredits_output_file << '\n';
		}
	}
	outputfileMutex.unlock();
}

void write_nonhom_paths(CKMCFile& file, std::vector<int>& shared_reads, int k, bool& queue_broken, std::vector<std::string>& smoothed_nonhom_portions, std::ofstream& hetpaths_queuecomplete0_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& hetedits_output_file, std::string& smoothed_read, int& read_number, int& first_nonhom_idx, int& last_nonhom_idx, std::string& before_first_nonhom_kmer, std::string& original_nonhom_portion, std::string& after_last_nonhom_kmer, const unordered_map<auto, vector<string>>& local_kmer_db, int& strict, std::mutex& outputfileMutex)
{	
	std::string left_part;
	if (!before_first_nonhom_kmer.empty())
	{
		left_part = before_first_nonhom_kmer.substr(1);
	}
	else
	{
		left_part = before_first_nonhom_kmer;
	}
	std::string right_part;
	if (!after_last_nonhom_kmer.empty())
	{
		right_part = after_last_nonhom_kmer.substr(0, after_last_nonhom_kmer.length()-1);
	}
	else
	{
		right_part = after_last_nonhom_kmer;
	}
	auto is_greater_than = [&](std::string smoothed_nonhom_portion1, std::string smoothed_nonhom_portion2)
	{
		std::vector<uint32_t> w;
		std::size_t found1 = smoothed_nonhom_portion1.find("-");
		//added this case to account for when anchors overlap
		if (found1!=std::string::npos)
		{
			string left_right = left_part + right_part.substr(smoothed_nonhom_portion1.size());
			if (shared_reads.size() > 30) {
				get_local_counts(left_right, local_kmer_db, w, k);
			} else {
                        	file.GetCountersForRead(left_right, w);
                	}
			//file.GetCountersForRead(left_part + right_part.substr(smoothed_nonhom_portion1.size()), w);
		}
		else
		{
			string left_right = left_part + smoothed_nonhom_portion1 + right_part;
			if (shared_reads.size() > 30) {
				get_local_counts(left_right, local_kmer_db, w, k);
			} else {
                                file.GetCountersForRead(left_right, w);
                        }
			//file.GetCountersForRead(left_part + smoothed_nonhom_portion1 + right_part, w);
		}
		//std::sort(w.begin(), w.end());
		float average1 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
		std::size_t found2 = smoothed_nonhom_portion2.find("-");
		//added this case to account for when anchors overlap
		if (found2!=std::string::npos)
		{
			string left_right = left_part + right_part.substr(smoothed_nonhom_portion2.size());
			if (shared_reads.size() > 30) {
				get_local_counts(left_right, local_kmer_db, w, k);
			} else {
                                file.GetCountersForRead(left_right, w);
                        }
			//file.GetCountersForRead(left_part + right_part.substr(smoothed_nonhom_portion2.size()), w);
		}
		else
		{
			string left_right = left_part + smoothed_nonhom_portion2 + right_part;
			if (shared_reads.size() > 30) {
				get_local_counts(left_right, local_kmer_db, w, k);
			} else {
                                file.GetCountersForRead(left_right, w);
                        }
			//file.GetCountersForRead(left_part + smoothed_nonhom_portion2 + right_part, w);
		}
		//std::sort(w.begin(), w.end());
		float average2 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
		return average1 > average2;
	};
	std::sort(smoothed_nonhom_portions.begin(), smoothed_nonhom_portions.end(), is_greater_than);
	std::ofstream* hetwrite_output_file;
	//we finished the search and presumably we have found two heterozygous paths
	//actually, let's relax the !queue_broken requirement, jk restoring this requirement
	//and when strict==1 we only smooth when there are exactly two paths
	//and when strict==0 we smooth when there are at least two paths
	//if ((!queue_broken) && (smoothed_nonhom_portions.size() == 2))
	bool end_condition;
	bool was_smoothed = false;
	std::size_t found;
	if (strict==1)
	{
		end_condition = ((!queue_broken) && (smoothed_nonhom_portions.size() == 2));
	}
	else
	{
		end_condition = ((!queue_broken) && (smoothed_nonhom_portions.size() >= 2));
	}
	if (end_condition)
	{
		found = smoothed_nonhom_portions[0].find("-");
		if (found!=std::string::npos)
		{
			smoothed_read += left_part.substr(0, left_part.size() - smoothed_nonhom_portions[0].size());
		}
		else
		{
			smoothed_read += left_part + smoothed_nonhom_portions[0];
		}
		//We have checked that the condition for smoothing was met
		//Now we check whether the path chosen actually differs from the original
		//If so, then there was a true smooth and we add this to hetedits_output_file
		if (smoothed_nonhom_portions[0] != original_nonhom_portion)
		{
			was_smoothed = true;
		}
	}
	//there is not exactly two paths, or the search wasn't finished.
	//We are currently not smoothing.
	else
	{
		found = original_nonhom_portion.find("-");
		if (found!=std::string::npos)
		{
			smoothed_read += left_part.substr(0, left_part.size() - original_nonhom_portion.size());
		}
		else
		{
			smoothed_read += left_part + original_nonhom_portion;
		}
		
	}
	if ((!queue_broken) && (smoothed_nonhom_portions.size() == 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete1_numpaths2_output_file;
	}
	if ((!queue_broken) && (smoothed_nonhom_portions.size() < 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete1_numpaths0to1_output_file;
	}
	if ((!queue_broken) && (smoothed_nonhom_portions.size() > 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete1_numpaths3plus_output_file;
	}
	if ((queue_broken) && (smoothed_nonhom_portions.size() == 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete0_numpaths2_output_file;
	}
	if ((queue_broken) && (smoothed_nonhom_portions.size() < 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete0_numpaths0to1_output_file;
	}
	if ((queue_broken) && (smoothed_nonhom_portions.size() > 2))
	{
		hetwrite_output_file = &hetpaths_queuecomplete0_numpaths3plus_output_file;
	}
	outputfileMutex.lock();
	std::string original_nonhom_block;
	found = original_nonhom_portion.find("-");
	if (found!=std::string::npos)
	{
		original_nonhom_block = before_first_nonhom_kmer + after_last_nonhom_kmer.substr(original_nonhom_portion.size());
	}
	else
	{
		original_nonhom_block = before_first_nonhom_kmer + original_nonhom_portion + after_last_nonhom_kmer;
	}
	*hetwrite_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
	*hetwrite_output_file << before_first_nonhom_kmer << " " << original_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
	if (was_smoothed)
	{
		hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
		hetedits_output_file << before_first_nonhom_kmer << " " << original_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
	}
	std::vector<uint32_t> w;
	if (shared_reads.size() > 30) {
		get_local_counts(original_nonhom_block, local_kmer_db, w, k);
	} else {
		file.GetCountersForRead(original_nonhom_block, w);
        }
	//file.GetCountersForRead(original_nonhom_block, w);
	for (int j=0; j < w.size(); j++)
	{
		*hetwrite_output_file << w.at(j) << " ";
		if (was_smoothed)
		{
			hetedits_output_file << w.at(j) << " ";
		}
	}
	*hetwrite_output_file << '\n';
	if (was_smoothed)
	{
		hetedits_output_file << '\n';
	}
	for (int l = 0; l < smoothed_nonhom_portions.size(); l++)
	{
		std::string smoothed_nonhom_portion = smoothed_nonhom_portions[l];
		std::string smoothed_nonhom_block;
		found = smoothed_nonhom_portion.find("-");
		if (found!=std::string::npos)
		{
			smoothed_nonhom_block = before_first_nonhom_kmer + after_last_nonhom_kmer.substr(smoothed_nonhom_portion.size());
		}
		else
		{
			smoothed_nonhom_block = before_first_nonhom_kmer + smoothed_nonhom_portion + after_last_nonhom_kmer;
		}
		*hetwrite_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_path" << l << '\n';
		*hetwrite_output_file << before_first_nonhom_kmer << " " << smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
		if ((l==0) && (was_smoothed))
		{
			hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_smoothed" << '\n';
			hetedits_output_file << before_first_nonhom_kmer << " " << smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
		}
		w.clear();
		if (shared_reads.size() > 30) {
			get_local_counts(smoothed_nonhom_block, local_kmer_db, w, k);
		} else {
                	file.GetCountersForRead(smoothed_nonhom_block, w);
        	}
		//file.GetCountersForRead(smoothed_nonhom_block, w);
		for (int j=0; j < w.size(); j++)
		{
			*hetwrite_output_file << w.at(j) << " ";
			if ((l==0) && (was_smoothed))
			{
				hetedits_output_file << w.at(j) << " ";
			}
		}
		*hetwrite_output_file << '\n';
		if ((l==0) && (was_smoothed))
		{
			hetedits_output_file << '\n';
		}
	}
	outputfileMutex.unlock();
}

std::string remove_err (CKMCFile& file, std::vector<int>& shared_reads, std::vector<uint32_t>& v, std::string& read, int& read_number, const unordered_map<auto, vector<string>>& local_kmer_db, std::ofstream& errpaths_queuecomplete0_numpaths0_output_file, std::ofstream& errpaths_queuecomplete0_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& errpaths_queuecomplete1_numpaths0_output_file, std::ofstream& errpaths_queuecomplete1_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& erredits_output_file, int& error_threshold, int& het_threshold, int& unique_threshold, int& max_nodes_to_search, double& distance_multiplier, int& k, std::mutex& outputfileMutex)
{
	
	//initialize variables
	std::string edited_read;
	//int k = 21;

	//iterate over counts to edit errors
	int previous_type = -1;
	int first_nonerror_idx;
	int last_nonerror_idx;
	int first_error_idx;
	int last_error_idx;
	std::string before_first_error_kmer;
	std::string after_last_error_kmer;
	
	for (int i = 0; i < v.size(); i++)
	{
		
		int current_type = get_type(v[i], error_threshold, het_threshold, unique_threshold);
		//std::cout << "Current type: " << current_type << "\n";
		//std::cout.flush();
		
		//if kmer is error
		if (current_type == 0)
		{	
			
			//if this is the first kmer of the read
			if (previous_type == -1)
			{
				first_error_idx = i;
				last_nonerror_idx = i-1;
			}
			//if previous kmer was an error, we are continuing the error block
			
			if (previous_type == 0)
			{
				;
			}
			//if previous kmer was not an error, we are leaving the nonerror block
			if (previous_type > 0)
			{
				//get kmer that is right before the first error kmer of error block
				first_error_idx = i;
				last_nonerror_idx = i-1;
				before_first_error_kmer = read.substr(i-1, k);
				
				std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + 1);
				edited_read += nonerror_portion;
			}
			previous_type = current_type;
		}
		//if kmer is nonerror
		
		if (current_type > 0)
		{
			
			//if this is the first kmer of the read
			if (previous_type == -1)
			{
				first_nonerror_idx = i;
				last_error_idx = i-1;
			}
			
			//if previous kmer was error, and we are at the beginning of the read
			if (previous_type == 0 && before_first_error_kmer.empty())
			{
				
				//The very beginning of the read is an error portion
				after_last_error_kmer = read.substr(i, k);
				int min_distance_of_path = 0;
				int max_distance_of_path = i;
				bool queue_broken = false;
				std::vector<std::string> edited_error_portions = get_paths(file,shared_reads, local_kmer_db, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_error_idx = i-1;
				first_nonerror_idx = i;
				std::string original_error_portion = read.substr(0, i);
				
				write_error_paths(file,shared_reads, k, queue_broken, edited_error_portions, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, local_kmer_db, outputfileMutex);
			}
			
			//if previous kmer was error, we have left the error block
			
			if (previous_type == 0 && !before_first_error_kmer.empty())
			{
				
				int number_of_error_kmers = i - first_error_idx;
				//If the position of after_last_error_kmer overlaps before_first_error_kmer
				//we keep progressing as if nothing has happened, waiting to find another non_error kmer
				//if (number_of_error_kmers < k)
				//{
				//	current_type = previous_type;
				//	continue;
				//}
				//get kmer that is right after the last error kmer of block
				after_last_error_kmer = read.substr(i, k);
				//int min_distance_of_path = k;
				int min_distance_of_path = std::min(number_of_error_kmers, k);
				int max_distance_of_path = ceil(distance_multiplier * number_of_error_kmers);
				bool queue_broken;
				
				std::vector<std::string> edited_error_portions = get_paths(file,shared_reads, local_kmer_db, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				
				last_error_idx = i-1;
				first_nonerror_idx = i;
				std::string original_error_portion;
				//If the after_last_error_kmer overlaps before_first_error_kmer
				if (last_error_idx - first_error_idx + 2 - k < 0)
				{
					
					for (int number_overlaps=0; number_overlaps < first_error_idx - last_error_idx + k - 2; number_overlaps++)
					{
						original_error_portion += "-";
					}
				}
				else
				{
					original_error_portion = read.substr(first_error_idx+k-1, last_error_idx - first_error_idx + 2 - k);
				}
				
				write_error_paths(file,shared_reads, k, queue_broken, edited_error_portions, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, local_kmer_db, outputfileMutex);
			}
			//if previous kmer is nonerror, we are continuing a non error block
			if (previous_type > 0)
			{
				;
				
			}
			previous_type = current_type;
		}
	}
	
	//We have reached the end of the read, let's make sure we have added the last bit of the read
	if (previous_type == 0)
	{
		//std::cout << previous_type << "\n";
                //std::cout.flush();
		//We have "left" the error portion of the read
		after_last_error_kmer = "";
		int min_distance_of_path = 0;
		int max_distance_of_path = v.size()-first_error_idx;
		bool queue_broken = false;

		std::vector<std::string> edited_error_portions = get_paths(file,shared_reads, local_kmer_db, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
		
		last_error_idx = v.size()-1;
		first_nonerror_idx = v.size();
		
		std::string original_error_portion = read.substr(first_error_idx+k-1);
		write_error_paths(file,shared_reads, k, queue_broken, edited_error_portions, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer, local_kmer_db, outputfileMutex);
		//std::cout << previous_type << "\n";
                //std::cout.flush();
	}
	
	if (previous_type > 0)
	{
		//std::cout << previous_type << "\n";
		//std::cout.flush();
		//We have "left" the nonerror portion of the read
		first_error_idx = v.size();
		last_nonerror_idx = v.size()-1;
		std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + k);
		edited_read += nonerror_portion;
		
	}
	
	edited_read += '\n';
	return edited_read;
}

std::string smooth_het (CKMCFile& file, std::vector<int>& shared_reads, std::vector<uint32_t>& v, std::string& read, int& read_number, const unordered_map<auto, vector<string>>& local_kmer_db, std::ofstream& hetpaths_queuecomplete0_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& hetedits_output_file, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold, int& max_nodes_to_search, double& distance_multiplier, int& strict, int& k, std::mutex& outputfileMutex)
{
	//initialize variables
	std::string smoothed_read;
	//int k = 21;

	//iterate over counts to smoothe het
	int previous_type = -1;
	int first_hom_idx;
	int last_hom_idx;
	int first_nonhom_idx;
	int last_nonhom_idx;
	std::string before_first_nonhom_kmer;
	std::string after_last_nonhom_kmer;
	for (int i = 0; i < v.size(); i++)
	{
		
		std::string previous_kmer;
		int previous_count;
		if (i == 0)
		{
			previous_kmer = "";
			previous_count = 0;
		}
		else
		{
			previous_kmer = read.substr(i-1, k);
			std::vector<uint32_t> previous_count_vector;
			//file.GetCountersForRead(previous_kmer, previous_count_vector);
			if (shared_reads.size() > 30) {
				get_local_counts(previous_kmer, local_kmer_db, previous_count_vector, k);
			} else {
				file.GetCountersForRead(previous_kmer, previous_count_vector);
			}
			previous_count = previous_count_vector[0];
		}
		std::string current_kmer = read.substr(i, k);
		std::vector<uint32_t> current_count_vector;
		if (shared_reads.size() > 30) {
			get_local_counts(current_kmer, local_kmer_db, current_count_vector , k);
		} else {
			file.GetCountersForRead(current_kmer, current_count_vector);
		}
		//file.GetCountersForRead(current_kmer, current_count_vector);
		int current_count = current_count_vector[0];
		std::string anchor_found;
		int current_type = get_type_het(file,shared_reads, previous_type, previous_kmer, current_kmer, previous_count, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold, anchor_found);
		//if kmer is nonhom
		if ((current_type == 0) || (current_type == 1) || (current_type == 3))
		{
			//if this is the first kmer of the read
			if (previous_type == -1)
			{
				first_nonhom_idx = i;
				last_hom_idx = i-1;
			}
			//if previous kmer was nonhom, we are continuing the nonhom block
			if ((previous_type == 0) || (previous_type == 1) || (previous_type == 3))
			{
				;
			}
			//if previous kmer was hom, we are leaving the hom block
			if (previous_type == 2)
			{
				//get kmer that is right before the first nonhom kmer of nonhom block
				first_nonhom_idx = i;
				last_hom_idx = i-1;
				before_first_nonhom_kmer = read.substr(i-1, k);
				std::string hom_portion = read.substr(first_hom_idx, last_hom_idx - first_hom_idx + 1);
				smoothed_read += hom_portion;
				
			}
			previous_type = current_type;
		}
		//if kmer is hom
		if (current_type == 2)
		{
			//if this is the first kmer of the read
			if (previous_type == -1)
			{
				first_hom_idx = i;
				last_nonhom_idx = i-1;
			}
			//if previous kmer was nonhom, and we are at the beginning of the read
			if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && before_first_nonhom_kmer.empty())
			{
				//The very beginning of the read is an nonhom portion
				after_last_nonhom_kmer = read.substr(i, k);
				int min_distance_of_path = 0;
				int max_distance_of_path = i;
				bool queue_broken = false;
				std::vector<std::string> smoothed_nonhom_portions = get_paths(file,shared_reads, local_kmer_db, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_nonhom_idx = i-1;
				first_hom_idx = i;
				std::string original_nonhom_portion = read.substr(0, i);
				write_nonhom_paths(file,shared_reads, k, queue_broken, smoothed_nonhom_portions, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, local_kmer_db, strict, outputfileMutex);
			}
			//if previous kmer was nonhom, we have left the nonhom block
			if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && !before_first_nonhom_kmer.empty())
			{
				if (anchor_found == "left")
				{
					after_last_nonhom_kmer = read.substr(i, k);
					std::string hom_portion = read.substr(first_hom_idx+1, i-first_hom_idx-1);
					
					last_nonhom_idx = i-1;
					first_hom_idx = i;
					smoothed_read += hom_portion;
				}
				else
				{
				int number_of_nonhom_kmers = i - first_nonhom_idx;
				//If the position of after_last_nonhom_kmer overlaps before_first_nonhom_kmer
				//we keep progressing as if nothing has happened, waiting to find another hom kmer
				//if ((anchor_found == "right") && (number_of_nonhom_kmers < k))
				//{
				//	
				//	current_type = previous_type;
				//	continue;
				//}
				//get kmer that is right after the last nonhom kmer of block
				after_last_nonhom_kmer = read.substr(i, k);
				//int min_distance_of_path = k;
				int min_distance_of_path = std::min(number_of_nonhom_kmers, k);
				int max_distance_of_path = ceil(distance_multiplier * number_of_nonhom_kmers);
				bool queue_broken;
				std::vector<std::string> smoothed_nonhom_portions = get_paths(file,shared_reads, local_kmer_db, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
				last_nonhom_idx = i-1;
				first_hom_idx = i;
				std::string original_nonhom_portion;
				//If the after_last_nonhom_kmer overlaps before_first_nonhom_kmer
				if (last_nonhom_idx - first_nonhom_idx + 2 - k < 0)
				{
					for (int number_overlaps=0; number_overlaps < first_nonhom_idx - last_nonhom_idx + k - 2; number_overlaps++)
					{
						original_nonhom_portion += "-";
					}
				}
				else
				{
					original_nonhom_portion = read.substr(first_nonhom_idx+k-1, last_nonhom_idx - first_nonhom_idx + 2 - k);
				}
				write_nonhom_paths(file,shared_reads, k, queue_broken, smoothed_nonhom_portions, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, local_kmer_db, strict, outputfileMutex);
				}
			}
			//if previous kmer is hom, we are continuing a hom block
			if (previous_type == 2)
			{
				if (anchor_found == "left")
				{
					std::string hom_portion = read.substr(first_hom_idx, i-first_hom_idx);
					smoothed_read += hom_portion;
					
					first_hom_idx=i;
					last_nonhom_idx=i-1;
				}
			}
			previous_type = current_type;
		}
	}
	//We have reached the end of the read, let's make sure we have added the last bit of the read
	if ((previous_type == 0) || (previous_type == 1) || (previous_type == 3))
	{
		//We have "left" the nonhom portion of the read
		//If we have a homozygous on left with which to anchor
		if (first_nonhom_idx > 0)
		{
			after_last_nonhom_kmer = "";
			int min_distance_of_path = 0;
			int max_distance_of_path = v.size()-first_nonhom_idx;
			bool queue_broken = false;
			std::vector<std::string> smoothed_nonhom_portions = get_paths(file,shared_reads, local_kmer_db, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, k, queue_broken);
			last_nonhom_idx = v.size()-1;
			first_hom_idx = v.size();
			std::string original_nonhom_portion = read.substr(first_nonhom_idx+k-1);
			write_nonhom_paths(file,shared_reads, k, queue_broken, smoothed_nonhom_portions, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, local_kmer_db, strict, outputfileMutex);
		}
		else
		{
			smoothed_read += read;
		}
	}
	if (previous_type == 2)
	{
		//We have "left" the hom portion of the read
		first_nonhom_idx = v.size();
		last_hom_idx = v.size()-1;
		std::string hom_portion = read.substr(first_hom_idx, last_hom_idx - first_hom_idx + k);
		smoothed_read += hom_portion;
		
	}
	//smoothed_read += '\n';
	
	/*std::vector<uint32_t> vs;

        if (read.length() > 500) {
		get_local_counts(smoothed_read, local_kmer_db, vs, k, shared_reads);
        } else {
		file.GetCountersForRead(smoothed_read, vs);
        }
	std::cout << ">read" << read_number << "_post\n";
        for (uint32_t b: vs){
		std::cout << vs[b] << ",";
        }*/
        


	smoothed_read += '\n';
	return smoothed_read;
}

std::string getFileExt (const std::string &s)
{
	size_t i = s.find_last_of('.');
	if (i != std::string::npos)
	{
		return(s.substr(i+1, s.length() - i));
	}
	return("");
}

void processRead(int& cov_t, int& min_t, std::vector <std::vector <int>>& shared_reads_db, std::vector<std::string>& seqs, std::hash<std::string>& h, int& m, std::unordered_map <int, vector <int>>& min_db, std::ifstream& input_file, int& line_num, int& num_lines_per_read, CKMCFile& file, int& error_threshold, double& allowed_err_fraction, std::ofstream& err_output_file, std::ofstream& errpaths_queuecomplete0_numpaths0_output_file, std::ofstream& errpaths_queuecomplete0_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& errpaths_queuecomplete1_numpaths0_output_file, std::ofstream& errpaths_queuecomplete1_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& erredits_output_file, int& het_threshold, int& unique_threshold, int& max_nodes_to_search, double& distance_multiplier, int& k, int& polish, double& allowed_rep_fraction, std::ofstream& het_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& hetedits_output_file, int& anchor_threshold, int& strict, std::mutex& inputfileMutex, std::mutex& outputfileMutex, std::mutex& coutMutex, int& verbose)
{
	std::string line;
	int thread_line_num;
	std::string header;
	inputfileMutex.lock();
	//std::cout << "Processing the very first read\n" ;
	//std::cout.flush();
	while (getline(input_file, line))
	{
		
		line_num++;
		thread_line_num = line_num;
		//if on line with the header
		if (thread_line_num % num_lines_per_read == 0)
		{
			//save header, but make sure it is a fasta header, not a fastq header
			header = ">" + line.substr(1);
		}
		//if on line with the read
		if (thread_line_num % num_lines_per_read == 1)
		{
			
			inputfileMutex.unlock();
			std::string read = line;
			int read_number = (thread_line_num+num_lines_per_read-1)/num_lines_per_read - 1;
			//std::cout << "Read: " << read_number << "\n";
			//std::cout.flush();
			std::time_t thread_result = std::time(nullptr);
			if (verbose == 1)
			{
				coutMutex.lock();
				std::cout << "Analyzing read/contig/scaffold number  " << read_number << " at " << std::asctime(std::localtime(&thread_result)) << '\n';
				coutMutex.unlock();
			}
			//if (read_number%10000==0)
			//{
				//std::cout << read_number << '\n';
			//}
			
			//Find the shared reads TODO
			//This will be a list of reads that we are confident overlap
			//First step, use the minimizer DB and find reads that are similar 
			//std::vector<std::uint32_t> shared_reads;
			//int cov_t = 35000000;
			//std::vector<uint32_t> v_f; //Forward 
                        //file.GetCountersForRead(read, v_f);
			
			//std::vector<uint32_t> v_r; //Reverse
			//std::string read_rev = read;
		       	//reverse_c(read_rev);
			//file.GetCountersForRead(read_rev, v_r);
			//cout << "Got the OG read counts (F and R) \n";
			//compare_q(read_number[x], kmer_db, seqs, que_names, min_db, h, pos_names, k, m, cov_t, pos_mins, pos_overlaps, read_starts, read_stops);
			//std::cout << "About to get the shared reads\n";
                        //std::cout.flush();
			//find_shared_reads(min_t, shared_reads, seqs, v_f, v_r, k, h, m, cov_t, read_number, min_db, read, read_rev);
			//coutMutex.lock();
			//std::cout << "Found shared reads: " << shared_reads.size() << "\t" << read_number << " \n";	
			//coutMutex.unlock();
			//Make a mini Kmer DB based on the local reads TODO
			//map<string, uint32_t> local_kmer_db;
			//std::cout << "About to get shared reads\n";
			//std::cout.flush();
			auto start_time = chrono::high_resolution_clock::now();
			std::vector<int> shared_reads = shared_reads_db[read_number] ;
			/*if (shared_reads_db.count(read_number)){
				shared_reads =  shared_reads_db[read_number];
			} else {
				shared_reads = {0};
			}*/

			//std::cout << "Found shared reads\n";
			//std::cout.flush();
			auto end_time = chrono::high_resolution_clock::now();
        		auto total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
			coutMutex.lock();
        		std::cout << "FOUND SHARED READS " << total_time.count() << " MILLISECONDS" << '\n';
        		std::cout.flush();
			coutMutex.unlock();
			
			//std::vector<int> shared_reads = shared_reads_db[read_number];
			//std::cout << "Number of shared reads: " << shared_reads.size() << "\n";
			//std::cout.flush();
			//local_kmer_db.clear();
			std::unordered_map<int, vector<string>> local_kmer_db;
			//   unordered_map<string, vector<uint32_t>>& kmer_db
			std::cout << "About to make the local db\n";
			std::cout.flush();
			
			//make_local_kmer_db(shared_reads, seqs, local_kmer_db, k);
			//std::cout << "Made the local DB\n";
                        //std::cout.flush();
			
			//cout << "Made local KMER DB \n";
			//Only use the kmer DB to get counts for the rest of the smoothing TODO
			//Convert to the counts, not the reads they map with 
			std::vector<uint32_t> v;
			start_time = chrono::high_resolution_clock::now();
			if (read.length() > 2000) {
				if (shared_reads.size() > 30) {
					make_local_kmer_db(shared_reads, seqs, local_kmer_db, k);
					get_local_counts(read, local_kmer_db, v, k);
				} else {
					file.GetCountersForRead(read, v);
				}
			} else {

				file.GetCountersForRead(read, v);
			}
			end_time = chrono::high_resolution_clock::now();
                        total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
			coutMutex.lock();
                        std::cout << "FOUND LOCAL COUNTS " << total_time.count() << " MILLISECONDS" << '\n';
                        std::cout.flush();
			coutMutex.unlock();

			//std::cout << "Found local counts: " << read_number << "\n";
                        //std::cout.flush();

			//get_counts(read, local_kmer_db, v, k);
			//cout << "Got local counts \n";
			//TODO Blocked this out while running the DMEL tests 
			//std::vector<uint32_t> w;
                        //file.GetCountersForRead(read, w);
			//coutMutex.lock();
			//std::cout << ">read" << read_number << "_pre\n";
			//for (uint32_t b: v){
			//	std::cout << b << ",";
			//}
			//std::cout << "\n";
			//coutMutex.unlock();
			//cout << "V: " << v.size() << "\n";
			//cout << "W: " << w.size() << "\n";

			/*if (read_number == 0){
				std::vector<uint32_t> w;
				file.GetCountersForRead(read, w);

				for (int i = 0; i < v.size(); i++){
					cout << v[i] << "\t" << w[i] << "\n";
			}}*/
			
			//get counters of kmers in read
			//std::vector<uint32_t> v;
			//file.GetCountersForRead(read, v);
			int num_kmers = v.size();
			
			//check if too many error kmers
			//if greater than allowed_err_fraction of kmers in the read are errors
			//discard the read to not waste time
			/*int num_err = 0;
			for (int i = 0; i < v.size(); i++)
			{
				
				if (v[i] <= error_threshold)
				{
					num_err++;
				}
			}
			//cout << "Number of errors " << num_err << "\n" ;
			//std::cout << "Num err: " << num_err << "\n" ;
			//std::cout << "Went through the err threshold\n";
			//continue;
			double err_fraction = static_cast<double>(num_err)/num_kmers;
			if (err_fraction > allowed_err_fraction)
			{
				//std::cout << "read number " << read_number << " has " << err_fraction << " percent error kmers, discarding." << '\n';
				continue;
			}*/
			//std::cout << "Made sure that the err fraction is acceptable \n";
			//Calculate the lambda value and re-evaluate the thresholds 
			int error_threshold_local; // = error_threshold ;//round(ceil(0.25 * local_lambda)) ; //round(ceil(0.25 * l))
                        int het_threshold_local ;// = het_threshold;//round(ceil(1.5 * local_lambda));
                        int anchor_threshold_local; // = anchor_threshold;//round(ceil(2.5 * local_lambda));
                        int unique_threshold_local; // = unique_threshold;//round(ceil(3.5 * local_lambda));

			if (shared_reads.size() > 30) {
				if (read.length() > 2000) {
					float local_lambda = (shared_reads.size()/5.96) + 0.936;
					std::cout << "Local lambda: " << local_lambda << "\n" ;
					std::cout.flush();
					error_threshold_local = round(ceil(0.25 * local_lambda)) ; //round(ceil(0.25 * l))
					het_threshold_local = round(ceil(1.5 * local_lambda));
					anchor_threshold_local = round(ceil(2.5 * local_lambda));
					unique_threshold_local = round(ceil(3.5 * local_lambda));
			
				} else {
					error_threshold_local = error_threshold ;//round(ceil(0.25 * local_lambda)) ; //round(ceil(0.25 * l))
                                	het_threshold_local = het_threshold;//round(ceil(1.5 * local_lambda));
                                	anchor_threshold_local = anchor_threshold;//round(ceil(2.5 * local_lambda));
                                	unique_threshold_local = unique_threshold;
				}
			} else {
				error_threshold_local = error_threshold ;//round(ceil(0.25 * local_lambda)) ; //round(ceil(0.25 * l))
                                het_threshold_local = het_threshold;//round(ceil(1.5 * local_lambda));
                                anchor_threshold_local = anchor_threshold;//round(ceil(2.5 * local_lambda));
                                unique_threshold_local = unique_threshold;//round(ceil(3.5 * local_lambda));

			}
			//std::cout << "Found lambdas: " << read_number << "\n";
                        //std::cout.flush();
			//remove errors from the read to get edited read
			start_time = chrono::high_resolution_clock::now();
			std::string edited_read = remove_err(file,shared_reads, v, read, read_number, local_kmer_db, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, error_threshold_local, het_threshold_local, unique_threshold_local, max_nodes_to_search, distance_multiplier, k, outputfileMutex);
			//std::cout << time_me.lap() << " Removed Error \n";
			//cout << "Removed ERR \n";
			end_time = chrono::high_resolution_clock::now();
                        total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
			
			coutMutex.lock();
                        std::cout << "REMOVED ERROR " << total_time.count() << " MILLISECONDS" << '\n';
                        std::cout.flush();
			
			//coutMutex.lock();

			//std::cout << "Error Removed: " << read_number << "\n";
			//std::cout.flush();
			coutMutex.unlock();
			//write header and edited read to err_output_file
			outputfileMutex.lock();
			err_output_file << header << '\n';
			err_output_file << edited_read;
			outputfileMutex.unlock();
			edited_read.pop_back();
			if (polish == 1)
			{
				continue;
			}

			//get counters of kmers in edited read
			//file.GetCountersForRead(edited_read, v); TODO
			v.clear();
			if (shared_reads.size() > 30) {
				get_local_counts(edited_read, local_kmer_db, v, k);
			} else {
				file.GetCountersForRead(edited_read, v);
			}
			num_kmers = v.size();

			//check if too many repetitive kmers
			//if greater than allowed_rep_fraction of kmers in the edited read are repetitive
			//discard the read to not waste time
			int num_rep = 0;
			for (int i = 0; i < v.size(); i++)
			{
				//std::cout << v[i] << ' ';
				if (v[i] > unique_threshold)
				{
					num_rep++;
				}
			}
			//std::cout << "Num rep: " << num_rep << '\n';
			//continue;
			double rep_fraction = static_cast<double>(num_rep)/num_kmers;
			if (rep_fraction > allowed_rep_fraction)
			{
				//std::cout << "read number " << read_number << " has " << rep_fraction << " percent repetitive kmers, discarding." << '\n';
				outputfileMutex.lock();
				het_output_file << header << '\n';
				het_output_file << edited_read << '\n';
				outputfileMutex.unlock();
				continue;
			}
			//coutMutex.lock();
			//std::cout << "About to remove het: " << read_number << "\n";
			//coutMutex.unlock();
			start_time = chrono::high_resolution_clock::now();
			//smoothe het from the edited read to get smoothed read
			std::string smoothed_read = smooth_het(file,shared_reads, v, edited_read, read_number, local_kmer_db, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, error_threshold_local, het_threshold_local, unique_threshold_local, anchor_threshold_local, max_nodes_to_search, distance_multiplier, strict, k, outputfileMutex);
			
			end_time = chrono::high_resolution_clock::now();
                        total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

                        coutMutex.lock();
                        std::cout << "HET ERROR " << total_time.count() << " MILLISECONDS" << '\n';
                        std::cout.flush();
                        coutMutex.unlock();

			//cout << "Removed Het \n";
			//std::cout << time_me.lap() << " Het removed \n";
			/*std::cout << ">read" << read_number << "_post\n";
                        for (uint32_t b: v){
                                std::cout << v[b] << ",";
                        }
                        std::cout << "\n";*/
			//write header and smoothed read to het_output_file
			outputfileMutex.lock();
			het_output_file << header << '\n';
			het_output_file << smoothed_read;
			outputfileMutex.unlock();
			inputfileMutex.lock();

			std::vector<uint32_t> vs;

        		/*if (read.length() > 2000) {
				if (shared_reads.size() > 30) {
                			get_local_counts(smoothed_read, local_kmer_db, vs, k, shared_reads);
        			} else {
					file.GetCountersForRead(smoothed_read, vs);
				}
			} else {
                		file.GetCountersForRead(smoothed_read, vs);
        		}*/
			//coutMutex.lock(); TODO
			//blocked out this while running the DMEL test. 
        		//std::cout << ">read" << read_number << "_post\n";
        		//for (uint32_t b: vs){
                	//	std::cout << b << ",";
        		//}
        		//std::cout << "\n";
			//coutMutex.unlock();

		}
	}
	inputfileMutex.unlock();
}

void read_fasta_slow(std::ifstream& fasta_file_que, std::vector<string>& seqs){
	int seq_cntr; seq_cntr=0;
        int name_cntr; name_cntr=0;
        
	std::ifstream &file(fasta_file_que);
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
                                name_cntr++;
			}
                }
                //file.close();
        }
	return void();
}

int main(int argc, char* argv[])
{
	//parse arguments
	
	int c;
	std::ifstream input_file;
	int num_lines_per_read;
	std::string kmcdb;
	std::string outdir;
	std::ofstream err_output_file;
	std::ofstream errpaths_queuecomplete0_numpaths0_output_file;
	std::ofstream errpaths_queuecomplete0_numpaths1to2_output_file;
	std::ofstream errpaths_queuecomplete0_numpaths3plus_output_file;
	std::ofstream errpaths_queuecomplete1_numpaths0_output_file;
	std::ofstream errpaths_queuecomplete1_numpaths1to2_output_file;
	std::ofstream errpaths_queuecomplete1_numpaths3plus_output_file;
	std::ofstream erredits_output_file;
	std::ofstream het_output_file;
	std::ofstream hetpaths_queuecomplete0_numpaths0to1_output_file;
	std::ofstream hetpaths_queuecomplete0_numpaths2_output_file;
	std::ofstream hetpaths_queuecomplete0_numpaths3plus_output_file;
	std::ofstream hetpaths_queuecomplete1_numpaths0to1_output_file;
	std::ofstream hetpaths_queuecomplete1_numpaths2_output_file;
	std::ofstream hetpaths_queuecomplete1_numpaths3plus_output_file;
	std::ofstream hetedits_output_file;
	
	int k = 0;
	double l = 0;
	int error_threshold = 0;
	int het_threshold = 0;
	int unique_threshold = 0;
	int anchor_threshold = 0;
	double allowed_err_fraction = 1;
	double allowed_rep_fraction = 1;
	int max_nodes_to_search = 1000;
	double distance_multiplier = 1.2;
	int strict = 1;
	int polish = 0;
	int num_threads = 1;
	int verbose = 0;
	int min_t = -1;
	int cov_t = 35;
	int m = 21;
	//std::cout << "Have not made it very far at all\n" ;	
	while ((c = getopt(argc, argv, "hi:j:o:k:M:C:z:l:g:m:a:u:t:pn:d:e:r:s:v")) != -1)
	{
		switch (c)
		{
			case 'h':
				fprintf(stderr, "Usage: %s -i input.fa/fq -j kmcdb -o outdir -k kmersize -l lambda [-t num_threads (default 1)] [-p (run in polish mode, i.e. run only error correction and not het smoothing)] [-n max_nodes_to_search (default 1000)] [-d distance_multiplier (default 1.2)] [-e allowed_err_fraction (default 1.0)] [-r allowed_rep_fraction (default 1.0)] [-s strict (0 or 1, default 1)] [-v (run in verbose mode, i.e. print time stamps for analyses)]\n", argv[0]);
				exit(EXIT_FAILURE);
			case 'i':
				//is the input fasta or fastq?
				//note: the output will be fasta format, since quality values
				//will not match once the read is edited and smoothed
				if ((getFileExt(optarg) == "fasta") || (getFileExt(optarg) == "fa"))
				{
					num_lines_per_read = 2;
				}
				else if ((getFileExt(optarg) == "fastq") || (getFileExt(optarg) == "fq"))
				{
					num_lines_per_read = 4;
				}
				else
				{
					fprintf(stderr, "Input filename must end in .fa .fasta .fq or .fastq.\n");
					exit(EXIT_FAILURE);
				}
				input_file.open(optarg);
				if (!input_file.is_open())
				{
					fprintf(stderr, "Please ensure %s exists.\n", optarg);
					exit(EXIT_FAILURE);
				}
				break;
			case 'j':
				kmcdb = optarg;
				break;
			case 'o':
				outdir = optarg;
				err_output_file.open(outdir+"/errremoved.fasta");
				errpaths_queuecomplete0_numpaths0_output_file.open(outdir+"/errpaths_queuecomplete0_numpaths0.fasta");
				errpaths_queuecomplete0_numpaths1to2_output_file.open(outdir+"/errpaths_queuecomplete0_numpaths1to2.fasta");
				errpaths_queuecomplete0_numpaths3plus_output_file.open(outdir+"/errpaths_queuecomplete0_numpaths3plus.fasta");
				errpaths_queuecomplete1_numpaths0_output_file.open(outdir+"/errpaths_queuecomplete1_numpaths0.fasta");
				errpaths_queuecomplete1_numpaths1to2_output_file.open(outdir+"/errpaths_queuecomplete1_numpaths1to2.fasta");
				errpaths_queuecomplete1_numpaths3plus_output_file.open(outdir+"/errpaths_queuecomplete1_numpaths3plus.fasta");
				erredits_output_file.open(outdir+"/erredits.fasta");
				het_output_file.open(outdir+"/hetremoved.fasta");
				hetpaths_queuecomplete0_numpaths0to1_output_file.open(outdir+"/hetpaths_queuecomplete0_numpaths0to1.fasta");
				hetpaths_queuecomplete0_numpaths2_output_file.open(outdir+"/hetpaths_queuecomplete0_numpaths2.fasta");
				hetpaths_queuecomplete0_numpaths3plus_output_file.open(outdir+"/hetpaths_queuecomplete0_numpaths3plus.fasta");
				hetpaths_queuecomplete1_numpaths0to1_output_file.open(outdir+"/hetpaths_queuecomplete1_numpaths0to1.fasta");
				hetpaths_queuecomplete1_numpaths2_output_file.open(outdir+"/hetpaths_queuecomplete1_numpaths2.fasta");
				hetpaths_queuecomplete1_numpaths3plus_output_file.open(outdir+"/hetpaths_queuecomplete1_numpaths3plus.fasta");
				hetedits_output_file.open(outdir+"/hetedits.fasta");
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'M':
				//similarity threshold (min_t)
                                min_t = atoi(optarg);
                                break;
			case 'C':
                                //similarity threshold (min_t)
                                cov_t = atoi(optarg);
                                break;
			case 'z':
                                //similarity threshold (min_t)
                                m = atoi(optarg);
                                break;
			case 'l':
				//only set values based on lambda if not already set by hidden parameters g, m, a, or u
				l = std::stod(optarg);
				if (error_threshold == 0)
				{
					error_threshold = round(ceil(0.25 * l));
				}
				if (het_threshold == 0)
				{
					het_threshold = round(ceil(1.5 * l));
				}
				if (unique_threshold == 0)
				{
					unique_threshold = round(ceil(3.5 * l));
				}
				if (anchor_threshold == 0)
				{
					anchor_threshold = round(ceil(2.5 * l));
				}
				break;
			case 'g':
				error_threshold = atoi(optarg);
				break;
			case 'm':
				het_threshold = atoi(optarg);
				break;
			case 'a':
				anchor_threshold = atoi(optarg);
				break;
			case 'u':
				unique_threshold = atoi(optarg);
				break;
			case 't':
				num_threads = atoi(optarg);
				break;
			case 'p':
				polish = 1;
				break;
			case 'n':
				max_nodes_to_search = atoi(optarg);
				break;
			case 'd':
				distance_multiplier = std::stod(optarg);
				break;
			case 'e':
				allowed_err_fraction = std::stod(optarg);
				break;
			case 'r':
				allowed_rep_fraction = std::stod(optarg);
				break;
			case 's':
				strict = atoi(optarg);
				break;
			case 'v':
				verbose = 1;
				break;
			case '?':
				fprintf(stderr, "Option -%c is invalid or requires an argument.\n", optopt);
			default:
				fprintf(stderr, "Usage: %s -i input.fa/fq -j kmcdb -o outdir -k kmersize -l lambda [-t num_threads (default 1)] [-p (run in polish mode, i.e. run only error correction and not het smoothing)] [-n max_nodes_to_search (default 1000)] [-d distance_multiplier (default 1.2)] [-e allowed_err_fraction (default 1.0)] [-r allowed_rep_fraction (default 1.0)] [-s strict (0 or 1, default 1)] [-v (run in verbose mode, i.e. print time stamps for analyses)]\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}
	//std::cout << "Checking the arrguments\n";
	//check that required arguments are given
	if (!input_file.is_open())
	{
		fprintf(stderr, "Please provide input file with -i argument.\n");
		exit(EXIT_FAILURE);
	}
	if (kmcdb.empty())
	{
		fprintf(stderr, "Please provide kmcdb with -j argument.\n");
		exit(EXIT_FAILURE);
	}
	if (outdir.empty())
	{
		fprintf(stderr, "Please provide output directory with -o argument.\n");
		exit(EXIT_FAILURE);
	}
	if (k==0)
	{
		fprintf(stderr, "Please provide kmer size with -k argument.\n");
		exit(EXIT_FAILURE);
	}
	if ((error_threshold == 0) || (het_threshold == 0) || (anchor_threshold == 0) || (unique_threshold == 0))
	{
		fprintf(stderr, "Please provide average kmer coverage with -l argument. (Or specify the thresholds manually).\n");
		exit(EXIT_FAILURE);
	}
	//std::cout << error_threshold << " " << het_threshold << " " << anchor_threshold << " " << unique_threshold << "\n";
	//load KMC database
	CKMCFile file;
	std::time_t result = std::time(nullptr);
	//std::cout << "Loading KMC database at " << std::asctime(std::localtime(&result)) << '\n';
	auto start_time = chrono::high_resolution_clock::now();
	file.OpenForRA(kmcdb);
	result = std::time(nullptr);
	auto end_time = chrono::high_resolution_clock::now();
	auto total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
	
	std::cout << "KMC database took " << total_time.count() << " MILLISECONDS" << '\n';
	std::cout.flush();	
	//Now we need to make the minimizer DB. This can be a global variable?  TODO
	//First we will have to go sequence by sequence, find the kmers in each sequence,
	//use the getcounts function to get the counts, then send this off to the minimizer 
	//db maker. 
	std::vector<std::string> seqs;
	//Time how long it takes to read in the file
	start_time = chrono::high_resolution_clock::now();
	read_fasta_slow(input_file, seqs);
	end_time = chrono::high_resolution_clock::now();
	total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
	std::cout << "Read the file " << total_time.count() << " MILLISECONDS\n";
	std::cout.flush();

	input_file.clear();
	input_file.seekg(0, std::ios::beg);
	//if (!input_file.is_open()) {
	//	cout << "OH NO!\n";
	//}
	std::hash<std::string> h;
	int seq_cntr=seqs.size();
	
	//int cov_t=40;
	//min_db is the database where we keep track of what minimizers we have found. 
	//key = minimizer
	//value = vector of read ids 
	std::unordered_map <int, vector <int>> min_db;
	
	//all_v_all is an unordered map where we keep track of the number of times 
	//different reads shared minimizers 
	//key = read id
	//value = unordered map of read ids 
	//	key = read id 
	//	value = number of times they overlap 
	//std::unordered_map <int, std::unordered_map<int, uint32_t>> all_v_all;
	std::unordered_map <float, uint32_t> all_v_all;

	//And the shared_reads_db keeps track of how many reads are shared 
	//key = read id 
	//value = vector of other read ids 
	//std::unordered_map<int, vector<int>> shared_reads_db;
	std::vector <std::vector <int>> shared_reads_db;
	//std::unordered_map<int, vector<int>>
        //TimerSK time_me;
	//cout << "Set up for minimizer db done" << '\n';
	//std::cout << time_me.lap() << " Set up for the minimizer DB is done\n";	
	//Now it's time to find the minimizers
	//std::cout << "About to make these three weird datasets\n";
	//std::cout.flush();
	start_time = chrono::high_resolution_clock::now();
	for (int x=0; x < seq_cntr; x++){
		//Maybe we should just add this right away? 
		shared_reads_db[x] = {x};
		//For each read we will need to find the minimizers, which will require the counts 
                std::vector<uint32_t> cnts_q;
                file.GetCountersForRead(seqs[x], cnts_q);
		//get_counts(seqs[x], kmer_db, cnts_q);
                //get the minimizers for that seq
		std::vector <int> line_mins;
		//std::cout << "About to find mins\n" ;
		//std::cout.flush();
		find_mins(seqs[x], cnts_q, line_mins, k, m, cov_t);
		//Now we want to go through the minimizers and add these the the minimizer db one at a time 
		//if (line_mins.empty()) {
		//	std::cout << "EMTPY\n" ;
		//	std::cout.flush();
		//}
		//std::cout << "Found minimizers: " << x << "\t" << line_mins.size() << "\n";
        	//std::cout.flush();
		//std::cout << x << "\n";
		//std::cout.flush();	
		for (int im=0; im < line_mins.size(); im++){	
			int this_m = line_mins[im];

			if (min_db.count(this_m)){
				//This minimizer has already been found. That's cool though. 
				//We'll need to look at every other value already in the vector 
				vector<int> tmp_vec = min_db[this_m];
				for (int y=0; y < tmp_vec.size(); y++) {
					//Alright, we now need to go through this vector (containing read ids) that share this minimizer 
					//First, check if this read id (tmp_vec[y]) is in the all_v_all keys 
					int tmp_rd_id = tmp_vec[y]; // This is the other read id that we are interested in. 
					//if (x >= tmp_rd_id){
						//Since tmp_rd_id is smaller, we are going to add the pointer to x in this read id. 
						//Make the float id number with the two rd numbers 
						std::string the_two_reads = std::to_string(tmp_rd_id) + "." + std::to_string(x) ;
						float this_is_the_id = std::stof(the_two_reads);
						std::cout << this_is_the_id << "\n";
						std::cout.flush();
						//std::unordered_map<int, uint32_t> this_map;
						//Check to see if the temp read is already in the all_v_all dataset 
						//All_v_all will keep track of how many times two reads overlap
						//Maybe concatanating the read numbers will be a better way of keeping track. Rather than having a map of a map? 
						//Or making them into a float? With the dot seperating them? Let's test this first. We can time it a few tmes. 
						if (all_v_all.count(this_is_the_id)) {
							/*this_map = all_v_all[tmp_rd_id];
							if (this_map.count(x)){
								++this_map[x];
							} else {
								this_map[x] = 1;
							}*/
							++all_v_all[this_is_the_id];
						} else {
							//std::unordered_map<int, uint32_t> this_map;
							//this_map[x] = 1;
							all_v_all[this_is_the_id] = 1;
						}
						//all_v_all[x] ;
						if (all_v_all[this_is_the_id] == m ){	
							//Now we have to update the shared read database, for both of these reads 
							//Read 1
							std::cout << x << "\t" << shared_reads_db.size() << "\n";
							std::cout.flush();
							//if (shared_reads_db.size() == (x+1)){ //(shared_reads_db.count(x)) {
							vector<int> another_tmp = shared_reads_db[x];
							if (another_tmp.back() != tmp_rd_id) {
								another_tmp.push_back(tmp_rd_id);
								shared_reads_db[x] = another_tmp;
							}
								/*	}
							} else {
								shared_reads_db[x] = {x, tmp_rd_id};
							}*/
							//Check if this read is already in the database 
							//Read 2
							//if (shared_reads_db.count(tmp_rd_id)) {
                                                        //vector<int> another_tmp = shared_reads_db[tmp_rd_id];
                                                        //if (another_tmp.back() != x){
							//	another_tmp.push_back(x);
                                                        //        shared_reads_db[tmp_rd_id] = another_tmp;
							//}
                                                        //} else {
                                                        //        shared_reads_db[tmp_rd_id] = {tmp_rd_id, x};
                                                        //}
						
						}
					/*} else {
						//std::unordered_map<int, uint32_t> this_map;
                                                std::string the_two_reads = std::to_string(x) + "." + std::to_string(tmp_rd_id) ;
						float this_is_the_id = std::stof(the_two_reads);
						if (all_v_all.count(this_is_the_id)) {
                                                        this_map = all_v_all[x];
                                                        if (this_map.count(tmp_rd_id)){
                                                                ++this_map[tmp_rd_id];
                                                        } else {
                                                                this_map[tmp_rd_id] = 1;
                                                        }
                                                        ++all_v_all[this_is_the_id] ; //= this_map;
                                                } else {
                                                        //std::unordered_map<int, uint32_t> this_map;
                                                        //this_map[tmp_rd_id] = 1;
                                                        all_v_all[this_is_the_id] = 1;
                                                }
						//++all_v_all[x];
						if (all_v_all[this_is_the_id] == m){
							//Now we update the shared read database for both reads, and we have to make sure the reads are in the db. 
							if (shared_reads_db.count(x)) {
                                                                vector<int> another_tmp = shared_reads_db[x];
                                                                if (another_tmp.back() != tmp_rd_id) {
                                                                       another_tmp.push_back(tmp_rd_id);
                                                                        shared_reads_db[x] = another_tmp;
                                                                }
                                                        } else {
                                                                shared_reads_db[x] = {x, tmp_rd_id};
                                                        }
                                                        //Check if this read is already in the database 
                                                        //Read 2
                                                        if (shared_reads_db.count(tmp_rd_id)) {
                                                                vector<int> another_tmp = shared_reads_db[tmp_rd_id];
                                                                if (another_tmp.back() != x){
                                                                        another_tmp.push_back(x);
                                                                        shared_reads_db[tmp_rd_id] = another_tmp;
                                                                }
                                                        } else {
                                                                shared_reads_db[tmp_rd_id] = {tmp_rd_id, x};
                                                        }
						} 

					}*/
				}
				/*if (shared_reads_db.size() == x){
						shared_reads_db[x] = {x};
				}*/
					//If yes, then we still need to check if our read is in the all_v_all key. 
					//If it isn't, check if our read id is in any of the all_v_all keys 
					//If neither is, then we'll add our read id to the all_v_all keys, and make the value a map
					//The map will contain the tmp_vec[y] as the key and the value will be one. 
				//}
				
				//And update the all vs all records to include this read 
				//And if the all vs all records are already above the threshold, then update the shared_reads_db 
				//Before we leave, we need to update the min_db 
				tmp_vec.push_back(x);
				min_db[this_m] = tmp_vec;
			} else {
				//Add this minimizer to the min db, which means adding this read id as well 
				//This is an easy case! 
				min_db[this_m] = {x};
			}
		}
        }
	end_time = chrono::high_resolution_clock::now();
        total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

        std::cout << "Making the three databases took " << total_time.count() << " MILLISECONDS" << '\n';
        std::cout.flush();

	//std::cout << time_me.lap() << " Made the min DB \n";
	//cout << "Made min DB" << std::asctime(std::localtime(&result)) << '\n';
	//global variables
	std::cout << "Made these three weird datasets\n";
        std::cout.flush();
	int line_num = -1;
	//std::cout << "About to make the local db";
	//std::cout.flush();
	//std::unordered_map<string, vector<uint32_t>> l_kmer_db;
	//make_kmer_db(seqs, local_kmer_db, k);

	std::cout << "Made the min DB \n" ;
	std::cout.flush();
	//std::cout << time_me.lap() << " Made the kmer DB \n";
	//cout << "Made KMER DB" << std::asctime(std::localtime(&result)) << '\n';
	
	//make locks for the threads
	std::mutex inputfileMutex;
	std::mutex outputfileMutex;
	std::mutex coutMutex;
	
	//process reads using threads
	std::vector<std::thread> threads;
	start_time = chrono::high_resolution_clock::now();
	for (int i = 0; i < num_threads; i++)
	{
		std::cout << i << "\n";
		std::cout.flush();
		threads.push_back(std::thread(processRead, std::ref(cov_t), std::ref(min_t), std::ref(shared_reads_db), std::ref(seqs), std::ref(h), std::ref(m), std::ref(min_db), std::ref(input_file), std::ref(line_num), std::ref(num_lines_per_read), std::ref(file), std::ref(error_threshold), std::ref(allowed_err_fraction), std::ref(err_output_file), std::ref(errpaths_queuecomplete0_numpaths0_output_file), std::ref(errpaths_queuecomplete0_numpaths1to2_output_file), std::ref(errpaths_queuecomplete0_numpaths3plus_output_file), std::ref(errpaths_queuecomplete1_numpaths0_output_file), std::ref(errpaths_queuecomplete1_numpaths1to2_output_file), std::ref(errpaths_queuecomplete1_numpaths3plus_output_file), std::ref(erredits_output_file), std::ref(het_threshold), std::ref(unique_threshold), std::ref(max_nodes_to_search), std::ref(distance_multiplier), std::ref(k), std::ref(polish), std::ref(allowed_rep_fraction), std::ref(het_output_file), std::ref(hetpaths_queuecomplete0_numpaths0to1_output_file), std::ref(hetpaths_queuecomplete0_numpaths2_output_file), std::ref(hetpaths_queuecomplete0_numpaths3plus_output_file), std::ref(hetpaths_queuecomplete1_numpaths0to1_output_file), std::ref(hetpaths_queuecomplete1_numpaths2_output_file), std::ref(hetpaths_queuecomplete1_numpaths3plus_output_file), std::ref(hetedits_output_file), std::ref(anchor_threshold), std::ref(strict), std::ref(inputfileMutex), std::ref(outputfileMutex), std::ref(coutMutex), std::ref(verbose)));
	}
	end_time = chrono::high_resolution_clock::now();
        total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

        std::cout << "Processing the reads took " << total_time.count() << " MILLISECONDS" << '\n';
        std::cout.flush();

	for (auto &th : threads)
	{
		th.join();
	}
	
	//close files
	input_file.close();
	err_output_file.close();
	errpaths_queuecomplete0_numpaths0_output_file.close();
	errpaths_queuecomplete0_numpaths1to2_output_file.close();
	errpaths_queuecomplete0_numpaths3plus_output_file.close();
	errpaths_queuecomplete1_numpaths0_output_file.close();
	errpaths_queuecomplete1_numpaths1to2_output_file.close();
	errpaths_queuecomplete1_numpaths3plus_output_file.close();
	erredits_output_file.close();
	het_output_file.close();
	hetpaths_queuecomplete0_numpaths0to1_output_file.close();
	hetpaths_queuecomplete0_numpaths2_output_file.close();
	hetpaths_queuecomplete0_numpaths3plus_output_file.close();
	hetpaths_queuecomplete1_numpaths0to1_output_file.close();
	hetpaths_queuecomplete1_numpaths2_output_file.close();
	hetpaths_queuecomplete1_numpaths3plus_output_file.close();
	hetedits_output_file.close();
	return 0;
}
