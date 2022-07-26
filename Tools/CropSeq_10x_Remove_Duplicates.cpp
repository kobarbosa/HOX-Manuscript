/*
 * 210126_10x_dedup_reads_for_Karina_CROPseq.cpp
 *
 *  Created on: Jan 26, 2021
 *      Author: anna
 */

/* MODIFIED SCRIPT:
 * 190223_sciRNA_remove_duplicates.cpp
 *
 *  Created on: Feb 23, 2019
 *      Author: anna
 *
 *      This script takes in a .sam file (processed) and outputs a .sam file w/ duplicates removed.
 *      It also counts # of duplicates corresponding to each index. It also prints out a .sam file of
 *      sequences not corresponding to any index.
 *
 *      Here, I saved indexes in csv format, and substring the length of the index -- there are hidden characters that are problematic.
 *
 *
 *index_list.open(argv[1]);
 *input_sam_file.open(argv[2]);
 *output_sam_file.open(argv[3]);
 *output_sam_tossed_seq_file.open(argv[4]);
 *reads_and_dup_metrics_file.open(argv[5]);
 *bad_index_counts_file.open(argv[6]);
 *
 *
 *Pseudocode:
 *
 *Input indexes into vector

For each line of sortedinput file:
  Get genomic position
  If RT index doesn't match any in map
    Check hamming distance w/ all other indexes
      if find one, change that index to the correct sequence
  If not,
     Record tossed sequence information

  If real index and Tn5 integration not within 2bp of last sequence:
     write into outfile
     *clear recent vector
     *put into vector
  If real index but integration within 2bp of last sequence:
     If UMI & Index also match, then add to total reads for that index (but not unique reads; don't put sequence into outfile.
     If different UMI & Index
       **check if you've already seen it!
       if not
         put into outrile
         put into vector
       if so
         don't put into outfile (count as seen)
 *
 *
 *
 *
 *
 *
 */


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <ostream>
#include <stdio.h>
#include <string.h>

using namespace std;

struct index_struct{
	string index_seq;
	int num_unique_reads;
	int num_total_reads;
};

void print_index_struct(index_struct is1){
	cout << is1.index_seq << endl;
}

int same_base_counter;
bool check_hamming_distance(string str1, string str2, int max_hamming_dist)
{
	same_base_counter = 0;
	for (int i = 0; i < str1.size(); i++)
	{
		if(str1.at(i) == str2.at(i))
		{
			same_base_counter++;
		}
	}

	if(str1.length() - same_base_counter <= max_hamming_dist)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool check_if_index_in_map(map<string, int> &bad_index_map, string current_index)
{
	if (bad_index_map.find(current_index) != bad_index_map.end())
	{
		return true;
	}
	else
	{
		return false;
	}

}

bool found_matching_sequence;
/* ***THIS WAS THE OLD FUNCTION!!!
bool check_if_string_in_vector(vector<string>umi_index_vector, string curr_umi_index)
{
	found_matching_sequence = false;
	for (int i = 0; i < umi_index_vector.size(); i++)
	{
		if (umi_index_vector[i] == curr_umi_index)
		{
			found_matching_sequence = true;
		}
	}
	return found_matching_sequence;
}
*/


bool check_if_string_in_vector(vector<string>umi_index_vector, string curr_umi_index)
{
	found_matching_sequence = false;
	for (int i = 0; i < umi_index_vector.size(); i++)
	{
		if (check_hamming_distance(umi_index_vector[i], curr_umi_index, 2) == true)
		{
			found_matching_sequence = true;
		}
		//if (umi_index_vector[i] == curr_umi_index)
		//{
		//	found_matching_sequence = true;
		//}
	}
	return found_matching_sequence;
}

vector<string> get_all_1_bp_mismatches(string gRNA_seq){
	vector<string> all_1_bp_mismatch_vect;
	string bases = "ATCGN";
	string temp_base;
	string temp_replacement_base;
	string new_gRNA_seq;
	for(int i = 0; i < gRNA_seq.size(); i++){
		//new_gRNA_seq = gRNA_seq;
		temp_base = gRNA_seq[i];
		for(int j = 0; j < bases.size(); j++){
			temp_replacement_base = bases[j];
			if(temp_base != temp_replacement_base){
				new_gRNA_seq = gRNA_seq;
				new_gRNA_seq.replace(i,1,temp_replacement_base);
				all_1_bp_mismatch_vect.push_back(new_gRNA_seq);
			}
		}
	}
	return(all_1_bp_mismatch_vect);
}

int main(int argc, char* argv[]){

	vector<index_struct> vector_of_indexes;
	map<string, int> index_to_pos_in_vect_map;
	map<string,string> mismatch_index_to_real_index_map;
	vector<string> temp_mismatch_vector;
	map<string, int> bad_index_map;
	index_struct temp_index_struct;

	ifstream(index_list);
	index_list.open(argv[1]);
	string index_list_line;

	while(getline(index_list, index_list_line)) //make vector if index structs
	{
		cout << index_list_line << endl;
		temp_index_struct.index_seq = index_list_line.substr(0,16);
		cout << temp_index_struct.index_seq << endl;
		temp_index_struct.num_total_reads = 0;
		temp_index_struct.num_unique_reads = 0;
		vector_of_indexes.push_back(temp_index_struct);
		index_to_pos_in_vect_map[index_list_line.substr(0,16)] = vector_of_indexes.size() - 1;
		temp_mismatch_vector = get_all_1_bp_mismatches(temp_index_struct.index_seq);
		for(int i = 0; i < temp_mismatch_vector.size(); i++){
			mismatch_index_to_real_index_map[temp_mismatch_vector[i]] = temp_index_struct.index_seq;
		}
		temp_index_struct = index_struct();
	}
	cout << "Num Cells: " << index_to_pos_in_vect_map.size() << endl;
	cout << "Num Cells w/ 1 Mismatch: " << mismatch_index_to_real_index_map.size() << endl;

	int last_genomic_pos = 0;
	string last_umi = "AAAAAAAA";
	string last_index = "AAAAAAAAAA";

	int curr_genomic_pos = 0;
	string curr_umi = "AAAAAAAA";
	string curr_index = "AAAAAAAAAA";


	ifstream(input_sam_file);
	input_sam_file.open(argv[2]);
	string input_sam_line;
	string sam_word_in_line;
	stringstream input_sam_line_ss;
	int word_counter = 0;

	bool ok_hamming_distance;
	int vector_pos_counter;
	//bool should_correct_index = false;

	ofstream(output_sam_file);
	output_sam_file.open(argv[3]);
	ofstream(output_sam_tossed_seq_file);
	output_sam_tossed_seq_file.open(argv[4]);

	vector<string> recent_umi_index_combos_vect;
	bool already_in_vector = false;

	map<string, string> recent_umi_index_combos_vect_MAP;

	int test_counter = 0;

	cout << "print this at least" << endl;

	int sam_file_line_counter = 0;

	while(getline(input_sam_file, input_sam_line))
	{
		sam_file_line_counter++;
		if(sam_file_line_counter%10000 == 0){
			cout << sam_file_line_counter << endl;
		}
		//cout << "in while loop!!" << endl;
		ok_hamming_distance = true;
		curr_umi = input_sam_line.substr(input_sam_line.find("_") + 1, 12);
		if (test_counter < 10)
		{
			cout << curr_umi << endl;
			test_counter++;
		}
		//curr_index = input_sam_line.substr(input_sam_line.find("_") + 10, 10);
		curr_index = input_sam_line.substr(0, 16);
		if (test_counter < 10)
		{
			cout << curr_index << endl;
			test_counter++;
		}
		word_counter = 0;
		input_sam_line_ss << input_sam_line;
		while(getline(input_sam_line_ss, sam_word_in_line, '\t')) //parse for genomic position
		{
			if (word_counter == 3)
			{
				curr_genomic_pos = stoi(sam_word_in_line);
			}
			word_counter++;
		}
		input_sam_line_ss.clear();

		//If index not in map, check to see if it's close enough to existing index;
		if(index_to_pos_in_vect_map.find(curr_index) == index_to_pos_in_vect_map.end())
		{
			ok_hamming_distance = false;
			if (mismatch_index_to_real_index_map.find(curr_index) != mismatch_index_to_real_index_map.end()){
				curr_index = mismatch_index_to_real_index_map[curr_index];
				ok_hamming_distance = true;
			}
			//ok_hamming_distance = true;
			/*ok_hamming_distance = false;
			vector_pos_counter = 0;
			while(ok_hamming_distance == false && vector_pos_counter < vector_of_indexes.size())
			{

				ok_hamming_distance = check_hamming_distance(vector_of_indexes.at(vector_pos_counter).index_seq, curr_index, 1);
				if (ok_hamming_distance == true)
				{
					curr_index = vector_of_indexes.at(vector_pos_counter).index_seq;
				}
				vector_pos_counter++;


			}
		*/
			if(ok_hamming_distance == false) //if index doesn't match map even w/ homing dist correction
			{
				output_sam_tossed_seq_file << input_sam_line << endl;
				if (check_if_index_in_map(bad_index_map, curr_index) == true) //UNTESTED
				{
					bad_index_map[curr_index]++;
				}
				else
				{
					bad_index_map[curr_index] = 1;
				}
			}
		}

		if (ok_hamming_distance == true && curr_genomic_pos - last_genomic_pos > 3)
		{
			vector_of_indexes.at(index_to_pos_in_vect_map[curr_index]).num_total_reads++;
			vector_of_indexes.at(index_to_pos_in_vect_map[curr_index]).num_unique_reads++;
			//output_sam_file << input_sam_line.substr(0, input_sam_line.find("_")) + "_" + curr_umi + "_" + curr_index + input_sam_line.substr(input_sam_line.find("_") + 20) << endl;
			output_sam_file << curr_index + "_" + curr_umi + input_sam_line.substr(input_sam_line.find("_") + 13) << endl;

			//recent_umi_index_combos_vect.clear();
			//recent_umi_index_combos_vect.push_back(curr_umi + "_" + curr_index);
			recent_umi_index_combos_vect_MAP.clear();
			recent_umi_index_combos_vect_MAP[curr_umi + "_" + curr_index] = curr_umi + "_" + curr_index;

		}
		else if (ok_hamming_distance == true && curr_genomic_pos - last_genomic_pos < 3)
		{
			if (check_hamming_distance(curr_index, last_index, 1) == true && check_hamming_distance(curr_umi, last_umi, 2) == true)
			{
				vector_of_indexes.at(index_to_pos_in_vect_map[curr_index]).num_total_reads++;
			}
			else
			{
				//already_in_vector = check_if_string_in_vector(recent_umi_index_combos_vect, curr_umi + "_" + curr_index);
				if(recent_umi_index_combos_vect_MAP.find(curr_umi + "_" + curr_index) != recent_umi_index_combos_vect_MAP.end()){
					already_in_vector = true;
				}
				//NEED TO ADD HAMMING DISTANCE HERE!!!!
				if (already_in_vector == false)
				{
					vector_of_indexes.at(index_to_pos_in_vect_map[curr_index]).num_total_reads++;
					vector_of_indexes.at(index_to_pos_in_vect_map[curr_index]).num_unique_reads++;
					//output_sam_file << input_sam_line.substr(0, input_sam_line.find("_")) + "_" + curr_umi + "_" + curr_index + input_sam_line.substr(input_sam_line.find("_") + 20) << endl;
					output_sam_file << curr_index + "_" + curr_umi + input_sam_line.substr(input_sam_line.find("_") + 13) << endl;
					//recent_umi_index_combos_vect.push_back(curr_umi + "_" + curr_index);
					recent_umi_index_combos_vect_MAP[curr_umi + "_" + curr_index] = curr_umi + "_" + curr_index;

				}
				else
				{
					vector_of_indexes.at(index_to_pos_in_vect_map[curr_index]).num_total_reads++;
					//cout << input_sam_line.substr(0, input_sam_line.find("_")) << "_" << curr_umi << "_" << curr_index << input_sam_line.substr(input_sam_line.find("_") + 20) << endl;
					already_in_vector = false;
				}

			}
		}
		last_index = curr_index;
		last_umi = curr_umi;
		last_genomic_pos = curr_genomic_pos;

		//cout << curr_umi << "+" << curr_index << "+" << curr_genomic_pos << endl;
	}

	ofstream(reads_and_dup_metrics_file);
	reads_and_dup_metrics_file.open(argv[5]);

	for(int i = 0; i < vector_of_indexes.size(); i++)
	{
		reads_and_dup_metrics_file << vector_of_indexes.at(i).index_seq << '\t' << vector_of_indexes.at(i).num_unique_reads << '\t' << vector_of_indexes.at(i).num_total_reads << endl;
	}

	//print bad index map
	ofstream(bad_index_counts_file);
	bad_index_counts_file.open(argv[6]);
	for(auto it = bad_index_map.begin(); it != bad_index_map.end(); it++)
	{
		bad_index_counts_file << it -> first << "\t" << it -> second << endl;
	}


	return 0;
}








