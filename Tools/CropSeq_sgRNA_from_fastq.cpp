/*
 * 201222_Karina_CROPseq_Quick_Count_gRNAs.cpp
 *
 *  Created on: Dec 22, 2020
 *      Author: anna
 *
 *      gRNA_file.open(argv[1]);
 *      	fastq_R1.open(argv[2]);
		fastq_R2.open(argv[3]);
		list_of_failed_seqs.open(argv[4]);
		main_output.open(argv[5]);
			string upstream_seq = argv[6];
	string downstream_seq = argv[7];

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

struct cell{
	string cell_barcode;
	//map<string, int> gRNA_to_UMI_count;
	map<string, map<string,int> > gRNA_to_UMI_to_count;

};

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

	//string upstream_seq = "GCTCTAAAAC";
	//string downstream_seq = "CGGTGTTTCG";

	string upstream_seq = argv[6];
	string downstream_seq = argv[7];

	//read in gRNAs to expect
	ifstream(gRNA_file);
	gRNA_file.open(argv[1]);

	map<string, string> full_gRNA_to_name;
	map<string, string> start_gRNA_to_name;
	map<string, string> end_gRNA_to_name;

	string line;
	string elem;
	stringstream ss;

	string temp_name;
	string temp_gRNA;
	int counter;

	int gRNA_counter = 0;

	vector<string> temp_vect_of_1_bp_mismatches;

	getline(gRNA_file, line);

	while(getline(gRNA_file, line)){
		gRNA_counter++;
		counter = 0;
		ss << line;
		while(getline(ss, elem, '\t')){
			if (counter == 0){
				temp_name = elem;
				counter++;
			}
			else if(counter == 1){
				temp_gRNA = elem;
				full_gRNA_to_name[temp_gRNA] = temp_name;
				temp_vect_of_1_bp_mismatches = get_all_1_bp_mismatches(temp_gRNA);
				for(int i = 0; i < temp_vect_of_1_bp_mismatches.size(); i++){
					full_gRNA_to_name[temp_vect_of_1_bp_mismatches[i]] = temp_name;
				}
				start_gRNA_to_name[temp_gRNA.substr(0,10)] = temp_name;
				end_gRNA_to_name[temp_gRNA.substr(10,10)] = temp_name;
				counter++;
			}
		}
		ss.clear();
	}


	cout << "Total # gRNAs" << gRNA_counter << endl;
	cout << "Full_map_size: " << full_gRNA_to_name.size() << endl;
	cout << "start_map_size: " << start_gRNA_to_name.size() << endl;
	cout << "end_map_size: " << end_gRNA_to_name.size() << endl;



	ifstream(fastq_R1);
	ifstream(fastq_R2);

	fastq_R1.open(argv[2]);
	fastq_R2.open(argv[3]);

	string line2;
	string temp_cell_barcode;
	string temp_umi;
	string temp_read;

	int temp_pos_upsteam_seq;
	int temp_pos_downstream_seq;

	string temp_seq_after_upstream_seq;
	string temp_seq_before_downstream_seq;
	//string temp_gRNA;
	string temp_gRNA_name;

	ofstream(list_of_failed_seqs);
	list_of_failed_seqs.open(argv[4]);

	counter = 0;

	cell temp_cell;
	vector<cell> vector_of_cells;
	map<string,int> map_of_cells_to_pos_in_vect;
	int temp_pos_of_cell;

	while(getline(fastq_R1, line)){
		getline(fastq_R2, line2);
		if (counter == 0){
			counter ++;
		}
		else if (counter == 1){
			counter++;
			temp_gRNA = "";
			temp_gRNA_name = "";
			temp_cell_barcode = "";
			temp_umi = "";
			temp_pos_upsteam_seq = line2.find(upstream_seq);
			temp_pos_downstream_seq = line2.find(downstream_seq);
			if (temp_pos_upsteam_seq != string::npos || temp_pos_downstream_seq != string::npos){ //if upstream of downstream seq found
				//cout << "Found a hit??" << endl;
				//cout << line2 << endl;
				if(temp_pos_upsteam_seq != string::npos){
					temp_seq_after_upstream_seq = line2.substr(temp_pos_upsteam_seq+upstream_seq.size());
					if(temp_seq_after_upstream_seq.size() >= 20){
						temp_gRNA = temp_seq_after_upstream_seq.substr(0,20);
						if (full_gRNA_to_name.find(temp_gRNA) != full_gRNA_to_name.end()){
							temp_gRNA_name = full_gRNA_to_name[temp_gRNA];
							//HERE WE CAN ADD HAMMING DISTANCE
						}
					}
					if (temp_gRNA_name == "" && temp_seq_after_upstream_seq.size() >= 10){
						temp_gRNA = temp_seq_after_upstream_seq.substr(0,10);
						if (start_gRNA_to_name.find(temp_gRNA) != start_gRNA_to_name.end()){
							temp_gRNA_name = start_gRNA_to_name[temp_gRNA];
							//HERE WE CAN ADD HAMMING DISTANCE
						}
					}
				}

				if (temp_gRNA_name == "" && temp_pos_downstream_seq != string::npos){
					temp_seq_before_downstream_seq = line2.substr(0,temp_pos_downstream_seq);
					if(temp_seq_before_downstream_seq.size() >=20){
						temp_gRNA = temp_seq_before_downstream_seq.substr(temp_seq_before_downstream_seq.size()-20);
						if (full_gRNA_to_name.find(temp_gRNA) != full_gRNA_to_name.end()){
							temp_gRNA_name = full_gRNA_to_name[temp_gRNA];
							//HERE WE CAN ADD HAMMING DISTANCE
						}
					}
					if(temp_gRNA_name == "" && temp_seq_before_downstream_seq.size() >=10){
						temp_gRNA = temp_seq_before_downstream_seq.substr(temp_seq_before_downstream_seq.size()-10);
						if (end_gRNA_to_name.find(temp_gRNA) != end_gRNA_to_name.end()){
							temp_gRNA_name = end_gRNA_to_name[temp_gRNA];
						}
					}
				}

				if (temp_gRNA_name == ""){ //still
					list_of_failed_seqs << line2 << endl;

				}
				else {
					temp_cell_barcode = line.substr(0,16);
					temp_umi = line.substr(15,12);
					//cout << "gRNA: " << temp_gRNA_name << endl;
					//cout << "cell: " << temp_cell_barcode << endl;
					//cout << "umi: " << temp_umi << endl;

					if(map_of_cells_to_pos_in_vect.find(temp_cell_barcode) == map_of_cells_to_pos_in_vect.end()){
						temp_cell.cell_barcode = temp_cell_barcode;
						temp_cell.gRNA_to_UMI_to_count[temp_gRNA_name][temp_umi] = 1;
						vector_of_cells.push_back(temp_cell);
						map_of_cells_to_pos_in_vect[temp_cell_barcode] = vector_of_cells.size()-1;
						temp_cell = cell();
					}
					else {
						temp_pos_of_cell = map_of_cells_to_pos_in_vect[temp_cell_barcode];
						if(vector_of_cells[temp_pos_of_cell].gRNA_to_UMI_to_count.find(temp_gRNA_name) == vector_of_cells[temp_pos_of_cell].gRNA_to_UMI_to_count.end()){
							vector_of_cells[temp_pos_of_cell].gRNA_to_UMI_to_count[temp_gRNA_name][temp_umi] = 1;
						}
						else { //if gRNA already in map
							if(vector_of_cells[temp_pos_of_cell].gRNA_to_UMI_to_count[temp_gRNA_name].find(temp_umi) == vector_of_cells[temp_pos_of_cell].gRNA_to_UMI_to_count[temp_gRNA_name].end()){
								vector_of_cells[temp_pos_of_cell].gRNA_to_UMI_to_count[temp_gRNA_name][temp_umi] = 1;
							}
							else {
								vector_of_cells[temp_pos_of_cell].gRNA_to_UMI_to_count[temp_gRNA_name][temp_umi]++;
							}
						}
					}
				}
			}
		}
		else if (counter == 2){
			counter++;
		}
		else if (counter == 3){
			counter = 0;
		}
	}

	ofstream(main_output);
	main_output.open(argv[5]);

	int temp_num_umis;

	for (int i = 0; i < vector_of_cells.size(); i++){
		for(auto it = vector_of_cells[i].gRNA_to_UMI_to_count.begin(); it != vector_of_cells[i].gRNA_to_UMI_to_count.end(); it++){
			temp_gRNA = it -> first;
			temp_num_umis = it -> second.size();
			main_output << vector_of_cells[i].cell_barcode << "\t" << temp_gRNA << "\t" << temp_num_umis << "\t";
			for(auto it2 = it -> second.begin(); it2 != it -> second.end(); it2++){
				main_output << it2 -> first << "_" << it2 -> second << ";";
			}
			main_output << endl;
		}
	}

	return 0;
}




