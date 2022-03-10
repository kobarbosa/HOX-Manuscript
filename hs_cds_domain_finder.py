#!/usr/bin/env python

import sys
import urllib
import urllib2
import json
import time
import xml.etree.ElementTree as element_tree
import textwrap

DELIMETER = "\t"


class UniprotRestClient(object):
	def __init__(self, server='http://www.uniprot.org/uniprot/', reqs_per_sec=15):
		self.server = server
		self.reqs_per_sec = reqs_per_sec
		self.req_count = 0
		self.last_req = 0


	def perform_rest_action(self, endpoint):
		if endpoint:
			try:
				endpoint = self.server + endpoint
				request = urllib2.Request(endpoint)
				response = urllib2.urlopen(request)
				content = response.read()
				print "Recovered data from endpoint: {0}   Status: {1}".format(endpoint, response.code)
				return content
			except urllib2.HTTPError, e:
				sys.stderr.write('Uniprot request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))
				return None


	def parse_xml_for_sequence(self, content):
		dom = element_tree.fromstring(content)
		entry_node = dom.find("{http://uniprot.org/uniprot}entry")
		
		sequence_node = entry_node.find("{http://uniprot.org/uniprot}sequence")
		sequence = sequence_node.text.strip().replace("\n", "")
		sequence_len = sequence_node.attrib["length"]
		return sequence, sequence_len


	def parse_xml_for_feature(self, content):
		dom = element_tree.fromstring(content)
		entry_node = dom.find("{http://uniprot.org/uniprot}entry")

		unwanted_feature_types = ("non-terminal residue","metal ion-binding site", "site", "modified residue", "splice variant", "sequence variant", "mutagenesis site")
		feature_nodes = entry_node.findall("{http://uniprot.org/uniprot}feature")
		for feature_node in feature_nodes:
			if feature_node.attrib["type"].strip() not in unwanted_feature_types:
				begin = feature_node.find("{http://uniprot.org/uniprot}location").find("{http://uniprot.org/uniprot}begin").attrib["position"]
				end = feature_node.find("{http://uniprot.org/uniprot}location").find("{http://uniprot.org/uniprot}end").attrib["position"]
				description = ""
				if "description" in feature_node.attrib:
					feature_node.attrib["description"].strip()
				yield feature_node.attrib["type"].strip(), description, begin.strip(), end.strip()


	def parse_gff_for_feature(self, content):
		unwanted_feature_types = ("Chain","Sequence conflict","Non-terminal residue","Metal binding", "Site", "Modified residue", "Natural variant", "Mutagenesis", "Alternative sequence", "Cross-link")
		for i, line in enumerate(content.split("\n")):
			if i <= 1:
				pass
			else:
				if line != '':
					feature_type = line.strip().split("\t")[2]
					begin = line.strip().split("\t")[3]
					end = line.strip().split("\t")[4]
					description = line.strip().split("\t")[8]
					if feature_type not in unwanted_feature_types:
						yield feature_type, description, begin, end




class EnsemblRestClient(object):
	def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
		self.server = server
		self.reqs_per_sec = reqs_per_sec
		self.req_count = 0
		self.last_req = 0

	def perform_rest_action(self, endpoint, hdrs=None, params=None):
		if hdrs is None:
			hdrs = {}

		if 'Content-Type' not in hdrs:
			hdrs['Content-Type'] = 'application/json'

		if params:
			endpoint += '?' + self.build_url_params(params)

		data = None

		# check if we need to rate limit ourselves
		if self.req_count >= self.reqs_per_sec:
			delta = time.time() - self.last_req
			if delta < 1:
				time.sleep(1 - deltagenes)
			self.last_req = time.time()
			self.req_count = 0

		try:
			request = urllib2.Request(self.server + endpoint, headers=hdrs)
			response = urllib2.urlopen(request)
			content = response.read()
			if content:
				data = json.loads(content)
				print "Recovered data from endpoint: {0}   Status: {1}".format(endpoint, response.code)
			self.req_count += 1

		except urllib2.HTTPError, e:
			# check if we are being rate limited by the server
			if e.code == 429:
				if 'Retry-After' in e.headers:
					retry = e.headers['Retry-After']
					time.sleep(float(retry))
					self.perform_rest_action(endpoint, hdrs, params)
			else:
				sys.stderr.write('Ensembl request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))
		return data


	def build_url_params(self, params):
		url_params = ""
		for i, key in enumerate(params.keys()):
			if i == 0:
				url_params += "{0}={1}".format(key, params[key])
			elif i > 0:
				url_params += ";{0}={1}".format(key, params[key])
		return url_params



	def get_variants(self, stable_id):
		variants = None
		variants = self.perform_rest_action(
				'/overlap/id/{0}'.format(stable_id),
				params={'feature': 'variation'}
		)
		return variants


	def get_ensembl_id(self, species, symbol):
		stable_id = None
		genes = self.perform_rest_action(
			'/xrefs/symbol/{0}/{1}'.format(species, symbol), 
			params={'object_type': 'gene'}
		)
		if genes:
			stable_id = genes[0]['id']
		return stable_id


	def get_cds(self, stable_id):
		cds = self.perform_rest_action(
				'/sequence/id/{0}'.format(stable_id),
				params={'type':'cds', 'content-type':'application/json', 'multiple_sequences':'1'}
			)
		return cds


	def get_protein(self, stable_id):
		cds = self.perform_rest_action(
				'/sequence/id/{0}'.format(stable_id),
				params={'type':'protein', 'content-type':'application/json', 'multiple_sequences':'1'}
			)
		return cds






def read_file(input_filename):
	uniprot_ensembl_id_map = {}
	file_map = {}
	with open(input_filename, "r") as input_file:
		for i, line in enumerate(input_file):
			if i == 0: pass
			else:
				line_parts = line.strip().split(DELIMETER)
				gene_id, ensembl_id, id_conversion, uniprot_id, uniprot_species_id, review_status, gene_description, gene_names, species, protein_len = line_parts

				uniprot_ensembl_id_map[uniprot_id] = ensembl_id
				file_map[uniprot_id] = (gene_id, id_conversion, uniprot_species_id, review_status, gene_description, gene_names, species, protein_len)
	print "Read {0} uniprot identifiers from the file!".format(str(len(file_map.keys())))
	return uniprot_ensembl_id_map, file_map


def get_subsequence(sequence, start, end):
	start = int(start) - 1
	end = int(end) - 1
	return sequence[start:end]


def get_uniprot_data(uniprot_id):
	uniprot_client = UniprotRestClient()
	content = uniprot_client.perform_rest_action(uniprot_id + ".xml")

	protein_seq, prot_len = uniprot_client.parse_xml_for_sequence(content)
	print protein_seq

	for feature_type, description, begin, end in uniprot_client.parse_xml_for_feature(content):
		print feature_type, description, begin, end, get_subsequence(protein_seq, begin, end)

	content = uniprot_client.perform_rest_action(uniprot_id + ".gff")

	for feature_type, description, begin, end in uniprot_client.parse_gff_for_feature(content):
		print feature_type, description, begin, end, get_subsequence(protein_seq, begin, end)
	return protein_seq


def match_protein_cds_sequences(protein_json, cds_json, uniprot_seq):
	matched_cds = None
	matched_protein_id = None
	matched_cds_id = None
	for i, entity in enumerate(protein_json):
		if entity['seq'] == uniprot_seq:
			matched_protein_id = entity['id']
			matched_cds = cds_json[i]['seq']
			matched_cds_id = cds_json[i]['id']
	return matched_cds, matched_cds_id, matched_protein_id


def find_domain_cds_sequence(protein_seq, domain_seq, cds_seq):
	domain_cds = None
	start = protein_seq.find(domain_seq)
	end = start + len(domain_seq)
	nuc_start = start*3
	nuc_end = end*3 
	return cds_seq[nuc_start:nuc_end]


def generate_fa_format(seq):
	limit = 50
	start = 0
	end = 50
	fa = ""
	while len(seq) > end:
		fa += seq[start:end] + "\n"
		start += limit
		end += limit
	fa += seq[start:len(seq)]
	return fa


def yield_domain_cds_for_id(uniprot_id, ensembl_id):

	uniprot_client = UniprotRestClient()
	content = uniprot_client.perform_rest_action(uniprot_id + ".xml")

	uniprot_protein_seq, prot_len = uniprot_client.parse_xml_for_sequence(content)
	if uniprot_protein_seq:
		print "Found uniprot protein of size: {0}".format(str(prot_len))

		client = EnsemblRestClient()
		cds_json = client.get_cds(ensembl_id)

		proteins_json = client.get_protein(ensembl_id)

		matched_cds, matched_cds_id, matched_protein_id = match_protein_cds_sequences(proteins_json, cds_json, uniprot_protein_seq)

		if not matched_cds:
			print "********* NO MATCHING CDS FOUND ***********"
		else:
			print "Found matching CDS sequence for protein"

			content = uniprot_client.perform_rest_action(uniprot_id + ".gff")

			for feature_type, description, begin, end in uniprot_client.parse_gff_for_feature(content):
				feature_seq = get_subsequence(uniprot_protein_seq, begin, end)
				domain_cds = find_domain_cds_sequence(uniprot_protein_seq, feature_seq, matched_cds)
				yield domain_cds, matched_cds_id, matched_protein_id, feature_type, description, feature_seq



def run_domain():

	input_filename = "/Users/kbalbosa/Documents/Lab - Personal /ubiquitin_library/human_ubiquitin_library.tab"
	uniprot_ensembl_id_map, input_map = read_file(input_filename)

	output_file = "/Users/kbalbosa/Documents/Lab - Personal /ubiquitin_library/hsubiquitin_libraryCDS.fa"

	gene_func_list = []
	with open(output_file, "w") as output:
		for uniprot_id in uniprot_ensembl_id_map.keys():

			ensembl_id = uniprot_ensembl_id_map[uniprot_id]
			gene_name = input_map[uniprot_id][0]

			print "working on gene: {0}".format(gene_name)

			feature_count = 0
			for domain_cds, matched_cds_id, matched_protein_id, feature_type, description, feature_seq in yield_domain_cds_for_id(uniprot_id, ensembl_id):
				if domain_cds and len(domain_cds) > 30:
					print matched_cds_id, matched_protein_id, feature_type, description, feature_seq
					feature_count +=1
					output.write(">{0}:{1}:{2}:{3}:{4}\n".format(gene_name, matched_cds_id, matched_protein_id, feature_type, description))
					output.write(generate_fa_format(domain_cds) + "\n")
					gene_func_list.append("{0}\t{1}\t{2}".format(gene_name, feature_type, description))
			print "Features found: {0}".format(feature_count)

			print "\n ************* \n"

	for i in gene_func_list:
		print i

	print "\nCOMPLETE!\n"


def run_cds_only():

	input_filename = "/Users/kbalbosa/Documents/Lab - Personal /ubiquitin_library/human_ubiquitin_library.tab"
	uniprot_ensembl_id_map, input_map = read_file(input_filename)

	output_file = "/Users/kbalbosa/Documents/Lab - Personal /ubiquitin_library/hsubiquitin_libraryCDS.fa"

	no_cds_list = []
	with open(output_file, "w") as output:
		for uniprot_id in uniprot_ensembl_id_map.keys():

			ensembl_id = uniprot_ensembl_id_map[uniprot_id]
			gene_name = input_map[uniprot_id][0]

			print "working on gene: {0}".format(gene_name)

			client = EnsemblRestClient()
			cds_json = client.get_cds(ensembl_id)

			cds_size = 0
			winning_index = None
			for i, entity in enumerate(cds_json):
				seq = entity['seq']
				if len(seq) > cds_size:
					winning_index = i

			sequence = cds_json[winning_index]['seq']
			transcript_id = cds_json[winning_index]['id']
			if sequence is None:
				no_cds_list.append(gene_name)
			else:
				output.write(">{0}:{1}:{2}\n".format(gene_name, ensembl_id, transcript_id))
				output.write(generate_fa_format(sequence) + "\n")

			print "\n ************* \n"

	for i in no_cds_list:
		print "Nothing found: {0}".format(i)

	print "\nCOMPLETE!\n"




if __name__ == '__main__':

	## run_domain takes the input and matches the protein/ensembl_cds sequences together with functional data
	run_domain()

	# this gets the largest cds sequence from ensembl
	#run_cds_only()
