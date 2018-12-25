#!/usr/bin/env python3

"""
This code is for comparing the estimated haplotype blocks with a grand truth.

usage: python3 compare.py head_truth_21 out_algorithm.txt

My first code in python.

sina majdiain
10 Dec

"""
from sys import argv
import numpy as np
from collections import Counter


def parse_hap_truth(file_name):
	""" extracting haplotype from haplotype truth file  in probhap website
	input: file name
	output: a dictionary, label:variantpostion , value: allele 0 or 1
	"""
	truth_file = open(file_name, "r")
	dic = {}
	if '.vcf' in file_name:
		for line in truth_file:
			if not line.startswith('#'):
				line_slice = line.strip().split("\t")
				GT = line_slice[9]
				if '|' in GT:
					variantpostion = int(line_slice[1])
					dic[variantpostion] = int(GT[0])
	else:
		for line in truth_file:  # 2	695745	1	0
			line_slice = line.strip().split("\t")
			label = int(line_slice[0])
			dic[label] = int(line_slice[3])
	return dic


def parse_vcf_freeb(file):
	""" extracting haplotype from haplotype truth file  in probhap website
	input: file name
	output: a dictionary, label:variantpostion , value: allele 0 or 1
	"""
	vcf_file = open(file, "r")
	list_position = []
	for line in vcf_file:
		if not line.startswith('#'):
			line_slice = line.strip().split("\t")
			list_position.append(int(line_slice[1]))
	return list_position


def parse_hap_estimated_blocks(file, list_position):

	""" extracting estimated haplotype blocks from  algorithm's output
	input: file name
	output: dictionary, label: block number,  value: a dictionary of idx and hap
	"""
	hap_block_file = open(file, "r")
	dic = {}
	dic_in = {}
	block_idx = 0
	start = True
	for line in hap_block_file:
			line = line.strip()
			if line.startswith('B'):  # new block started
				block_idx += 1
				if not start:
					dic[block_idx] = dic_in
				start = False
				dic_in = {}
			elif not line.startswith('*'):  # & ('-' not in line)
				line_slice = line.split("\t")  # 4	0	1	chr1	801628	C ...#hapcut2
				if '-' not in line_slice[1]:
					allele = int(line_slice[1])
					idx = int(line_slice[0])
					if len(list_position):
						var_position = list_position[idx - 1]
					else:
						var_position = idx  # for the current matlab output the indx is 1 or 2
						# var_position = int(line_slice[4])  # hapcut

					dic_in[var_position] = allele
	dic[block_idx] = dic_in
	return dic


def compare_dics(dic_truth, dic_block_estimated):
	""" comaring estimated block of haplotype with truth haplotype
	input: two dicinary
	output: some metrics
	"""
	counts = Counter()
	allele_origin = []
	allele_origin_dic = {}
	dic_block_estimated_in_truth = {}
	for idx_estimated, allele_estimated in dic_block_estimated.items():
		if idx_estimated in dic_truth:
			dic_block_estimated_in_truth[idx_estimated] = allele_estimated
			counts['block_length_pure'] += 1  # number of variants which is in both estimated and the truth
			if allele_estimated == dic_truth[idx_estimated]:  # truth consists of one allele which assumed as paternal
				allele_origin_dic[idx_estimated] = 'P'
				allele_origin.append('P')  # paternal
				counts['isfrom_P'] += 1
			else:
				allele_origin_dic[idx_estimated] = 'M'
				allele_origin.append('M')
				counts['isfrom_M'] += 1
	num_correct_alleles = max(counts['isfrom_P'], counts['isfrom_M'])

	for indx in range(len(allele_origin)):
		if indx != 0:  # for calculating number of switches, first allele is ignored.
			if allele_origin[indx] != allele_origin[indx-1]:
				counts['switch'] += 1
		if (indx != 0) and (indx != len(allele_origin)-1):  # for calculating short/long switches, first+second allele are ignored.
			if (allele_origin[indx] != allele_origin[indx-1]) and (allele_origin[indx] != allele_origin[indx+1]):
				counts['short_switch'] += 1  # if MPM or PMP happens, we call it short switch
			if (allele_origin[indx] != allele_origin[indx-1]) and (allele_origin[indx] == allele_origin[indx+1]):
				counts['long_switch'] += 1  # if MPP or PMM happens, we call it long switch

		position_estimated_in_truth = list(dic_block_estimated_in_truth.keys())

	if len(allele_origin_dic) == len(dic_block_estimated_in_truth):
		errorr = 1
		pass
	span_block_adjusted_list = []
	dic_blocks_no_switch = {}
	if (len(dic_block_estimated_in_truth) > 1) and (counts['switch'] > 0): # create blocks with no switch.
		start = True
		for idx_pos, position  in enumerate(position_estimated_in_truth):
			# for idx_estimated, allele_estimated in dic_block_estimated_in_truth.items():
			if start:  # for first of idx_estimated in this block
				dic_blocks_no_switch = {}
				dic_in = {}
				dic_in[position] = dic_block_estimated_in_truth[position]
				start = False
			else:  	# for rest of idx_estimated in this block
				if allele_origin_dic[position] == allele_origin_dic[position_estimated_in_truth[idx_pos-1]]:  # No switch
					dic_in[position] = dic_block_estimated_in_truth[position]
				else:  # close the dic_in and store it and create a new dic_in
					if len(dic_in) > 1:
						counts['Number_of_dicin'] += 1
						dic_blocks_no_switch[counts['Number_of_dicin']] = dic_in
					dic_in = {}
					dic_in[position] = dic_block_estimated_in_truth[position]
		# age akhari  in	 else bashe chi dobar dic_in add shode >> moshkeli  pish nemiad, hade aghal tool yek mishe
		if len(dic_in) > 1:
			counts['Number_of_dicin'] += 1  # for the last one
			dic_blocks_no_switch[counts['Number_of_dicin']] = dic_in

		# for ech block_in calculate span_block_adjusted
		if len(dic_blocks_no_switch):
			list_postion_truth = np.array(list(dic_truth.keys()))
			for idx_dic, dic_in_s in dic_blocks_no_switch.items():
				if len(dic_in_s):
					list_postion_block_in = list(dic_in_s.keys())
					span_block_in = list_postion_block_in[-1] - list_postion_block_in[0]
					start_position_in = list_postion_truth > list_postion_block_in[0] - 1
					end_position_in = list_postion_truth < list_postion_block_in[-1] + 1
					number_snps_truth_within_block_in = sum(
						[i and j for i, j in zip(start_position_in, end_position_in)])

					propotion_phased_in = float(number_snps_truth_within_block_in) / len(dic_in_s)
					span_block_adjusted_in = span_block_in * propotion_phased_in
				else:
					span_block_adjusted_in = 0
				span_block_adjusted_list.append(span_block_adjusted_in)

		# a =  dic_blocks_no_switch
	list_postion_block = list(dic_block_estimated.keys())
	if len(list_postion_block):
		span_block = list_postion_block[-1] - list_postion_block[0]  # duitama 2011
		list_postion_truth = np.array(list(dic_truth.keys()))
		start_position = list_postion_truth > list_postion_block[0]-1
		end_position = list_postion_truth < list_postion_block[-1]+1
		number_snps_truth_within_block = sum([i and j for i, j in zip(start_position, end_position)])
		if counts['block_length_pure'] != 0:
			propotion_phased = float(number_snps_truth_within_block)/counts['block_length_pure']
		else:
			propotion_phased = 0
		span_block_adjusted = span_block * propotion_phased
	else:
		span_block_adjusted = 0



	return [num_correct_alleles, counts['block_length_pure'], counts['switch'], counts['short_switch'], counts['long_switch'], span_block_adjusted, span_block_adjusted_list]


if __name__ == "__main__":
	address = argv[1]
	filename_hap_truth = address+argv[2]
	filename_hap_matlab = address+argv[3]
	if len(argv) > 4:
		filename_vcf = address+argv[4]
		list_position = parse_vcf_freeb(filename_vcf)
	else:
		list_position = []

	dic_hap_truth = parse_hap_truth(filename_hap_truth)
	dic_hap_estimated = parse_hap_estimated_blocks(filename_hap_matlab, list_position)


	blocks_num_correct_alleles = []
	blocks_length_pure = []
	blocks_num_switch = []
	blocks_num_switch_short = []
	blocks_num_switch_long = []
	span_blocks_adjusted = []
	span_block_adjusted_list_all = []
	for num_block, dic_block in dic_hap_estimated.items():  # enumerate
		block_num_correct_alleles, block_length_pure, block_num_switch, block_num_switch_short, block_num_switch_long, span_block_adjusted, span_block_adjusted_list = compare_dics(dic_hap_truth, dic_block)

		blocks_num_correct_alleles.append(block_num_correct_alleles)
		blocks_length_pure.append(block_length_pure)
		blocks_num_switch.append(block_num_switch)
		blocks_num_switch_short.append(block_num_switch_short)
		blocks_num_switch_long.append(block_num_switch_long)
		span_blocks_adjusted.append(span_block_adjusted)
		span_block_adjusted_list_all += span_block_adjusted_list
	print('number of block estimated', len(dic_hap_estimated))

	blocks_rr = [float(x)/y for x, y in zip(blocks_num_correct_alleles, blocks_length_pure) if y != 0]
	print('mean of rr all blocks is ', np.round(np.mean(blocks_rr), 4))
	print('number of block estimated containing estimation', len(blocks_length_pure))


	print('sum of block lengh is', sum(blocks_length_pure))
	print('mean of block lengh is', np.round(np.mean(blocks_length_pure), 4))
	print('coverage of genome is ', np.round(float(sum(blocks_length_pure))/len(dic_hap_truth),4))

	blocks_swer = [float(x) / y for x, y in zip(blocks_num_switch, blocks_length_pure) if y != 0]
	print('mean of swer  all blocks is ', np.round(np.mean(blocks_swer), 4))
	blocks_swer_short = [float(x) / y for x, y in zip(blocks_num_switch_short, blocks_length_pure) if y != 0]
	print('mean of swer short  all blocks is ', np.round(np.mean(blocks_swer_short), 4))
	blocks_swer_long = [float(x) / y for x, y in zip(blocks_num_switch_long, blocks_length_pure) if y != 0]
	print('mean of swer long all blocks is ', np.round(np.mean(blocks_swer_long), 4))
	print('sum of short switches ', sum(blocks_num_switch_short))
	print('sum of long switches ', sum(blocks_num_switch_long))
	print('sum of switches ', sum(blocks_num_switch))
	blocks_length_pure_sorted = np.sort(blocks_length_pure)
	# print(' N50 is ', blocks_length_pure_sorted[round(float(len(blocks_length_pure_sorted))/2)+1])

	#print(span_blocks_adjusted)
	span_blocks_adjusted_sorted = np.sort(span_blocks_adjusted)
	span_blocks_adjusted_sorted_revr = span_blocks_adjusted_sorted[::-1]
	print('AN50 is ', span_blocks_adjusted_sorted_revr[round(len(span_blocks_adjusted_sorted_revr)/2)]) # think more  ??

	span_block_adjusted_list_all_sorted = np.sort(span_block_adjusted_list_all)
	span_blocks_adjusted_sorted_all_revr = span_block_adjusted_list_all_sorted[::-1]
	print('QAN50 is ', span_blocks_adjusted_sorted_all_revr[round(len(span_blocks_adjusted_sorted_all_revr)/2)]) # think more  ??


	#  Remove blocks with length of less than two