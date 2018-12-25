#!/usr/bin/env python3

"""
This code is for comparing the estimated diploid haplotype blocks with a grand truth.

Fosmid  > indx
python3  compare.py "/SinaMc/code1new/fosmid/folder/"  "chr21.valid.master"  "chr21_hap_opt.txt" "chr21.hap"


Ashkenazim > variant
python3  compare.py "/SinaMc/code1new/Ashkenazim/" "HG002_phased_22.vcf" "HG002_hap_opt_60.txt" "haplotype_hapcut_60x" "HG002_22.vcf"






usage: python3  compare.py  adress_files truth estimated_haplotype [vcffile]
For comparing based on  variation position, assign the  vcffile as the fourth

Sina Majdiain
Start: 10 Dec

"""
from sys import argv
import numpy as np
from collections import Counter


def parse_hap_truth(file_name):
	""" extracting haplotype from haplotype truth file just heterozygous variants.
	input: file name
	output: a dictionary, key: position of variant  value: one allele (0 or 1)

	Since it is heterozygous, one allel is enough for comapring.
	The input haplotype if contain   .vcf is considered as this option
	a vcf file inwhich GT id contain | . This used in HapCHAT paper and also GIAB project. Also mentioned in part 1.4.2  of
	"The Variant Call Format (VCF) Version 4.1 Specification"

	20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51
	20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50

	else it is considered  as a Kuleshov (probhap) format for grand truth:

	2	695745	1	0
	3	766409	0	1
	4	801628	0	1
	"""
	truth_file = open(file_name, "r")
	dic_hap_truth = {}
	if '.vcf' in file_name:
		for line in truth_file:
			if not line.startswith('#'):  # GT:GQ:DP:HQ 0|0:48:1:51,51
				slices_line = line.strip().split("\t")
				gt_value = slices_line[9]
				if '|' in gt_value:  # just keep phased variants. some variant may not phased.
					postion_variant = int(slices_line[1])
					dic_hap_truth[postion_variant] = int(gt_value[0])
	else:
		for line in truth_file:  # 2	695745	1	0
			slices_line = line.strip().split("\t")
			idx_variant = int(slices_line[0])   # for kuleshov data, I don't know the position of all variants,
			# I used the indx instead of postion of variant.
			dic_hap_truth[idx_variant] = int(slices_line[2])
			# postion_variant = int(slices_line[1])  # if we know all positon we can use this,
			# dic_hap_truth[postion_variant] = int(slices_line[2])
	return dic_hap_truth


def parse_position_vcf(file):
	""" extracting position vcf file

	input: vcf file name
	output: a list of  variant postion
	"""
	vcf_file = open(file, "r")
	list_postion_variant = []
	for line in vcf_file:
		if not line.startswith('#'):  # 22	17056280	.	C	A	.	.	AC=2;AF=1.00;AN=2
			slices_line = line.strip().split("\t")
			postion_variant = int(slices_line[1])
			list_postion_variant.append(postion_variant)
	return list_postion_variant


def parse_hap_estimated_blocks(file, list_position):

	""" extracting estimated haplotype blocks from  algorithm's output file

	input: file name
	output: dictionary, label: block number,  value: a dictionary of idx and hap

	the inptut can be like this
	BLOCK	2	7
	2	0
	3	1

	or

	BLOCK: offset: 4 len: 14 phased: 14 SPAN: 25449 fragments 14
	4	0	1	chr1	801628	C	T	0/1:69:28	0	.	31.00
	5	0	1	chr1	805678	A	T	0/1:21:70	0	.	100.00
	"""
	hap_blocks = open(file, "r")
	set_indx = set()
	dic_blocks = {}  # dic_blocks = {1: dic_block1, 2: dic_block2}
	dic_block = {}  # dic_block = {2:0 , 3:1 , 28:1}
	block_idx = 0
	start = True
	for line in hap_blocks:
			line = line.strip()
			if line.startswith('B'):  # new block started
				if not start:
					block_idx += 1
					dic_blocks[block_idx] = dic_block
				start = False
				dic_block = {}
			elif not line.startswith('*'):
				slices_line = line.split("\t")
				if '-' not in slices_line[1]:
					allele = int(slices_line[1])
					idx = int(slices_line[0])
					if not len(list_position):  # for comparing based on index, not variation position
						dic_block[idx] = allele
						set_indx.add(idx)
					else:  # for comparing based on variation position
						var_position = list_position[idx - 1]
						set_indx.add(var_position)
						dic_block[var_position] = allele
						# For future  if len(slices_line) > 3:  # hapcut2: 4	0	1	chr1	801628	C ...
						# var_position = int(line_slice[4])
	block_idx += 1
	dic_blocks[block_idx] = dic_block
	return dic_blocks, set_indx


def compare_dics(dic_truth, dic_block_estimated):
	""" comaring estimated block of haplotype with truth haplotype

	input: two dictionaries
	output: some metrics
	"""
	counts = Counter()
	origin_allele = []
	origin_allele_dic = {}
	dic_block_estimated_truth = {}
	for idx_estimated, allele_estimated in dic_block_estimated.items(): # idx_estimated can be either indx or position
		if idx_estimated in dic_truth:
			dic_block_estimated_truth[idx_estimated] = allele_estimated  # create a dic of  variants in both estimated and truth
			if allele_estimated == dic_truth[idx_estimated]:  # truth consists of one allele which assumed as paternal. hetro
				origin_allele_dic[idx_estimated] = 'P'
				origin_allele.append('P')  # paternal
				counts['isfrom_P'] += 1
			else:
				origin_allele_dic[idx_estimated] = 'M'
				origin_allele.append('M')
				counts['isfrom_M'] += 1
	num_correct_alleles = max(counts['isfrom_P'], counts['isfrom_M'])

	if num_correct_alleles > 1:
		for indx_truth in range(len(origin_allele)):
			if indx_truth != 0:  # for calculating number of switches, first allele is ignored.
				if origin_allele[indx_truth] != origin_allele[indx_truth-1]:
					counts['switch'] += 1
			if (indx_truth != 0) and (indx_truth != len(origin_allele)-1):  # for calculating sh/lo switches, first+last allele are ignored.
				if (origin_allele[indx_truth] != origin_allele[indx_truth-1]) and (origin_allele[indx_truth] != origin_allele[indx_truth+1]):
					counts['short_switch'] += 1  # if MPM or PMP happens, we call it short switch
				if (origin_allele[indx_truth] != origin_allele[indx_truth-1]) and (origin_allele[indx_truth] == origin_allele[indx_truth+1]):
					counts['long_switch'] += 1  # if MPP or PMM happens, we call it long switch
	metrics = [num_correct_alleles, counts['switch'], counts['short_switch'], counts['long_switch']]
	return [dic_block_estimated_truth, origin_allele_dic, metrics]





def compare_dics_common(dic_truth, dic_block_estimated, set_idx2):
	""" comaring estimated block of haplotype with truth haplotype
	dic_block_estimated is from algorithm one
	input: two dictionaries
	output: some metrics
	"""
	counts = Counter()
	origin_allele = []
	origin_allele_dic = {}
	dic_block_estimated_truth = {}
	for idx_estimated, allele_estimated in dic_block_estimated.items():
		if (idx_estimated in dic_truth) & (idx_estimated in set_idx2):
			dic_block_estimated_truth[idx_estimated] = allele_estimated  # create a dic of  variants in both estimated and truth
			if allele_estimated == dic_truth[idx_estimated]:  # truth consists of one allele which assumed as paternal. hetro
				origin_allele_dic[idx_estimated] = 'P'
				origin_allele.append('P')  # paternal
				counts['isfrom_P'] += 1
			else:
				origin_allele_dic[idx_estimated] = 'M'
				origin_allele.append('M')
				counts['isfrom_M'] += 1
	num_correct_alleles = max(counts['isfrom_P'], counts['isfrom_M'])

	if num_correct_alleles > 1:
		for indx_truth in range(len(origin_allele)):
			if indx_truth != 0:  # for calculating number of switches, first allele is ignored.
				if origin_allele[indx_truth] != origin_allele[indx_truth-1]:
					counts['switch'] += 1
			if (indx_truth != 0) and (indx_truth != len(origin_allele)-1):  # for calculating sh/lo switches, first+last allele are ignored.
				if (origin_allele[indx_truth] != origin_allele[indx_truth-1]) and (origin_allele[indx_truth] != origin_allele[indx_truth+1]):
					counts['short_switch'] += 1  # if MPM or PMP happens, we call it short switch
				if (origin_allele[indx_truth] != origin_allele[indx_truth-1]) and (origin_allele[indx_truth] == origin_allele[indx_truth+1]):
					counts['long_switch'] += 1  # if MPP or PMM happens, we call it long switch
	metrics = [num_correct_alleles, counts['switch'], counts['short_switch'], counts['long_switch']]
	return [dic_block_estimated_truth, origin_allele_dic, metrics]






# def new(dic_truth, dic_block_estimated):
#
#
# 	if len(allele_origin_dic) == len(dic_block_estimated_in_truth):
# 		errorr = 1
# 		pass
# 	span_block_adjusted_list = []
# 	dic_blocks_no_switch = {}
# 	if (len(dic_block_estimated_in_truth) > 1) and (counts['switch'] > 0): # create blocks with no switch.
# 		start = True
# 		for idx_pos, position  in enumerate(position_estimated_in_truth):
# 			# for idx_estimated, allele_estimated in dic_block_estimated_in_truth.items():
# 			if start:  # for first of idx_estimated in this block
# 				dic_blocks_no_switch = {}
# 				dic_in = {}
# 				dic_in[position] = dic_block_estimated_in_truth[position]
# 				start = False
# 			else:  	# for rest of idx_estimated in this block
# 				if allele_origin_dic[position] == allele_origin_dic[position_estimated_in_truth[idx_pos-1]]:  # No switch
# 					dic_in[position] = dic_block_estimated_in_truth[position]
# 				else:  # close the dic_in and store it and create a new dic_in
# 					if len(dic_in) > 1:
# 						counts['Number_of_dicin'] += 1
# 						dic_blocks_no_switch[counts['Number_of_dicin']] = dic_in
# 					dic_in = {}
# 					dic_in[position] = dic_block_estimated_in_truth[position]
# 		# age akhari  in	 else bashe chi dobar dic_in add shode >> moshkeli  pish nemiad, hade aghal tool yek mishe
# 		if len(dic_in) > 1:
# 			counts['Number_of_dicin'] += 1  # for the last one
# 			dic_blocks_no_switch[counts['Number_of_dicin']] = dic_in
#
# 		# for ech block_in calculate span_block_adjusted
# 		if len(dic_blocks_no_switch):
# 			list_postion_truth = np.array(list(dic_truth.keys()))
# 			for idx_dic, dic_in_s in dic_blocks_no_switch.items():
# 				if len(dic_in_s):
# 					list_postion_block_in = list(dic_in_s.keys())
# 					span_block_in = list_postion_block_in[-1] - list_postion_block_in[0]
# 					start_position_in = list_postion_truth > list_postion_block_in[0] - 1
# 					end_position_in = list_postion_truth < list_postion_block_in[-1] + 1
# 					number_snps_truth_within_block_in = sum(
# 						[i and j for i, j in zip(start_position_in, end_position_in)])
#
# 					propotion_phased_in = float(number_snps_truth_within_block_in) / len(dic_in_s)
# 					span_block_adjusted_in = span_block_in * propotion_phased_in
# 				else:
# 					span_block_adjusted_in = 0
# 				span_block_adjusted_list.append(span_block_adjusted_in)
#
# 		# a =  dic_blocks_no_switch
# 	list_postion_block = list(dic_block_estimated.keys())
# 	if len(list_postion_block):
# 		span_block = list_postion_block[-1] - list_postion_block[0]  # duitama 2011
# 		list_postion_truth = np.array(list(dic_truth.keys()))
# 		start_position = list_postion_truth > list_postion_block[0]-1
# 		end_position = list_postion_truth < list_postion_block[-1]+1
# 		number_snps_truth_within_block = sum([i and j for i, j in zip(start_position, end_position)])
# 		if counts['block_length_pure'] != 0:
# 			propotion_phased = float(number_snps_truth_within_block)/counts['block_length_pure']
# 		else:
# 			propotion_phased = 0
# 		span_block_adjusted = span_block * propotion_phased
# 	else:
# 		span_block_adjusted = 0
#
#
#
# 	return [num_correct_alleles, , counts['switch'], counts['short_switch'], counts['long_switch'], span_block_adjusted, span_block_adjusted_list]


if __name__ == "__main__":
	address = argv[1]
	filename_hap_truth = address+argv[2]
	filename_hap_estimated_alg1 = address+argv[3]
	filename_hap_estimated_alg2 = address+argv[4]

	if len(argv) < 6:  # for comparing based on  variation position, assign the  filename_vcf # i know that hapcut may contain
		list_position = []
	else:
		filename_vcf = address + argv[5]
		list_position = parse_position_vcf(filename_vcf)

	dic_hap_truth = parse_hap_truth(filename_hap_truth)
	dic_hap_estimated_alg1, set_idx1 = parse_hap_estimated_blocks(filename_hap_estimated_alg1, list_position)
	dic_hap_estimated_alg2, set_idx2 = parse_hap_estimated_blocks(filename_hap_estimated_alg2, list_position)

	# set_position_algortihm2 = set(dic_hap_estimated_alg2.keys)



	dic_metrics1 = {'length_pure':[]}
	dic_metrics1['correct_alleles'] = []
	dic_metrics1['switch'] = []
	dic_metrics1['switch_short'] = []
	dic_metrics1['switch_long'] = []

	dic_metrics2 = {'length_pure': []}
	dic_metrics2['correct_alleles'] = []
	dic_metrics2['switch'] = []
	dic_metrics2['switch_short'] = []
	dic_metrics2['switch_long'] = []

	dic_metrics_cm = {'length_pure': []}
	dic_metrics_cm['correct_alleles'] = []
	dic_metrics_cm['switch'] = []
	dic_metrics_cm['switch_short'] = []
	dic_metrics_cm['switch_long'] = []
	# span_blocks_adjusted = []
	# span_block_adjusted_list_all = []

	for num_block, dic_block in dic_hap_estimated_alg1.items():
		[dic_block_estimated_truth, allele_origin_dic, metrics] = compare_dics(dic_hap_truth, dic_block)
		[dic_block_estimated_truth_cm, allele_origin_dic_cm, metrics_cm]= compare_dics_common(dic_hap_truth, dic_block, set_idx2)
		# metrics = [num_correct_alleles, counts['switch'], counts['short_switch'], counts['long_switch']]
		dic_metrics1['length_pure'].append(len(dic_block_estimated_truth))
		dic_metrics1['correct_alleles'].append(metrics[0])
		dic_metrics1['switch'].append(metrics[1])
		dic_metrics1['switch_short'].append(metrics[2])
		dic_metrics1['switch_long'].append(metrics[3])

		dic_metrics_cm['length_pure'].append(len(dic_block_estimated_truth_cm))
		dic_metrics_cm['correct_alleles'].append(metrics_cm[0])
		dic_metrics_cm['switch'].append(metrics_cm[1])
		dic_metrics_cm['switch_short'].append(metrics_cm[2])
		dic_metrics_cm['switch_long'].append(metrics_cm[3])

	for num_block, dic_block in dic_hap_estimated_alg2.items():
		[dic_block_estimated_truth, allele_origin_dic, metrics] = compare_dics(dic_hap_truth, dic_block)
		dic_metrics2['length_pure'].append(len(dic_block_estimated_truth))
		dic_metrics2['correct_alleles'].append(metrics[0])
		dic_metrics2['switch'].append(metrics[1])
		dic_metrics2['switch_short'].append(metrics[2])
		dic_metrics2['switch_long'].append(metrics[3])


	print('******** First  OPT ')
	blocks_rr1 = [float(x)/y for x, y in zip(dic_metrics1['correct_alleles'], dic_metrics1['length_pure']) if y != 0]
	print('mean of rr all blocks is ', np.round(np.mean(blocks_rr1), 4))
	blocks_swer1 = [float(x) / y for x, y in zip(dic_metrics1['switch'], dic_metrics1['length_pure']) if y != 0]
	print('mean of swer  all blocks is ', np.round(np.mean(blocks_swer1), 4))
	blocks_swer_short1 = [float(x) / y for x, y in zip(dic_metrics1['switch_short'], dic_metrics1['length_pure']) if y != 0]
	print('mean of swer short  all blocks is ', np.round(np.mean(blocks_swer_short1), 4))
	blocks_swer_long1 = [float(x) / y for x, y in zip(dic_metrics1['switch_long'], dic_metrics1['length_pure']) if y != 0]
	print('mean of swer long all blocks is ', np.round(np.mean(blocks_swer_long1), 4))

	print('******** Second  CUT')
	blocks_rr2 = [float(x) / y for x, y in zip(dic_metrics2['correct_alleles'], dic_metrics2['length_pure']) if y != 0]
	print('mean of rr all blocks is ', np.round(np.mean(blocks_rr2), 4))
	blocks_swer2 = [float(x) / y for x, y in zip(dic_metrics2['switch'], dic_metrics2['length_pure']) if y != 0]
	print('mean of swer  all blocks is ', np.round(np.mean(blocks_swer2), 4))
	blocks_swer_short2 = [float(x) / y for x, y in zip(dic_metrics2['switch_short'], dic_metrics2['length_pure']) if y != 0]
	print('mean of swer short  all blocks is ', np.round(np.mean(blocks_swer_short2), 4))
	blocks_swer_long2 = [float(x) / y for x, y in zip(dic_metrics2['switch_long'], dic_metrics2['length_pure']) if y != 0]
	print('mean of swer long all blocks is ', np.round(np.mean(blocks_swer_long2), 4))

	print ('******** Common for OPT')
	blocks_rr_cm = [float(x) / y for x, y in zip(dic_metrics_cm['correct_alleles'], dic_metrics_cm['length_pure']) if y != 0]
	print('mean of rr all blocks is ', np.round(np.mean(blocks_rr_cm), 4))
	blocks_swer_cm = [float(x) / y for x, y in zip(dic_metrics_cm['switch'], dic_metrics_cm['length_pure']) if y != 0]
	print('mean of swer  all blocks is ', np.round(np.mean(blocks_swer_cm), 4))
	blocks_swer_short_cm = [float(x) / y for x, y in zip(dic_metrics_cm['switch_short'], dic_metrics_cm['length_pure']) if y != 0]
	print('mean of swer short  all blocks is ', np.round(np.mean(blocks_swer_short_cm), 4))
	blocks_swer_long_cm = [float(x) / y for x, y in zip(dic_metrics_cm['switch_long'], dic_metrics_cm['length_pure']) if y != 0]
	print('mean of swer long all blocks is ', np.round(np.mean(blocks_swer_long_cm), 4))


	# block_num_correct_alleles, block_length_pure, block_num_switch, block_num_switch_short, block_num_switch_long, span_block_adjusted, span_block_adjusted_list



		# blocks_num_correct_alleles.append(block_num_correct_alleles)
		# blocks_length_pure.append(block_length_pure)
		# blocks_num_switch.append(block_num_switch)
		# blocks_num_switch_short.append(block_num_switch_short)
		# blocks_num_switch_long.append(block_num_switch_long)
		# span_blocks_adjusted.append(span_block_adjusted)
		# span_block_adjusted_list_all += span_block_adjusted_list
	# print('number of block estimated', len(dic_hap_estimated))
	#

# print('number of block estimated containing estimation', len(blocks_length_pure))
	#
	#
	# print('sum of block lengh is', sum(blocks_length_pure))
	# print('mean of block lengh is', np.round(np.mean(blocks_length_pure), 4))
	# print('coverage of genome is ', np.round(float(sum(blocks_length_pure))/len(dic_hap_truth),4))
	#

	# print('sum of short switches ', sum(blocks_num_switch_short))
	# print('sum of long switches ', sum(blocks_num_switch_long))
	# print('sum of switches ', sum(blocks_num_switch))
	# blocks_length_pure_sorted = np.sort(blocks_length_pure)
	# # print(' N50 is ', blocks_length_pure_sorted[round(float(len(blocks_length_pure_sorted))/2)+1])
	#
	# #print(span_blocks_adjusted)
	# span_blocks_adjusted_sorted = np.sort(span_blocks_adjusted)
	# span_blocks_adjusted_sorted_revr = span_blocks_adjusted_sorted[::-1]
	# print('AN50 is ', span_blocks_adjusted_sorted_revr[round(len(span_blocks_adjusted_sorted_revr)/2)]) # think more  ??
	#
	# span_block_adjusted_list_all_sorted = np.sort(span_block_adjusted_list_all)
	# span_blocks_adjusted_sorted_all_revr = span_block_adjusted_list_all_sorted[::-1]
	# print('QAN50 is ', span_blocks_adjusted_sorted_all_revr[round(len(span_blocks_adjusted_sorted_all_revr)/2)]) # think more  ??


	#  Remove blocks with length of less than two
