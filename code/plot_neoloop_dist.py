"""
HiChIP analysis - (Neo)loop distribution on different classes of focal amplifications (ecDNA/BFB/Complex/Linear)
Usage: python(3) plot_neoloop_dist.py <HiChIP to WGS map (hic2sample_all.txt)>
		<List of control samples, for each WGS (control_list_all.tsv)>
		<Directory of case neoloop calls>
		<Directory of control neoloop calls>
		<Directory of AC classification bed files>
		<AC classification summaries *_amplicon_classification_profiles.tsv>
		<Directory of AC classification cycle files>
		<Output (figure)>
		
Output: num_neoloops_per_mbp.tsv and the violin plot 
"""
import sys
import os
from os import listdir
import numpy as np
import scipy.stats
from scipy.stats import ranksums
from scipy.stats import ttest_ind

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.box(False)
from pylab import rcParams
rcParams['figure.figsize'] = [16, 12]
rcParams['pdf.fonttype'] = 42
mpl.rc('xtick', labelsize = 30) 
mpl.rc('ytick', labelsize = 30) 


def check_int(chr1, chr2, s1, s2, s3, s4, sint):
	r1 = -1
	r2 = -1
	for si in sint:
		mid1 = (s1 + s2) / 2
		mid2 = (s3 + s4) / 2
		if chr1 == si[0] and mid1 >= si[1] and mid1 <= si[2]:
			r1 = sint.index(si)
		if chr2 == si[0] and mid2 >= si[1] and mid2 <= si[2]:
			r2 = sint.index(si)
	return r1, r2


def check_int_overlap(sint1, sint2):
	for si1 in sint1:
		for si2 in sint2:
			if si1[0] == si2[0] and si1[1] <= si2[2] and si2[1] <= si1[2]:
				return True
	return False 


def search_cont(loop, cont_loop_list):
	mid_case_1 = (int(loop[1]) + int(loop[2])) / 2
	mid_case_2 = (int(loop[4]) + int(loop[5])) / 2
	for cont_loop in cont_loop_list:
		mid_cont_1 = (int(cont_loop[1]) + int(cont_loop[2])) / 2
		mid_cont_2 = (int(cont_loop[4]) + int(cont_loop[5])) / 2
		if cont_loop[0] == loop[0] and cont_loop[3] == loop[3] and mid_cont_1 == mid_case_1 and mid_case_2 == mid_cont_2:
			return True
		if cont_loop[0] == loop[3] and cont_loop[3] == loop[0] and mid_cont_1 == mid_case_2 and mid_case_1 == mid_cont_2:
			return True
	return False


if __name__ == '__main__':

	# Obtain the list of WGS samples and their matching HiChIP from file 
	sample_hic = []
	fp = open(sys.argv[1], 'r') # Path to "hic2sample_all.txt"
	for line in fp:
		line = line.strip()
		s = line.split()
		sample_hic.append([s[1], s[0]])
	fp.close()

	# List of control samples for each WGS sample
	control_list = dict()
	fp = open(sys.argv[2], 'r') # Path to "control_list_all.tsv"
	for line in fp:
		line = line.strip()
		s = line.split()
		try:
			control_list[s[0]].append(s[1])
		except:
			control_list[s[0]] = [s[1]]
	fp.close()

	# Read in neo loops
	neo_loop_dict = dict()
	neo_loop_dict_control = dict()
	for [wgs, hic] in sample_hic:
		for fn in listdir(sys.argv[3]): # Directory of case neoloop calls
			if hic in fn and 'neo-loops' in fn:
				fp = open(sys.argv[3] + fn, 'r')
				neo_loop_dict[(wgs, hic)] = []
				neo_loop_dict_control[(wgs, hic)] = []
				for line in fp:
					line = line.strip()
					s = line.split()
					if wgs in s[-1] or s[-1][0] == 'A' or s[-1][0] == 'C':
						neo_loop_dict[(wgs, hic)].append(s)
				fp.close()
		for cont_fn in listdir(sys.argv[4]): # Directory of control neoloop calls
			if hic in cont_fn and 'neo-loops' in cont_fn:
				fp_cont = open(sys.argv[4] + cont_fn, 'r')
				for line in fp_cont:
					line = line.strip()
					s = line.split()
					if s[-1] in control_list[wgs]:
						neo_loop_dict_control[(wgs, hic)].append(s)
				fp_cont.close()	
	print("Identified all %d WGS samples having NeoLoop results available." %len(neo_loop_dict))

	# Read AA classification results
	amplicon_dict = dict()
	for fn in listdir(sys.argv[5]):
		a_idx = fn.index('amplicon')
		sample_name = fn[:a_idx - 1]
		s = fn.split('_')
		amplicon_idx = 0
		for i in range(len(s)):
			if 'amplicon' in s[i]:
				amplicon_idx = i
				break
		assert amplicon_idx != 0
		amplicon_id = s[amplicon_idx]
		amplicon_type = s[amplicon_idx + 1]
		if sample_name not in amplicon_dict:
			amplicon_dict[sample_name] = dict()
		if amplicon_id not in amplicon_dict[sample_name]:
			amplicon_dict[sample_name][amplicon_id] = dict()
		fp = open(sys.argv[5] + fn, 'r')
		for line in fp:
			line = line.strip()
			s = line.split()
			try:
				amplicon_dict[sample_name][amplicon_id][amplicon_type].append([s[0], int(s[1]), int(s[2])])
			except:
				amplicon_dict[sample_name][amplicon_id][amplicon_type] = [[s[0], int(s[1]), int(s[2])]]
		fp.close()
	print ("Identified %d WGS samples having valid focal amplifications." %(len(amplicon_dict)))
	
	# Make up 'No amp/Invalid' classifications by AC
	fp = open(sys.argv[6], 'r') # AC *_amplicon_classification_profiles.tsv
	for line in fp:
		line = line.strip()
		s = line.split('\t')
		if s[0] != 'sample_name':
			if s[0] not in amplicon_dict:
				amplicon_dict[s[0]] = dict()
				assert s[2] == 'No amp/Invalid'
			if s[1] not in amplicon_dict[s[0]]:
				amplicon_dict[s[0]][s[1]] = dict()
				assert s[2] == 'No amp/Invalid'
				amplicon_dict[s[0]][s[1]][s[2]] = []
			if s[2] == 'No amp/Invalid':
				fp_cycle = open(sys.argv[7] + s[0] + "_" + s[1] + "_annotated_cycles.txt", 'r')
				for line_ in fp_cycle:
					line_ = line_.strip()
					t = line_.split()
					if t[0] == 'Segment':
						overlap = 0
						for intrvli in range(len(amplicon_dict[s[0]][s[1]][s[2]])):
							intrvl = amplicon_dict[s[0]][s[1]][s[2]][intrvli]
							if t[2] == intrvl[0] and int(t[3]) <= intrvl[2] and intrvl[1] <= int(t[4]):
								amplicon_dict[s[0]][s[1]][s[2]][intrvli][1] = min(amplicon_dict[s[0]][s[1]][s[2]][intrvli][1], int(t[3]))
								amplicon_dict[s[0]][s[1]][s[2]][intrvli][2] = max(amplicon_dict[s[0]][s[1]][s[2]][intrvli][2], int(t[4]))
								overlap = 1
								break
						if overlap == 0:
							amplicon_dict[s[0]][s[1]][s[2]].append([t[2], int(t[3]), int(t[4])])
				fp_cycle.close()
	fp.close()
	print ("Identified %d WGS samples in total having focal amplifications." %(len(amplicon_dict)))
	
	# Prioritize ecDNA > BFB > Complex > Linear
	amplicon_type_idx = {'ecDNA': 0, 'BFB': 1, 'Complex non-cyclic': 2, 'Linear amplification': 3, 'unknown' : 4, 'No amp/Invalid' : 5}
	for sn in amplicon_dict.keys():
		for aid in amplicon_dict[sn].keys():
			if len(amplicon_dict[sn][aid]) > 1:
				del_list = []
				for at_i1 in range(len(amplicon_dict[sn][aid])):
					for at_i2 in range(at_i1 + 1, len(amplicon_dict[sn][aid])):
						at1, sint1 = list(amplicon_dict[sn][aid].items())[at_i1]
						at2, sint2 = list(amplicon_dict[sn][aid].items())[at_i2]
						if check_int_overlap(sint1, sint2):
							if amplicon_type_idx[at1] > amplicon_type_idx[at2]:
								del_list.append(at1)
							if amplicon_type_idx[at1] < amplicon_type_idx[at2]:
								del_list.append(at2)
				for at in del_list:
					del amplicon_dict[sn][aid][at]

	# Count number of loops
	amplicon_type_idx = {'ecDNA': 0, 'BFB': 1, 'Complex non-cyclic': 2, 'Linear amplification': 3, 'No amp/Invalid' : 4}
	fp_w = open("num_neoloops_per_mbp.tsv", 'w')
	fp_w.write("Sample\tHiChIP sample\tAmplicon\tAmplicon classification\tLength of amplified intervals\tNum loops within amplicon not in control\tTotal num loops within amplicon\tNum loops across amplicon not in control\tTotal num loops across amplicon\n")
	neoloop_dist_1 = [[], [], [], [], []] 
	neoloop_dist_2 = [[], [], [], [], []]
	for (wgs, hic) in neo_loop_dict.keys():
		sample_amplified = 0
		for fn in listdir(sys.argv[7]):
			if wgs in fn:
				sample_amplified = 1
				break
		if sample_amplified == 1:
			for amplicon_id in amplicon_dict[wgs].keys():
				for amplicon_type in amplicon_dict[wgs][amplicon_id].keys():
					if amplicon_type == 'unknown':
						continue
					total_int_len = 0
					num_loops = [0, 0, 0, 0]
					for interval in amplicon_dict[wgs][amplicon_id][amplicon_type]:
						total_int_len += (interval[2] - interval[1])
					cont_loop_reduced = []
					for cont_loop in neo_loop_dict_control[(wgs, hic)]:
						if cont_loop[-1] not in control_list[wgs]:
							continue
						try:
							rc1, rc2 = check_int(cont_loop[0], cont_loop[3], int(cont_loop[1]), int(cont_loop[2]), 
									int(cont_loop[4]), int(cont_loop[5]), amplicon_dict[wgs][amplicon_id][amplicon_type])
						except:
							rc1, rc2 = -1, -1
						if rc1 >= 0 or rc2 >= 0:
							cont_loop_reduced.append(cont_loop)
					for loop in neo_loop_dict[(wgs, hic)]:
						try:
							r1, r2 = check_int(loop[0], loop[3], int(loop[1]), int(loop[2]), 
								int(loop[4]), int(loop[5]), amplicon_dict[wgs][amplicon_id][amplicon_type])
						except:
							r1, r2 = -1, -1
						loop_in_control = search_cont(loop, cont_loop_reduced)
						if r1 >= 0 and r2 >=0:
							if not loop_in_control:
								num_loops[0] += 1
							num_loops[1] += 1
							neo_loop_ = 0
							for j in range(len(loop[-1].split(','))):
								if j % 3 == 2 and loop[-1].split(',')[j] == '1':
									neo_loop_ = 1
									break
							if neo_loop_ == 1:
								if not loop_in_control:
									num_loops[2] += 1
								num_loops[3] += 1

						"""
						elif r1 >= 0:
							if not loop_in_control:
								num_loops[2] += 1
							num_loops[3] += 1
						elif r2 >= 0:
							if not loop_in_control:
								num_loops[2] += 1
							num_loops[3] += 1
						"""
					#if total_int_len > 500000:
					neoloop_dist_1[amplicon_type_idx[amplicon_type]].append(num_loops[0] * 1000000.0 / total_int_len)
					neoloop_dist_2[amplicon_type_idx[amplicon_type]].append(num_loops[2] * 1000000.0 / total_int_len)
					fp_w.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n" %(wgs, hic, amplicon_id, amplicon_type, total_int_len, 
								num_loops[0], num_loops[1], num_loops[2], num_loops[3]))
	fp_w.close()
	print ("%d amplicons passed the length filter." %(sum([len(neoloop_dist_1[i]) for i in range(5)])))

	# Violin plots
	fig, axs = plt.subplots(1, 1)
	vp = axs.violinplot(neoloop_dist_1, [1, 3, 5, 7, 9], widths = 1.0, showmeans = False,
                    showmedians=False, showextrema = False, bw_method = "silverman")
	for vpp in vp['bodies']:
		vpp.set_facecolor('none')
		vpp.set_edgecolor('#282724')
		vpp.set_alpha(1)
		vpp.set_linewidth(2)
	medianprops = {'linewidth': 4, 'color': '#747473', 'solid_capstyle': 'butt'}
	boxprops = {'linewidth': 2, 'color': '#747473'}

	axs.boxplot(neoloop_dist_1, positions = [1, 3, 5, 7, 9], showfliers = False, showcaps = False, medianprops = medianprops,
		whiskerprops = boxprops, boxprops = boxprops)
	COLOR_SCALE = ['#4A708B', '#4BAD4A', '#F9AB64', '#698AC6', '#FCFDF3']
	jitter = 0.04
	x_data = []
	for i in range(len(neoloop_dist_1)):
		x_ = []
		for d in range(len(neoloop_dist_1[i])):
			x_.append(i * 2 + 1)
		x_data.append(x_)
	
	x_jittered = [x_ + scipy.stats.t(df = 6, scale = 0.04).rvs(len(x_)) for x_ in x_data]
	i = 0
	for x, y, color in zip(x_jittered, neoloop_dist_1, COLOR_SCALE):
		if i < 4:
			axs.scatter(x, y, s = 100, color = color, alpha=0.5)
		else:
			axs.scatter(x, y, s = 100, color = color, edgecolors = 'k', alpha=0.5)
		i += 1

	max_y = max(max(neoloop_dist_1[0]), max(neoloop_dist_1[1]), max(neoloop_dist_1[2]), max(neoloop_dist_1[3]),
			max(neoloop_dist_1[4]))
	pvals = [ranksums(neoloop_dist_1[0], neoloop_dist_1[i], 'greater') for i in range(1, len(neoloop_dist_1))]
	axs.plot([1, 1, 3, 3], [max_y * 1.05, max_y * 1.07, max_y * 1.07, max_y * 1.05], '-k', linewidth = 1)
	axs.plot([1, 1, 5, 5], [max_y * 1.13, max_y * 1.15, max_y * 1.15, max_y * 1.13], '-k', linewidth = 1)
	axs.plot([1, 1, 7, 7], [max_y * 1.21, max_y * 1.23, max_y * 1.23, max_y * 1.21], '-k', linewidth = 1)
	axs.plot([1, 1, 9, 9], [max_y * 1.29, max_y * 1.31, max_y * 1.31, max_y * 1.29], '-k', linewidth = 1)
	axs.text(2, max_y * 1.072, "P=%.2e" %pvals[0][1], fontsize=30, va="bottom", ha="center")
	axs.text(3, max_y * 1.152, "P=%.2e" %pvals[1][1], fontsize=30, va="bottom", ha="center")
	axs.text(4, max_y * 1.232, "P=%.2e" %pvals[2][1], fontsize=30, va="bottom", ha="center")
	axs.text(5, max_y * 1.312, "P=%.2e" %pvals[3][1], fontsize=30, va="bottom", ha="center")
	plt.xticks([1, 3, 5, 7, 9], ['ecDNA', 'BFB', 'Complex', 'Linear', 'No amp'])
	plt.xlim([0, 10])
	axs.set_ylabel('#Loops both ends in amplicon/Mb', fontsize = 30)

	plt.tight_layout()
	plt.savefig(sys.argv[8])
	plt.close()
	

	

	
	 
	

	
				
	
	



