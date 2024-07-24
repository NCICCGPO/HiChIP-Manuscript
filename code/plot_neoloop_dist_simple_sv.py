"""
HiChIP analysis - (Neo)loop distribution on different types of simple SVs (INV/DUP/DEL/TRANSLOC)
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
import math
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
	hic2sample = dict()
	fp = open(sys.argv[1], 'r') # Path to "hic2sample_all.txt"
	for line in fp:
		line = line.strip()
		s = line.split()
		sample_hic.append([s[1], s[0]])
		try:
			hic2sample[s[0]].append(s[1])
		except:
			hic2sample[s[0]] = [s[1]]
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
	for hic in hic2sample.keys():
		for fn in listdir(sys.argv[3]): # Directory of case neoloop calls
			if hic in fn and 'neo-loops' in fn:
				fp = open(sys.argv[3] + fn, 'r')
				neo_loop_dict[hic] = []
				neo_loop_dict_control[hic] = []
				for line in fp:
					line = line.strip()
					s = line.split()
					if s[-1][0] == 'A' or s[-1][0] == 'C':
						neo_loop_dict[hic].append(s)
				fp.close()
		for cont_fn in listdir(sys.argv[4]): # Directory of control neoloop calls
			if hic in cont_fn and 'neo-loops' in cont_fn:
				fp_cont = open(sys.argv[4] + cont_fn, 'r')
				for line in fp_cont:
					line = line.strip()
					s = line.split()
					if s[-1] in control_list[hic2sample[hic][0]] and (s[-2][0] == 'A' or s[-2][0] == 'C'):
						neo_loop_dict_control[hic].append(s[:-1])
				fp_cont.close()
	print("Identified all %d samples having NeoLoop results available." %len(neo_loop_dict))

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

	# Extract simple SV assemblies not overlapping with 
	sv_dict = dict()
	for hic in hic2sample.keys():
		for fn in listdir(sys.argv[8]): # Directory of NeoloopFinder assemblies
			if hic in fn:
				if hic not in sv_dict:
					sv_dict[hic] = dict()
					fp_assemblies = open(sys.argv[8] + fn, 'r')
					for line_ in fp_assemblies:
						line_ = line_.strip()
						t = line_.split()
						if t[0][0] == 'C':
							sv_type = t[1].split(',')[0]
							sv_int_1 = ['chr' + t[1].split(',')[1], min(int(t[1].split(',')[2]), int(t[2].split(',')[1])), max(int(t[1].split(',')[2]), int(t[2].split(',')[1]))]
							sv_int_2 = ['chr' + t[3].split(',')[0], min(int(t[1].split(',')[5]), int(t[3].split(',')[1])), max(int(t[1].split(',')[5]), int(t[3].split(',')[1]))]
							overlap_flag = False
							bp = ['chr' + t[1].split(',')[1], int(t[1].split(',')[2]), t[1].split(',')[3], 'chr' + t[1].split(',')[4], int(t[1].split(',')[5]), t[1].split(',')[6]]
							for wgs in hic2sample[hic]:
								if wgs in amplicon_dict:
									for amplicon_id in amplicon_dict[wgs].keys():
										for amplicon_type in amplicon_dict[wgs][amplicon_id].keys():
											if check_int_overlap([sv_int_1, sv_int_2], amplicon_dict[wgs][amplicon_id][amplicon_type]):
												overlap_flag = True
												break
							if not overlap_flag:
								"""
								if bp[0] == sv_int_1[0] and bp[1] == sv_int_1[1]:
									sv_int_1[2] = min(sv_int_1[2], int(math.ceil(sv_int_1[1] / 10000)) * 10000 + 500000)
								else:
									sv_int_1[1] = max(sv_int_1[1], int(math.floor(sv_int_1[2] / 10000)) * 10000 - 500000)
								if bp[3] == sv_int_2[0] and bp[4] == sv_int_2[1]:
									sv_int_2[2] = min(sv_int_2[2], int(math.ceil(sv_int_2[1] / 10000)) * 10000 + 500000)
								else:
									sv_int_2[1] = max(sv_int_2[1], int(math.floor(sv_int_2[2] / 10000)) * 10000 - 500000)
								"""
								try:
									sv_dict[hic][sv_type].append((sv_int_1, sv_int_2))
								except:
									sv_dict[hic][sv_type] = [(sv_int_1, sv_int_2)]
							#else:
							#	print (sv_int_1, sv_int_2)
					fp_assemblies.close()
	print ("Identified %d HiChIP samples with NeoLoopFinder assemblies." %(len([sv_dict[hic] for hic in sv_dict.keys() if len(sv_dict[hic]) > 0 ])))
	
	neoloop_dist_1 = [[], [], [], []] # Inverison, duplication, deletion, translocation
	fp_w = open("num_neoloops_per_mbp_sv.tsv", 'w')
	fp_w.write("HiChIP sample\tSeg1\tSeg2\tSV type\tLength\tNum loops within SV assembly not in control\tTotal num loops within SV assembly\tNum neoloops within SV assembly not in control\tTotal num neoloops within SV assembly\n")
	sv_type_idx = {"inversion": 0, "duplication": 1, "deletion":2, "translocation":3}
	sv_type_idx_ = {0: "inversion", 1: "duplication", 2: "deletion", 3: "translocation"}
	for hic in sv_dict.keys():
		num_loops = []
		total_int_len = []
		sv_info = []
		for sv_type in sv_dict[hic].keys():
			for (interval1, interval2) in sv_dict[hic][sv_type]:
				total_int_len.append((interval1[2] - interval1[1]) + (interval2[2] - interval2[1]))
				num_loop = [0, 0, 0, 0]
				cont_loop_reduced = []
				for cont_loop in neo_loop_dict_control[hic]:
					try:
						rc1, rc2 = check_int(cont_loop[0], cont_loop[3], int(cont_loop[1]), int(cont_loop[2]), 
									int(cont_loop[4]), int(cont_loop[5]), [interval1, interval2])
					except:
						rc1, rc2 = -1, -1
					if rc1 >= 0 or rc2 >= 0:
						cont_loop_reduced.append(cont_loop)
				for loop in neo_loop_dict[hic]:
					try:
						r1, r2 = check_int(loop[0], loop[3], int(loop[1]), int(loop[2]), 
								int(loop[4]), int(loop[5]), [interval1, interval2])
					except:
						r1, r2 = -1, -1
					loop_in_control = search_cont(loop, cont_loop_reduced)
					if r1 >= 0 and r2 >=0:
						neo_loop_ = 0
						for j in range(len(loop[-1].split(','))):
							if j % 3 == 2 and loop[-1].split(',')[j] == '1':
								neo_loop_ = 1
								break
						if not loop_in_control:
							num_loop[0] += 1
						num_loop[1] += 1
						if neo_loop_ == 1:
							if not loop_in_control:
								num_loop[2] += 1
							num_loop[3] += 1
				num_loops.append(num_loop)
				sv_info.append((sv_type, interval1, interval2))
		for i in range(len(total_int_len)):
			neoloop_dist_1[sv_type_idx[sv_info[i][0]]].append(num_loops[i][0] * 1000000.0 / total_int_len[i])
			fp_w.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n" %(hic, sv_info[i][1], sv_info[i][2], sv_info[i][0], 
					total_int_len[i], num_loops[i][0], num_loops[i][1], num_loops[i][2], num_loops[i][3]))
	fp_w.close()

	# Violin plots
	fig, axs = plt.subplots(1, 1)
	vp = axs.violinplot(neoloop_dist_1, [1, 3, 5, 7], widths = 1.0, showmeans = False,
				showmedians=False, showextrema = False, bw_method = "silverman")
	for vpp in vp['bodies']:
		vpp.set_facecolor('none')
		vpp.set_edgecolor('#282724')
		vpp.set_alpha(1)
		vpp.set_linewidth(2)

	medianprops = {'linewidth': 4, 'color': '#747473', 'solid_capstyle': 'butt'}
	boxprops = {'linewidth': 2, 'color': '#747473'}

	axs.boxplot(neoloop_dist_1, positions = [1, 3, 5, 7], showfliers = False, showcaps = False, medianprops = medianprops,
		whiskerprops = boxprops, boxprops = boxprops)
	COLOR_SCALE = ['#5E67AA', '#9284B4', '#E9EDC9', '#CCD5AE']
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

	max_y = max(max(neoloop_dist_1[0]), max(neoloop_dist_1[1]), max(neoloop_dist_1[2]), max(neoloop_dist_1[3]))
	pvals = [ranksums(neoloop_dist_1[0], neoloop_dist_1[i], 'greater') for i in range(1, len(neoloop_dist_1))]
	#pvals = [ttest_ind(neoloop_dist_0_[0], neoloop_dist_0_[i], equal_var=False) for i in range(1, len(neoloop_dist_0_))]
	axs.plot([1, 1, 3, 3], [max_y * 1.05, max_y * 1.07, max_y * 1.07, max_y * 1.05], '-k', linewidth = 1)
	axs.plot([1, 1, 5, 5], [max_y * 1.13, max_y * 1.15, max_y * 1.15, max_y * 1.13], '-k', linewidth = 1)
	axs.plot([1, 1, 7, 7], [max_y * 1.21, max_y * 1.23, max_y * 1.23, max_y * 1.21], '-k', linewidth = 1)
	axs.text(2, max_y * 1.072, "P=%.2e" %pvals[0][1], fontsize=30, va="bottom", ha="center")
	axs.text(3, max_y * 1.152, "P=%.2e" %pvals[1][1], fontsize=30, va="bottom", ha="center")
	axs.text(4, max_y * 1.232, "P=%.2e" %pvals[2][1], fontsize=30, va="bottom", ha="center")
	plt.xticks([1, 3, 5, 7], ['inversion', 'duplication', 'deletion', 'translocation'])
	plt.xlim([0, 8])
	#axs[0, 0].set_yticks([])
	axs.set_ylabel('#Loops both ends in SV assembly/Mb', fontsize = 30)
	
	
	
	plt.tight_layout()
	#plt.savefig("neoloop_distribution.png", dpi = 600)
	plt.savefig("neoloop_distribution_sv.pdf")
	plt.close()
	

