"""
HiChIP analysis - For each WGS sample, identify the list of all HiChIP which do not have an amplicon overlap with that WGS 
Usage: python(3) extract_control_samples.py <HiChIP to WGS map (hic2sample_all.txt)>
		<AC summary of all WGS samples (*_amplicon_classification_profiles.tsv)>
		<Directory containing ALL cycle files (can use AC generated *_annotated_cycles_files)>
		<Output filename (e.g. control_list_all.tsv)>
Output: a *.tsv file in the following format
	WGS1	CONTROL1
	WGS1	CONTROL2
	...
	WGS1	CONTROLN
	WGS2	CONTROL1
	...
	WGS2	CONTROLN
	...
	WGSN	CONTROLN
"""
import sys
from os import listdir


if __name__ == '__main__':

	hic_list = []
	sample2hic = dict()
	fp = open(sys.argv[1], 'r')
	for line in fp:
		line = line.strip()
		s = line.split()
		if s[0] not in hic_list:
			hic_list.append(s[0])
		sample2hic[s[1]] = s[0]
	fp.close()

	amplified_intervals_all = dict()
	fp = open(sys.argv[2], 'r')
	for line in fp:
		line = line.strip()
		s = line.split()
		if '_' in s[1]:
			continue
		sample = s[0]
		amplicon_id = s[1]
		aa_intervals = []
		for fn_cycle in listdir(sys.argv[3]):
			if sample in fn_cycle and amplicon_id in fn_cycle:
				fp_cycle = open(sys.argv[3] + fn_cycle, 'r')
				for line in fp_cycle:
					line = line.strip()
					tokens = line.split()
					if tokens[0] == 'Segment':
						aa_intervals.append([tokens[2], int(tokens[3]), int(tokens[4])])		
				fp_cycle.close()
		# exclude samples which do not have hichip file
		if sample in sample2hic.keys():
			try:
				amplified_intervals_all[sample] += aa_intervals
			except:
				amplified_intervals_all[sample] = aa_intervals
	fp.close()

	fp_w = open(sys.argv[4], 'w')
	for sample in sample2hic.keys():
		hic_flags = dict()
		for hic in hic_list:
			hic_flags[hic] = 1
		if sample in amplified_intervals_all:
			aa_intervals = amplified_intervals_all[sample]
			for sample_ in sample2hic.keys():
				if sample2hic[sample_] != sample2hic[sample] and sample_ in amplified_intervals_all:
					aa_intervals_ = amplified_intervals_all[sample_]
					for intrvl in aa_intervals:
						for intrvl_ in aa_intervals_:
							if intrvl[0] == intrvl_[0] and intrvl[1] <= intrvl_[2] and intrvl_[1] <= intrvl[2]:
								if hic_flags[sample2hic[sample_]] == 1:
									hic_flags[sample2hic[sample_]] = 0
									break		
						if hic_flags[sample2hic[sample_]] == 0:
							break
		for hic in hic_flags.keys():
			if hic != sample2hic[sample] and hic_flags[hic] == 1:
				fp_w.write("%s\t%s\n" %(sample, hic))
	fp_w.close()
	

	
