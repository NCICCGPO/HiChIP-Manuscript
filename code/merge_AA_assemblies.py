"""
HiChIP analysis - Merge assemble_aa_neoloop.py output from the same HiChIP sample
 
Usage: python(3) merge_AA_assemblies.py <HiChIP to WGS map (hic2sample_all.txt)>
		<Directory containing ALL assemble_aa_neoloop.py output files (*_AA_assemblies_*.txt)>
		<Output directory>
Output: 
	*_H3K27ac_AA_assemblies.txt in output directory, for all HiChIP samples
"""
import sys
from os import listdir

if __name__ == '__main__':

	hic2sample = dict()
	fp = open(sys.argv[1], 'r')
	for line in fp:
		line = line.strip()
		s = line.split()
		try:
			hic2sample[s[0]].append(s[1])
		except:
			hic2sample[s[0]] = [s[1]]
	fp.close()
	
	fdir = sys.argv[2] + "/"
	for hic in hic2sample.keys():
		if len(hic2sample[hic]) == 1:
			for fn in listdir(fdir):
				if hic in fn:
					fp_r = open(fdir + fn, 'r')
					fp_w = open(sys.argv[3] + "/" + hic + "_H3K27ac_AA_assemblies.txt", 'w')
					for line in fp_r:
						fp_w.write("%s" %line)
					fp_r.close()
					fp_w.close()
		else:
			fp_w = open(sys.argv[3] + "/" + hic + "_H3K27ac_AA_assemblies.txt", 'w')
			aa_assemblies = []
			nlf_assemblies = []
			for fn in listdir(fdir):
				if hic in fn:
					fp_r = open(fdir + fn, 'r')
					for line in fp_r:
						line = line.strip()
						if line[0] == 'C' or line[0] == 'A':
							if line not in nlf_assemblies:
								nlf_assemblies.append(line)
						else:
							aa_assemblies.append(line)
					fp_r.close()
			for line in aa_assemblies:
				fp_w.write("%s\n" %line)
			for line in nlf_assemblies:
				fp_w.write("%s\n" %line)
			fp_w.close()



	