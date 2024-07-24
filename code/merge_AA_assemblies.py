from os import listdir

if __name__ == '__main__':

	hic2sample = dict()
	fp = open("/nucleus/projects/GDC/Significant_Loops/hic2sample_all.txt", 'r')
	for line in fp:
		line = line.strip()
		s = line.split()
		try:
			hic2sample[s[0]].append(s[1])
		except:
			hic2sample[s[0]] = [s[1]]
	fp.close()
	
	fdir = "/nucleus/projects/GDC/new_data/AA_assemblies_sep_/"
	for hic in hic2sample.keys():
		if len(hic2sample[hic]) == 1:
			for fn in listdir(fdir):
				if hic in fn:
					fp_r = open(fdir + fn, 'r')
					fp_w = open("/nucleus/projects/GDC/new_data/AA_assemblies_merged_/" + hic + "_H3K27ac_AA_assemblies.txt", 'w')
					for line in fp_r:
						fp_w.write("%s" %line)
					fp_r.close()
					fp_w.close()
		else:
			fp_w = open("/nucleus/projects/GDC/new_data/AA_assemblies_merged_/" + hic + "_H3K27ac_AA_assemblies.txt", 'w')
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



	