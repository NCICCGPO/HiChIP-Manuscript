import os
import sys
from os import listdir
import math

def interval_overlap(segment, aintervals):
	if segment[0] not in aintervals:
		return False
	else:
		for ainterval in aintervals[segment[0]]:
			if segment[1] < ainterval[1] and ainterval[0] < segment[2]:
					return True
		return False


def check_foldbacks(cycle):
	foldbacks = []
	for sidx in range(-1, len(cycle) - 1):
		if cycle[sidx][0] == cycle[sidx + 1][0] and cycle[sidx][1] < cycle[sidx + 1][2] and cycle[sidx + 1][1] < cycle[sidx][2]:
			foldbacks.append(sidx) 
	return foldbacks


def neoloop_breakpoint(seg_, nextseg_):
	b_ = [] 
	if seg_[-1] == '+' and nextseg_[-1] == '-':
		b_ = ['inversion', seg_[0], str(seg_[2]), '+', nextseg_[0], str(nextseg_[2]), '+']
	elif seg_[-1] == '-' and nextseg_[-1] == '+':
		b_ = ['inversion', seg_[0], str(seg_[1]), '-', nextseg_[0], str(nextseg_[1]), '-']
	elif seg_[-1] == '+' and nextseg_[-1] == '+':
		if seg_[2] < nextseg_[1]:
			b_ = ['deletion', seg_[0], str(seg_[2]), '+', nextseg_[0], str(nextseg_[1]), '-']
		else:
			b_ = ['duplication', seg_[0], str(seg_[2]), '+', nextseg_[0], str(nextseg_[1]), '-']
	else:
		if nextseg_[2] < seg_[1]:
			b_ = ['deletion', seg_[0], str(seg_[1]), '-', nextseg_[0], str(nextseg_[2]), '+']
		else:
			b_ = ['duplication', seg_[0], str(seg_[1]), '-', nextseg_[0], str(nextseg_[2]), '+']
	if seg_[0] != nextseg_[0]:
		b_[0] = 'translocation'
	breakpoint = '' + b_[0] + ',' + b_[1] + ',' + b_[2] + ',' + b_[3] + ',' + b_[4] + ',' + b_[5] + ',' + b_[6]
	return breakpoint


# Local assembly, extend from the left end of cycle[sidx] to right until repeat
def assemble_right(cycle, sidx):
	#print 'right', sidx
	assert sidx < len(cycle)
	al, ar = -1, -1
	alchr, archr = '', ''
	aintervals = dict()
	if cycle[sidx][-1] == '+':
		al = int(math.floor(cycle[sidx][1] / 10000.0)) * 10000
	else:
		al = int(math.ceil(cycle[sidx][2] / 10000.0)) * 10000
	
	blist = []
	sidx_ = sidx
	seg = cycle[sidx_]
	alchr = seg[0]
	nextseg = cycle[(sidx_ + 1) % len(cycle)]
	try:
		aintervals[seg[0]].append([seg[1], seg[2]])
	except:
		aintervals[seg[0]] = [[seg[1], seg[2]]]
	#print aintervals, nextseg
	while not interval_overlap(nextseg, aintervals):
		try:
			aintervals[nextseg[0]].append([nextseg[1], nextseg[2]])
		except:
			aintervals[nextseg[0]] = [[nextseg[1], nextseg[2]]]
		breakpoint = neoloop_breakpoint(seg, nextseg)
		blist.append(breakpoint)
		seg = nextseg 
		sidx_ += 1
		nextseg = cycle[(sidx_ + 1) % len(cycle)]

	if len(blist) == 0:
		return ''
	else:
		archr = seg[0]
		if seg[-1] == '+':
			ar = int(math.ceil(seg[2] / 10000.0)) * 10000
		else:
			ar = int(math.floor(seg[1] / 10000.0)) * 10000
		assembly_ = ''
		for breakpoint in blist:
			assembly_ += (breakpoint + '\t')
		assembly_ends = "%s,%d\t%s,%d" %(alchr, al, archr, ar)
		return assembly_ + assembly_ends


# Local assembly, extend from the right end of cycle[sidx] to left until repeat
# Note that traverse is still from left to right
def assemble_left(cycle, sidx):
	#print 'left', sidx
	assert sidx < len(cycle)
	al, ar = -1, -1
	alchr, archr = '', ''
	aintervals = dict()
	if cycle[sidx][-1] == '+':
		ar = int(math.ceil(cycle[sidx][2] / 10000.0)) * 10000
	else:
		ar = int(math.floor(cycle[sidx][1] / 10000.0)) * 10000
	
	blist = []
	sidx_ = sidx
	seg = cycle[sidx_]
	archr = seg[0]
	nextseg = cycle[sidx_ - 1]
	try:
		aintervals[seg[0]].append([seg[1], seg[2]])
	except:
		aintervals[seg[0]] = [[seg[1], seg[2]]]
	while not interval_overlap(nextseg, aintervals):
		try:
			aintervals[nextseg[0]].append([nextseg[1], nextseg[2]])
		except:
			aintervals[nextseg[0]] = [[nextseg[1], nextseg[2]]]
		breakpoint = neoloop_breakpoint(nextseg, seg)
		blist.append(breakpoint)
		seg = nextseg 
		sidx_ -= 1
		nextseg = cycle[sidx_ - 1]

	if len(blist) == 0:
		return ''
	else:
		alchr = seg[0]
		if seg[-1] == '+':
			al = int(math.floor(seg[1] / 10000.0)) * 10000
		else:
			al = int(math.ceil(seg[2] / 10000.0)) * 10000
		assembly_ = ''
		for breakpoint in blist[::-1]:
			assembly_ += (breakpoint + '\t')
		assembly_ends = "%s,%d\t%s,%d" %(alchr, al, archr, ar)
		return assembly_ + assembly_ends


# Local assembly, break the sequence at cycle[sidx]
# Potential caveat: will not check if two segments overlap 
def assemble_split_seg(cycle, sidx):
	assert sidx < len(cycle)
	mid_l = int(math.floor((cycle[sidx][1] + cycle[sidx][2]) / 20000.0)) * 10000
	mid_r = int(math.ceil((cycle[sidx][1] + cycle[sidx][2]) / 20000.0)) * 10000
	if cycle[sidx][-1] == '+':
		mid_l, mid_r = mid_r, mid_l	
	
	blist = []
	num_break_visited = 0
	sidx_ = sidx					
	while num_break_visited < len(cycle):
		seg = cycle[(sidx_) % len(cycle)]
		nextseg = cycle[(sidx_ + 1) % len(cycle)]
		breakpoint = neoloop_breakpoint(seg, nextseg)
		blist.append(breakpoint)
		sidx_ += 1
		num_break_visited += 1
	
	if len(blist) == 0:
		return ''
	assembly_ = ''
	for breakpoint in blist:
		assembly_ += (breakpoint + '\t')
	assembly_ends = "%s,%d\t%s,%d" %(cycle[sidx][0], mid_l, cycle[sidx][0], mid_r)
	return (assembly_ + assembly_ends)


# Local assembly, from the left end of cycle[sidx1] to the right end of cycle[sidx2]
def assemble_twoidx(cycle, sidx1, sidx2):
	assert sidx1 < len(cycle)
	assert sidx2 < len(cycle)
	if sidx1 == sidx2:
		return ''
	al, ar = -1, -1
	alchr, archr = '', ''
	aintervals = dict()
	if cycle[sidx1][-1] == '+':
		al = int(math.floor(cycle[sidx1][1] / 10000.0)) * 10000
	else:
		al = int(math.ceil(cycle[sidx1][2] / 10000.0)) * 10000

	blist = []
	sidx_ = sidx1
	seg = cycle[sidx_]
	alchr = seg[0]
	nextseg = cycle[(sidx_ + 1) % len(cycle)]
	try:
		aintervals[seg[0]].append([seg[1], seg[2]])
	except:
		aintervals[seg[0]] = [[seg[1], seg[2]]]
	if sidx1 > sidx2:
		sidx2_ = sidx2 + len(cycle)
	else:
		sidx2_ = sidx2
	while (not interval_overlap(nextseg, aintervals)) and sidx_ <= sidx2_:
		try:
			aintervals[nextseg[0]].append([nextseg[1], nextseg[2]])
		except:
			aintervals[nextseg[0]] = [[nextseg[1], nextseg[2]]]
		breakpoint = neoloop_breakpoint(seg, nextseg)
		blist.append(breakpoint)
		seg = nextseg 
		sidx_ += 1
		nextseg = cycle[(sidx_ + 1) % len(cycle)]

	if len(blist) == 0:
		return ''
	else:
		archr = seg[0]
		if seg[-1] == '+':
			ar = int(math.ceil(seg[2] / 10000.0)) * 10000
		else:
			ar = int(math.floor(seg[1] / 10000.0)) * 10000
		assembly_ = ''
		for breakpoint in blist:
			assembly_ += (breakpoint + '\t')
		assembly_ends = "%s,%d\t%s,%d" %(alchr, al, archr, ar)
		return assembly_ + assembly_ends


if __name__ == '__main__':

	# Read in AA cycle file
	aa_cycles = dict()
	cycle_classifications = dict()
	aa_paths = dict()
	path_classifications = dict()
	aa_segments = dict()
	cycle_flags = dict()
	path_flags = dict()

	print ("Cycle files:")
	for fn in listdir(sys.argv[1]):
		if 'cycles' in fn and sys.argv[3] in fn:
			print (fn)
			amplicon_id = fn[fn.find("amplicon"): fn.find("annotated_cycles") - 1]
			fp = open(sys.argv[1] + '/' + fn, 'r')
			cycles = [[], []] # Cycles / Paths
			cycle_ids = [[], []] # Cycles / Paths
			cycle_class = [[], []]
			flags = [[], []]
			segments = []
			for line in fp:
				line = line.strip()
				t = line.split()
				if t[0] == 'Segment':
					chr = t[2][3:]
					s_ = int(t[3])
					e_ = int(t[4])
					segments.append([chr, s_, e_])
				if t[0][:5] == 'Cycle':
					if "CycleClass" not in t[0]:
						print ("Please input AmpliconClassifier modified cycle file.")
						os.abort()
					t1 = t[0].split('=')
					t1s = t1[-1].split(',')
					t2 = t[0].split(';')[-2].split('=')[-1]
					#print t1
					if "Invalid" not in t[0]:
						if t1s[0] != '0+' and t1s[-1] != '0+' and t1[-1] not in cycles[0]:
							cycles[0].append(t1[-1])
							cycle_ids[0].append('Cycle' + t1[1].split(';')[0])
							flags[0].append(1)
							cycle_class[0].append(t2)
						if t1s[0] == '0+' and t1s[-1] == '0+' and t1[-1] not in cycles[1]:
							cycles[1].append(t1[-1])
							cycle_ids[1].append('Cycle' + t1[1].split(';')[0])
							flags[1].append(1)
							cycle_class[1].append(t2)
									
			fp.close()
			aa_segments[amplicon_id] = segments
			if amplicon_id not in aa_cycles:
				aa_cycles[amplicon_id] = dict()
				cycle_classifications[amplicon_id] = dict()
				cycle_flags[amplicon_id] = dict()
				aa_paths[amplicon_id] = dict()
				path_classifications[amplicon_id] = dict()
				path_flags[amplicon_id] = dict()
			for cycle_id in cycle_ids[0]:
				aa_cycles[amplicon_id][cycle_id] = cycles[0][cycle_ids[0].index(cycle_id)]
				cycle_classifications[amplicon_id][cycle_id] = cycle_class[0][cycle_ids[0].index(cycle_id)]
				cycle_flags[amplicon_id][cycle_id] = flags[0][cycle_ids[0].index(cycle_id)]
			for cycle_id in cycle_ids[1]:
				aa_paths[amplicon_id][cycle_id] = cycles[1][cycle_ids[1].index(cycle_id)]
				path_classifications[amplicon_id][cycle_id] = cycle_class[1][cycle_ids[1].index(cycle_id)]
				path_flags[amplicon_id][cycle_id] = flags[1][cycle_ids[1].index(cycle_id)]

	# Extract actual segments
	# Cycle = [[chr, start, end, strand] ... ]
	for amplicon_id in aa_cycles.keys():
		for cycle_id in aa_cycles[amplicon_id].keys():
			cycle = aa_cycles[amplicon_id][cycle_id]
			cycle_ = []
			for seg in cycle.split(','):
				cycle_.append(aa_segments[amplicon_id][int(seg[:-1]) - 1] + [seg[-1]])
			aa_cycles[amplicon_id][cycle_id] = cycle_
	for amplicon_id in aa_paths.keys():
		for path_id in aa_paths[amplicon_id].keys():
			path = aa_paths[amplicon_id][path_id]
			path_ = []
			for seg in path.split(','):
				if seg != '0+':
					path_.append(aa_segments[amplicon_id][int(seg[:-1]) - 1] + [seg[-1]])
			aa_paths[amplicon_id][path_id] = path_

	# Initial cycle check
	# Filter out cases (i) single sequence, length < 50000, and (ii) 2 fold-back sequences
	# Other cases left undetermined 
	for amplicon_id in aa_cycles.keys():
		for cycle_id in aa_cycles[amplicon_id].keys():
			cycle = aa_cycles[amplicon_id][cycle_id]
			if len(cycle) == 1:
				if cycle[0][2] - cycle[0][1] < 50000:
					cycle_flags[amplicon_id][cycle_id] = 0
			elif len(cycle) == 2:
				if cycle[0][0] == cycle[1][0] and cycle[0][1] < cycle[1][2] and cycle[1][1] < cycle[0][2]:
					cycle_flags[amplicon_id][cycle_id] = 0
			elif len(cycle) >= 3:
				cycle_flags[amplicon_id][cycle_id] = 2
	for amplicon_id in aa_paths.keys():
		for path_id in aa_paths[amplicon_id].keys():
			path = aa_paths[amplicon_id][path_id]
			if len(path) == 1:
				path_flags[amplicon_id][path_id] = 0
			elif len(path) == 2:
				if path[0][0] == path[1][0] and path[0][1] < path[1][2] and path[1][1] < path[0][2]:
					path_flags[amplicon_id][path_id] = 0
			elif len(path) >= 3:
				path_flags[amplicon_id][path_id] = 2

	# Read in NeoLoop breakpoints
	old_assemblies = [[], []]
	fp = open(sys.argv[2], 'r')
	for line in fp:
		line = line.strip()
		#tokens = line.split()
		if line[0] == 'A':
			old_assemblies[0].append(line)
		if line[0] == 'C':
			old_assemblies[1].append(line)
	fp.close()
	
	# Produce AA assemblies
	aa_cycle_assemblies = []
	aa_path_assemblies = []
	for amplicon_id in aa_cycles.keys():
		for cycle_id in aa_cycles[amplicon_id].keys():
			cycle = aa_cycles[amplicon_id][cycle_id]
			print ("Cycle = %s; Processing flag = %d" %(cycle, cycle_flags[amplicon_id][cycle_id]))
			if cycle_flags[amplicon_id][cycle_id] == 2:
				foldbacks = check_foldbacks(cycle)
				#print cycle, foldbacks
				if len(foldbacks) == 2:
					a1 = assemble_twoidx(cycle, foldbacks[0] + 1, foldbacks[1])
					a2 = assemble_twoidx(cycle, foldbacks[1] + 1, foldbacks[0])
					print("New assembly: %s" %a1)
					print ("New assembly: %s" %a2)
					if len(a1) > 0:
						aa_cycle_assemblies.append([amplicon_id, cycle_id, cycle_classifications[amplicon_id][cycle_id], a1])
					if len(a2) > 0:
						aa_cycle_assemblies.append([amplicon_id, cycle_id, cycle_classifications[amplicon_id][cycle_id], a2])						
				elif len(foldbacks) == 1:
					a1 = assemble_right(cycle, foldbacks[0] + 1)
					a2 = assemble_left(cycle, foldbacks[0])
					print("New assembly: %s" %a1)
					print ("New assembly: %s" %a2)
					if len(a1) > 0:
						aa_cycle_assemblies.append([amplicon_id, cycle_id, cycle_classifications[amplicon_id][cycle_id], a1])
					if len(a2) > 0:
						aa_cycle_assemblies.append([amplicon_id, cycle_id, cycle_classifications[amplicon_id][cycle_id], a2])
				elif len(foldbacks) == 0:
					for sidx in range(len(cycle)):
						if cycle[sidx][2] - cycle[sidx][1] >= 50000:
							#a1 = assemble_split_seg(cycle, sidx)
							a1 = assemble_right(cycle, sidx)
							print("New assembly: %s" %a1)
							if len(a1) > 0:
								aa_cycle_assemblies.append([amplicon_id, cycle_id, cycle_classifications[amplicon_id][cycle_id], a1]) 
						else:
							print("New assembly: ")
			elif cycle_flags[amplicon_id][cycle_id] == 1:
				#print cycle
				for sidx in range(len(cycle)):
					if cycle[sidx][2] - cycle[sidx][1] >= 50000:
						#a1 = assemble_split_seg(cycle, sidx)
						a1 = assemble_right(cycle, sidx)
						print("New assembly: %s" %a1)
						if len(a1) > 0:
							aa_cycle_assemblies.append([amplicon_id, cycle_id, cycle_classifications[amplicon_id][cycle_id], a1]) 
					else:
						print("New assembly: ")
	for amplicon_id in aa_paths.keys():
		for path_id in aa_paths[amplicon_id].keys():
			path = aa_paths[amplicon_id][path_id]
			print ("Path = %s; Processing flag = %d" %(path, path_flags[amplicon_id][path_id]))
			if path_flags[amplicon_id][path_id] == 2:
				a1 = assemble_right(path, 0)
				a2 = assemble_left(path, len(path) - 1)
				print("New assembly: %s" %a1)
				if a1 != a2:
					print("New assembly: %s" %a2)
				if len(a1) > 0:
					aa_path_assemblies.append([amplicon_id, path_id, path_classifications[amplicon_id][path_id], a1])
				if len(a2) > 0 and a1 != a2:
					aa_path_assemblies.append([amplicon_id, path_id, path_classifications[amplicon_id][path_id], a2])
			elif path_flags[amplicon_id][path_id] == 1:
				assert len(path) == 2
				a1 = assemble_twoidx(path, 0, 1)
				if len(a1) > 0:
					print("New assembly: %s" %a1)
					aa_path_assemblies.append([amplicon_id, path_id, path_classifications[amplicon_id][path_id], a1]) 
				else:
					print("New assembly: ")


	fp_w = open(sys.argv[4], 'w')
	for i in range(len(aa_cycle_assemblies)):
		if aa_cycle_assemblies[i][2] == 'ecDNA-like':
			fp_w.write(sys.argv[3] + "_%s_%s_ecDNA%d\t%s\n" 
					%(aa_cycle_assemblies[i][0], aa_cycle_assemblies[i][1], i, aa_cycle_assemblies[i][3]))
		elif aa_cycle_assemblies[i][2] == 'Linear':
			fp_w.write(sys.argv[3] + "_%s_%s_Linear%d\t%s\n" 
					%(aa_cycle_assemblies[i][0], aa_cycle_assemblies[i][1], i, aa_cycle_assemblies[i][3]))
		elif aa_cycle_assemblies[i][2] == 'Rearranged':
			fp_w.write(sys.argv[3] + "_%s_%s_Rearranged%d\t%s\n" 
					%(aa_cycle_assemblies[i][0], aa_cycle_assemblies[i][1], i, aa_cycle_assemblies[i][3]))
		elif aa_cycle_assemblies[i][2] == 'BFB-like':
			fp_w.write(sys.argv[3] + "_%s_%s_BFB%d\t%s\n" 
					%(aa_cycle_assemblies[i][0], aa_cycle_assemblies[i][1], i, aa_cycle_assemblies[i][3]))
		else:
			os.abort()
	for i in range(len(aa_path_assemblies)):
		if aa_path_assemblies[i][2] == 'ecDNA-like':
			fp_w.write(sys.argv[3] + "_%s_%s_ecDNA%d\t%s\n" 
					%(aa_path_assemblies[i][0], aa_path_assemblies[i][1], i + len(aa_cycle_assemblies), aa_path_assemblies[i][3]))
		elif aa_path_assemblies[i][2] == 'Linear':
			fp_w.write(sys.argv[3] + "_%s_%s_Linear%d\t%s\n" 
					%(aa_path_assemblies[i][0], aa_path_assemblies[i][1], i + len(aa_cycle_assemblies), aa_path_assemblies[i][3]))
		elif aa_path_assemblies[i][2] == 'Rearranged':
			fp_w.write(sys.argv[3] + "_%s_%s_Rearranged%d\t%s\n" 
					%(aa_path_assemblies[i][0], aa_path_assemblies[i][1], i + len(aa_cycle_assemblies), aa_path_assemblies[i][3]))
		elif aa_path_assemblies[i][2] == 'BFB-like':
			fp_w.write(sys.argv[3] + "_%s_%s_BFB%d\t%s\n" 
					%(aa_path_assemblies[i][0], aa_path_assemblies[i][1], i + len(aa_cycle_assemblies), aa_path_assemblies[i][3]))
		else:
			os.abort()
	for assembly in old_assemblies[0]:
		fp_w.write("%s\n" %assembly)
	for assembly in old_assemblies[1]:
		fp_w.write("%s\n" %assembly)
	fp_w.close()
	#print aa_assemblies




	
	