### some genes do not have stop and start codons. MOST are on short scafs, but some seem like just not annotated.
### don't look like they have UTR annot
### have a look at Tps_LRv5b_scf3	AUGUSTUS	gene	6247362	6500000	.	+	.	ID=Tps_017089;Alias=g_7918;Name=Tps_017089
### why does it end early?


import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections


try:
	opts, args = getopt.getopt(sys.argv[1:], 'g:u:m:n:Gh')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


orig_gff_file_name = None
utr_gff_file_name  = None
max_utr_gap = 500
GAP_OUT = False
merge_n = 0

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** UTR_GUSHR_fix.py | Written by DJP, 12/06/22 in Python 3.5 in Porto ****\n")
		print("This script fixes the UTR from the Braker/GUSHR pipeline")
		print("Basically takes the original gff (no UTR) and adds UTR regions (trimmed to a max gap size) to it from the gff with UTR")
		print("This get around the issue of some genes have HUGE gaps in the UTR covering half a chromosome")		
		print("\n**** Usage ****\n")
		print("\n*********ALL gffs must be run through fix_braker_gtf.py with -G option before running this*********\n")
		print("python3 UTR_GUSHR_fix.py -g [original gff] -u [utr gff] [options] \n\n")
		print("\n**** Options ****\n")
		print("-m\tmax UTR gap size. default 500")
		print("-n\tmerge UTR gaps - merges gaps in UTR annotation. default 0")
		print("-G\tOutput a file with all UTR gap sizes. default off.")
		
		sys.exit(2)
	elif opt in ('-g'):
		orig_gff_file_name = arg
	elif opt in ('-u'):
		utr_gff_file_name  = arg
	elif opt in ('-m'):
		max_utr_gap  = int(arg)
	elif opt in ('-n'):
		merge_n = int(arg)
	elif opt in ('-G'):
		GAP_OUT = True
	else:
		print("i dont know")
		sys.exit(2)


####################################
### start and stop codons for mRNAs from orig file

all_mRNA_IDs = set()
all_mRNA_IDs_with_start_codon = set()
all_mRNA_IDs_with_stop_codon  = set()
mRNA_start_codon_dict = {}
mRNA_stop_codon_dict  = {}
mRNA_scaf_dict = {}


orig_gff_file = open(orig_gff_file_name)
for line in orig_gff_file:
	if not line.startswith("#"):
		line = line.rstrip("\n").split("\t")
		feature = line[2]
		ID = line[8].split("ID=")[1].split(";")[0].strip()
		loc = (line[0], line[3], line[4])
		strand_c = line[6]
		if feature == "mRNA":
			all_mRNA_IDs.add(ID)
			mRNA_scaf_dict[ID] = line[0]
		if feature == "start_codon":
			Parent = line[8].split("Parent=")[1].split(";")[0].strip()
			all_mRNA_IDs_with_start_codon.add(Parent)
			mRNA_start_codon_dict[Parent] = loc
		if feature == "stop_codon":
			Parent = line[8].split("Parent=")[1].split(";")[0].strip()
			all_mRNA_IDs_with_stop_codon.add(Parent)
			mRNA_stop_codon_dict[Parent] = loc

	
orig_gff_file.close()		

print("N mRNA in orig gff = " + str(len(all_mRNA_IDs)))
print("N mRNA in orig gff with start codon = " + str(len(all_mRNA_IDs_with_start_codon)))
print("N mRNA in orig gff with stop  codon = " + str(len(all_mRNA_IDs_with_stop_codon)))

# a4 = all_mRNA_IDs.difference(all_mRNA_IDs_with_start_codon) ## in a1 not in a2
# a5 = all_mRNA_IDs.difference(all_mRNA_IDs_with_stop_codon) ## in a1 not in a2
# # print(len(a4))
# # print(len(a5))

# for el in list(a4) + list(a5):
# 	print(el)
# 	scaf = mRNA_scaf_dict.get(el)
# 	print(scaf)

set_list = [all_mRNA_IDs_with_start_codon,all_mRNA_IDs_with_stop_codon]
mRNA_with_start_and_stop = set.intersection(*set_list)	

print("N mRNA in orig gff with start and stop codons = " + str(len(mRNA_with_start_and_stop)))

mRNA_pos_to_ID_dict = {}
all_mRNA_pos = set()
for m in mRNA_with_start_and_stop:
	start_pos = mRNA_start_codon_dict.get(m)
	stop_pos  = mRNA_stop_codon_dict.get(m)
	mRNA_pos_to_ID_dict[m] = (start_pos[0], start_pos[1], start_pos[2], stop_pos[1], stop_pos[2])
	all_mRNA_pos.add((start_pos[0], start_pos[1], start_pos[2], stop_pos[1], stop_pos[2]))
#print(mRNA_pos_to_ID_dict)

#####################################################################################################
### get UTRs from UTR file

#all_mRNA_loc = set()
five_prime_mRNAs_all = set()
three_prime_mRNAs_all = set()
#mRNA_loc_dict = {}
five_UTR_dict = {}
three_UTR_dict = {}


utr_gff_file = open(utr_gff_file_name)
for line in utr_gff_file:
	if not line.startswith("#"):
		line = line.rstrip("\n").split("\t")
		feature = line[2]
		ID = line[8].split("ID=")[1].split(";")[0].strip()
		
		loc = (line[0], int(line[3]), int(line[4]))
		if feature == "5'-UTR":
			Parent = line[8].split("Parent=")[1].split(";")[0].strip()
			if Parent not in five_prime_mRNAs_all:
				five_prime_mRNAs_all.add(Parent)
				five_UTR_dict[Parent] = [loc]
			else:
				rec = five_UTR_dict.get(Parent)
				rec.append(loc)
				five_UTR_dict[Parent] = rec
			
		if feature == "3'-UTR":
			Parent = line[8].split("Parent=")[1].split(";")[0].strip()
			if Parent not in three_prime_mRNAs_all:
				three_prime_mRNAs_all.add(Parent)
				three_UTR_dict[Parent] = [loc]
			else:
				rec = three_UTR_dict.get(Parent)
				rec.append(loc)
				three_UTR_dict[Parent] = rec

utr_gff_file.close()

# print(mRNA_loc_dict )
# print(five_UTR_dict)
# print(three_UTR_dict)


#############################################################################################################
### start and stop codons for mRNAs from UTR file

UTR_all_mRNA_IDs = set()
UTR_all_mRNA_IDs_with_start_codon = set()
UTR_all_mRNA_IDs_with_stop_codon  = set()
UTR_mRNA_start_codon_dict = {}
UTR_mRNA_stop_codon_dict  = {}
UTR_mRNA_scaf_dict = {}

utr_gff_file = open(utr_gff_file_name)
for line in utr_gff_file:
	if not line.startswith("#"):
		line = line.rstrip("\n").split("\t")
		feature = line[2]
		ID = line[8].split("ID=")[1].split(";")[0].strip()
		loc = (line[0], line[3], line[4])
		strand_c = line[6]
		if feature == "mRNA":
			UTR_all_mRNA_IDs.add(ID)
			UTR_mRNA_scaf_dict[ID] = line[0]
		if feature == "start_codon":
			Parent = line[8].split("Parent=")[1].split(";")[0].strip()
			UTR_all_mRNA_IDs_with_start_codon.add(Parent)
			UTR_mRNA_start_codon_dict[Parent] = loc
		if feature == "stop_codon":
			Parent = line[8].split("Parent=")[1].split(";")[0].strip()
			UTR_all_mRNA_IDs_with_stop_codon.add(Parent)
			UTR_mRNA_stop_codon_dict[Parent] = loc

	
utr_gff_file.close()		


# 
# print(len(UTR_all_mRNA_IDs))
# print(len(UTR_all_mRNA_IDs_with_start_codon))
# print(len(UTR_all_mRNA_IDs_with_stop_codon))

# a4 = UTR_all_mRNA_IDs.difference(UTR_all_mRNA_IDs_with_start_codon) ## in a1 not in a2
# a5 = UTR_all_mRNA_IDs.difference(UTR_all_mRNA_IDs_with_stop_codon) ## in a1 not in a2
# print(len(a4))
# print(len(a5))

# for el in list(a4) + list(a5):
# 	print(el)
# 	scaf = UTR_mRNA_scaf_dict.get(el)
# 	print(scaf)

set_list = [UTR_all_mRNA_IDs_with_start_codon,UTR_all_mRNA_IDs_with_stop_codon]
UTR_mRNA_with_start_and_stop = set.intersection(*set_list)	
# print(len(UTR_mRNA_with_start_and_stop))


#### get IDs with some UTR annot

all_mRNA_with_UTR = three_prime_mRNAs_all | five_prime_mRNAs_all 

UTR_mRNA_pos_to_ID_dict = {}
all_UTR_mRNA_pos = set()

for m in UTR_mRNA_with_start_and_stop:
	start_pos = UTR_mRNA_start_codon_dict.get(m)
	stop_pos  = UTR_mRNA_stop_codon_dict.get(m)
	mRNA_UTR_pos = (start_pos[0], start_pos[1], start_pos[2], stop_pos[1], stop_pos[2])
	if m in all_mRNA_with_UTR:
		if mRNA_UTR_pos not in all_UTR_mRNA_pos:
			all_UTR_mRNA_pos.add(mRNA_UTR_pos)
			UTR_mRNA_pos_to_ID_dict[mRNA_UTR_pos] = [m]
		else:
			rec = UTR_mRNA_pos_to_ID_dict.get(mRNA_UTR_pos)
			rec.append(m)
			UTR_mRNA_pos_to_ID_dict[mRNA_UTR_pos] = rec


#############################################################
### UTR_clip

def merge_small_gaps_s2b(keep_utr_list_s):
	kdict = {}
	fixed_UTR = []
	if len(keep_utr_list_s) == 1:
		fixed_UTR = [keep_utr_list_s[0]]
	else:
		merge_key = []
		for i in range(1, len(keep_utr_list_s)):
			diff_to_next = int(keep_utr_list_s[i][1]) - int(keep_utr_list_s[i-1][2])
			if diff_to_next <= merge_n :
				merge_key.append("m")
			else:
				merge_key.append("k")
		merge_key.append(merge_key[-1])
		
		k_n = 0
		m_n = 0
		seen_k = set()
		dict_order = []
		for s in range(0, len(merge_key)):
			#print(keep_utr_list_s[s])
			if merge_key[s] == "k":
				k_n = k_n + 1
				kdict["k" + str(k_n)] = [keep_utr_list_s[s]]
				dict_order.append("k" + str(k_n))
				m_n = m_n + 1
			if merge_key[s] == "m":
				if "m" + str(m_n) not in seen_k:
					dict_order.append("m" + str(m_n))
					kdict["m" + str(m_n)] = [keep_utr_list_s[s]]
					seen_k.add("m" + str(m_n))
				else:
					rec = kdict.get("m" + str(m_n))
					rec.append(keep_utr_list_s[s])
					kdict["m" + str(m_n)] = rec
		
		# print(merge_key)
		# print(kdict)
		# print(dict_order)

		for el in dict_order:
			rec = kdict.get(el)
			if el.startswith("k"):
				fixed_UTR.append(rec[0]) 
			else:
				all_coods = []
				scaf = None
				for r in rec:
					scaf = r[0] 
					all_coods.append(r[1])
					all_coods.append(r[2])
				
					
				fixed_UTR.append((scaf, min(all_coods), max(all_coods))) 	
				#print(all_coods)	
	
	return(fixed_UTR)


def merge_small_gaps_b2s(keep_utr_list_s):
	kdict = {}
	fixed_UTR = []
	#print(keep_utr_list_s)
	if len(keep_utr_list_s) == 1:
		fixed_UTR = [keep_utr_list_s[0]]
	else:
		merge_key = []
		for i in range(1, len(keep_utr_list_s)):
			diff_to_next = int(keep_utr_list_s[i-1][1]) - int(keep_utr_list_s[i][2])
			if diff_to_next <= merge_n :
				merge_key.append("m")
			else:
				merge_key.append("k")
		merge_key.append(merge_key[-1])
		
		k_n = 0
		m_n = 0
		seen_k = set()
		dict_order = []
		for s in range(0, len(merge_key)):
			#print(keep_utr_list_s[s])
			if merge_key[s] == "k":
				k_n = k_n + 1
				kdict["k" + str(k_n)] = [keep_utr_list_s[s]]
				dict_order.append("k" + str(k_n))
				m_n = m_n + 1
			if merge_key[s] == "m":
				if "m" + str(m_n) not in seen_k:
					dict_order.append("m" + str(m_n))
					kdict["m" + str(m_n)] = [keep_utr_list_s[s]]
					seen_k.add("m" + str(m_n))
				else:
					rec = kdict.get("m" + str(m_n))
					rec.append(keep_utr_list_s[s])
					kdict["m" + str(m_n)] = rec
		
		# print(merge_key)
		# print(kdict)
		# print(dict_order)

		for el in dict_order:
			rec = kdict.get(el)
			if el.startswith("k"):
				fixed_UTR.append(rec[0]) 
			else:
				all_coods = []
				scaf = None
				for r in rec:
					scaf = r[0] 
					all_coods.append(r[1])
					all_coods.append(r[2])
				
					
				fixed_UTR.append((scaf, min(all_coods), max(all_coods))) 	
				#print(all_coods)	
	return(fixed_UTR)




def three_UTR_clip(UTR_list, strand):
	fixed_UTR = []
	#print(UTR_list)
	if len(UTR_list) == 1:
		fixed_UTR.append(UTR_list[0])
	else:
		fixed_UTR = []
		if strand == "+":
			UTR_list_s = sorted(UTR_list, key=lambda x: x[1])
			#print(UTR_list_s)
			fixed_UTR = [UTR_list_s[0]]
			for i in range(1, len(UTR_list_s)):
				diff_to_next = int(UTR_list_s[i][1]) - int(UTR_list_s[i-1][2])
				#print(diff_to_next)
				if diff_to_next <= max_utr_gap:
					fixed_UTR.append(UTR_list_s[i])
				else:
					break

		elif strand == "-":
			UTR_list_s = sorted(UTR_list, key=lambda x: x[1], reverse = True)
			#print(UTR_list_s)
			fixed_UTR = [UTR_list_s[0]]
			for i in range(1, len(UTR_list_s)):
				diff_to_next =  int(UTR_list_s[i-1][1])  - int(UTR_list_s[i][2])
				#print(diff_to_next)
				if diff_to_next <= max_utr_gap:
					fixed_UTR.append(UTR_list_s[i])
				else:
				 	break
		
		else:
			print("No strand info")
			sys.exit(2)

	if strand == "+":
		fixed_UTR = merge_small_gaps_s2b(fixed_UTR)
	
	if strand == "-":
		fixed_UTR = merge_small_gaps_b2s(fixed_UTR)				
	return(fixed_UTR)

# print("ggg")
# print(three_UTR_clip(three_UTR_dict.get("Tps_015122-RB"), "+"))
# print("222")
# print(three_UTR_clip(three_UTR_dict.get("Tps_015221-RB"), "-"))
# print("xxx")
# print(three_UTR_clip(three_UTR_dict.get("Tps_000001-RB"), "+"))

def five_UTR_clip(UTR_list, strand):
	fixed_UTR = []
	#print(UTR_list)
	if len(UTR_list) == 1:
		fixed_UTR.append(UTR_list[0])
	else:
		if strand == "-":
			UTR_list_s = sorted(UTR_list, key=lambda x: x[1])
			#print(UTR_list_s)
			fixed_UTR = [UTR_list_s[0]]
			for i in range(1, len(UTR_list_s)):
				diff_to_next = int(UTR_list_s[i][1]) - int(UTR_list_s[i-1][2])
				#print(diff_to_next)
				if diff_to_next <= max_utr_gap:
					fixed_UTR.append(UTR_list_s[i])
				else:
					break
		elif strand == "+":
			UTR_list_s = sorted(UTR_list, key=lambda x: x[1], reverse = True)
			#print(UTR_list_s)
			fixed_UTR = [UTR_list_s[0]]
			for i in range(1, len(UTR_list_s)):
				diff_to_next =  int(UTR_list_s[i-1][1])  - int(UTR_list_s[i][2])
				#print(diff_to_next)
				if diff_to_next <= max_utr_gap:
					fixed_UTR.append(UTR_list_s[i])
				else:
				 	break
		else:
			print("No strand info")
			sys.exit(2)
	
	## remove small gaps
	
	#print(fixed_UTR)
	
	if strand == "-":
		fixed_UTR = merge_small_gaps_s2b(fixed_UTR)
	
	if strand == "+":
		fixed_UTR = merge_small_gaps_b2s(fixed_UTR)			
	
	return(fixed_UTR)


# print("pp999")
# print(five_UTR_clip(five_UTR_dict.get("Tps_015275-RB"), "+"))
# print("QQ999")
# print(five_UTR_clip(five_UTR_dict.get("Tps_015221-RB"), "-"))
# print("xxx")
# print(five_UTR_clip(five_UTR_dict.get("Tps_000001-RB"), "+"))
# 


# # Tps_LRv5b_scf2	AUGUSTUS	mRNA	150502225	150568970	.	+	.	ID=Tps_015275-RB;Alias=anno2.Tps_LRv5b_scf2+_file_1_g30576.t2;Parent=Tps_015275;Name=Tps_015275-RB;gene_id=Tps_015275
# # Tps_LRv5b_scf2	AnnotationFinalizer	5'-UTR	150502225	150502584	.	+	.	ID=Tps_015275-RB:5'-UTR;Parent=Tps_015275-RB;gene_id=Tps_015275
# # Tps_LRv5b_scf2	AnnotationFinalizer	5'-UTR	150502621	150506304	.	+	.	ID=Tps_015275-RB:5'-UTR;Parent=Tps_015275-RB;gene_id=Tps_015275
# 
#

def merge_ranges(my_ranges):
	ranges_a = []
	scaf     = my_ranges[0][0] 
	for i in my_ranges:
		#print(i)
		ranges_a.append((i[1], i[2]))
	#print(ranges_a)
	b = []
	for begin,end in sorted(ranges_a ):
		if b and b[-1][1] >= begin - 1:
			b[-1][1] = max(b[-1][1], end)
		else:
			b.append([begin, end])
	
	#print(b)
	out_ranges = []
	for s in b:
		out_ranges.append((scaf, s[0], s[1]))
	return(out_ranges)




if GAP_OUT == True:

	gap_file = open(orig_gff_file_name.replace(".gff", "_gaps" + ".txt"), "w")
	gap_file.write("mRNA_id\tutrgap\n")
	every_gap = []
	
	def three_UTR_gaps(UTR_list, strand):
		gaps = []
		if len(UTR_list) == 1:
			gaps.append(0)
		else:
			if strand == "+":
				UTR_list_s = sorted(UTR_list, key=lambda x: x[1])
				start_UTR = UTR_list_s[0][1]
				end_UTR   = UTR_list_s[0][2]
				for i in range(1, len(UTR_list_s)):
					diff_to_next = int(UTR_list_s[i][1]) - int(UTR_list_s[i-1][2])
					gaps.append(diff_to_next)
			if strand == "-":
				#print(UTR_list)
				UTR_list_s = sorted(UTR_list, key=lambda x: x[1], reverse = True)
				#print(UTR_list_s)
				start_UTR = UTR_list_s[0][1]
				end_UTR   = UTR_list_s[0][2]
				for i in range(1, len(UTR_list_s)):
					diff_to_next =  int(UTR_list_s[i-1][1])  - int(UTR_list_s[i][2])
					gaps.append(diff_to_next)
		return(gaps)
	

	def five_UTR_gaps(UTR_list, strand):
		gaps = []
		if len(UTR_list) == 1:
			gaps.append(0)
		else:
			if strand == "-":
				UTR_list_s = sorted(UTR_list, key=lambda x: x[1])
				start_UTR = UTR_list_s[0][1]
				end_UTR   = UTR_list_s[0][2]
				for i in range(1, len(UTR_list_s)):
					diff_to_next = int(UTR_list_s[i][1]) - int(UTR_list_s[i-1][2])
					gaps.append(diff_to_next)
			if strand == "+":
				#print(UTR_list)
				UTR_list_s = sorted(UTR_list, key=lambda x: x[1], reverse = True)
				#print(UTR_list_s)
				start_UTR = UTR_list_s[0][1]
				end_UTR   = UTR_list_s[0][2]
				for i in range(1, len(UTR_list_s)):
					diff_to_next =  int(UTR_list_s[i-1][1])  - int(UTR_list_s[i][2])
					gaps.append(diff_to_next)
		
		return(gaps)
	
	orig_gff_file = open(orig_gff_file_name)
	
	for line in orig_gff_file:
		if not line.startswith("#"):
			line_o = line
			line = line.rstrip("\n").split("\t")
			feature = line[2]
			ID = line[8].split("ID=")[1].split(";")[0].strip()
			mRNA_loc = mRNA_pos_to_ID_dict.get(ID)
			strand_c = line[6]
			
			if feature == "mRNA":
				gene_ID = ID.split("-")
				if len(gene_ID) != 2:
					print("gene name error. Exit")
					sys.exit(2)
				gene_ID = gene_ID[0]
				# print(line)
				# print(mRNA_loc)
				
				UTR_mRNA_ID = UTR_mRNA_pos_to_ID_dict.get(mRNA_loc)
				if UTR_mRNA_ID != None:
					if len(UTR_mRNA_ID) == 1:
						UTR_mRNA_ID = UTR_mRNA_ID[0]
					elif len(UTR_mRNA_ID) > 1:
						five_test  = set()
						three_test = set()
						for iso in UTR_mRNA_ID:
							five_UTR_rec_a  = five_UTR_dict.get(iso) 
							three_UTR_rec_a = three_UTR_dict.get(iso)
							if five_UTR_rec_a != None:
								five_test.add(tuple(five_UTR_rec_a))
							if three_UTR_rec_a != None:
								three_test.add(tuple(three_UTR_rec_a))
						
						#print(len(three_test))
						if len(five_test) <= 1 and len(three_test) <= 1:
							UTR_mRNA_ID = UTR_mRNA_ID[0]
						else:
							UTR_mRNA_ID = UTR_mRNA_ID[0]
						
					if UTR_mRNA_ID != None:
						five_UTR_rec  = five_UTR_dict.get(UTR_mRNA_ID) 
						three_UTR_rec = three_UTR_dict.get(UTR_mRNA_ID) 
						if five_UTR_rec != None:
							#print(five_UTR_rec)
							five_UTR_rec_fix  = five_UTR_gaps(five_UTR_rec, strand_c)
							for gi in five_UTR_rec_fix:
								gap_file.write(ID + "\t" + str(gi) + "\n")
								every_gap.append(gi)
							#print(five_UTR_rec_fix)

						if three_UTR_rec != None:
							#print(three_UTR_rec)
							three_UTR_rec_fix = three_UTR_gaps(three_UTR_rec, strand_c)
							for gi in three_UTR_rec_fix:
								gap_file.write(ID + "\t" + str(gi) + "\n")
								every_gap.append(gi)
			
						
	import statistics
	print("Median gap size: " + str(statistics.median(every_gap)))


		


 
#### add UTR to orig transcripts

orig_gff_file = open(orig_gff_file_name)
temp_fixgff_added_UTR = open(orig_gff_file_name.replace(".gff", "_temp_fix_1.gff"), "w")

for line in orig_gff_file:
	if not line.startswith("#"):
		line_o = line
		line = line.rstrip("\n").split("\t")
		feature = line[2]
		ID = line[8].split("ID=")[1].split(";")[0].strip()
		mRNA_loc = mRNA_pos_to_ID_dict.get(ID)
		strand_c = line[6]
		temp_fixgff_added_UTR.write(line_o)
		
		
		
		if feature == "mRNA":
			gene_ID = ID.split("-")
			if len(gene_ID) != 2:
				print("gene name error. Exit")
				sys.exit(2)
			gene_ID = gene_ID[0]
			# print(line)
			# print(mRNA_loc)
			
			UTR_mRNA_ID = UTR_mRNA_pos_to_ID_dict.get(mRNA_loc)
			#print(UTR_mRNA_ID)
			if UTR_mRNA_ID != None:
				if len(UTR_mRNA_ID) == 1:
					UTR_mRNA_ID = UTR_mRNA_ID[0]
					
				elif len(UTR_mRNA_ID) > 1:
					five_test  = set()
					three_test = set()
					for iso in UTR_mRNA_ID:
						five_UTR_rec_a  = five_UTR_dict.get(iso) 
						three_UTR_rec_a = three_UTR_dict.get(iso)
						if five_UTR_rec_a != None:
							five_test.add(tuple(five_UTR_rec_a))
						if three_UTR_rec_a != None:
							three_test.add(tuple(three_UTR_rec_a))
					
					#print(len(three_test))
					## if multi transcripts but onlt one has UTR - keep UTR
					if len(five_test) <= 1 and len(three_test) <= 1:
						UTR_mRNA_ID = UTR_mRNA_ID[0]
					## else merge them together.
					else:
						all_five_recs = []
						all_three_recs = []
						for iso in UTR_mRNA_ID:
							five_UTR_rec_a  = five_UTR_dict.get(iso)
							three_UTR_rec_a = three_UTR_dict.get(iso)
							if five_UTR_rec_a != None:
								for xi in five_UTR_rec_a:
									all_five_recs.append(xi)
							if three_UTR_rec_a != None:
								for yi in three_UTR_rec_a:
									all_three_recs.append(yi)
								
							# print(five_UTR_rec_a)
							# print(three_UTR_rec_a)
					
						UTR_mRNA_ID = UTR_mRNA_ID[0] + "_merged"
						#print(UTR_mRNA_ID)
						if all_five_recs  != []:
							merge_iso_rec = merge_ranges(all_five_recs)
							five_UTR_dict[UTR_mRNA_ID] = merge_iso_rec
						else:
							five_UTR_dict[UTR_mRNA_ID] = None
							
						if all_three_recs  != []:
							merge_iso_rec = merge_ranges(all_three_recs)
							three_UTR_dict[UTR_mRNA_ID] = merge_iso_rec
						else:
							three_UTR_dict[UTR_mRNA_ID] = None
						
						#UTR_mRNA_ID = UTR_mRNA_ID[0]
					
				if UTR_mRNA_ID != None:
					five_UTR_rec  = five_UTR_dict.get(UTR_mRNA_ID) 
					three_UTR_rec = three_UTR_dict.get(UTR_mRNA_ID)
					if five_UTR_rec != None:
						#print(five_UTR_rec)
						five_UTR_rec_fix  = five_UTR_clip(five_UTR_rec, strand_c)
						for x in five_UTR_rec_fix:
							temp_fixgff_added_UTR.write(line[0] + "\t" + "AnnotationFinalizer" + "\t" + "5'-UTR" + "\t" +
													str(x[1]) + "\t" + str(x[2]) + "\t" + "." + "\t" + strand_c + "\t" + "." +
													"\tID=" + ID + ":5'-UTR;Parent=" + ID + ";gene_id=" + gene_ID + "\n")	
					
					if three_UTR_rec != None:
						#print(three_UTR_rec)
						three_UTR_rec_fix = three_UTR_clip(three_UTR_rec, strand_c)
						#print(three_UTR_rec_fix)
						for x in three_UTR_rec_fix:
							#print(x)
							temp_fixgff_added_UTR.write(line[0] + "\t" + "AnnotationFinalizer" + "\t" + "3'-UTR" + "\t" +
													str(x[1]) + "\t" + str(x[2]) + "\t" + "." + "\t" + strand_c + "\t" + "." +
													"\tID=" + ID + ":3'-UTR;Parent=" + ID + ";gene_id=" + gene_ID + "\n")				
					
# # Tps_LRv5b_scf2	AnnotationFinalizer	5'-UTR	136904479	136904585	.	+	.	ID=Tps_015122-RB:5'-UTR;Parent=Tps_015122-RB;gene_id=Tps_015122
# # Tps_LRv5b_scf2	AnnotationFinalizer	3'-UTR	136907901	136907912	.	+	.	ID=Tps_015122-RB:3'-UTR;Parent=Tps_015122-RB;gene_id=Tps_015122

temp_fixgff_added_UTR.close()


##########################################################################
### Fix mRNA and gene coords to match UTR


### genes

temp_fixgff_added_UTR = open(orig_gff_file_name.replace(".gff", "_temp_fix_1.gff"))
feat_want = set(["mRNA", "gene", "5'-UTR", "3'-UTR"])
new_gene_coords_dict = {}
seen_gene = set()

for line in temp_fixgff_added_UTR:
	line = line.rstrip("\n").split("\t")
	feature = line[2]
	start_c = int(line[3])
	end_c   = int(line[4])


	if feature in feat_want:
		gene_ID = None
		if feature == "gene":	
			gene_ID = line[8].split("ID=")[1].split(";")[0]
		else:
			gene_ID = line[8].split("gene_id=")[1].split(";")[0]
		
		if gene_ID not in seen_gene:
			seen_gene.add(gene_ID)
			new_gene_coords_dict[gene_ID] = [start_c, end_c]
		else:
			rec = new_gene_coords_dict.get(gene_ID)
			rec = rec + [start_c, end_c]
			new_gene_coords_dict[gene_ID] = rec

temp_fixgff_added_UTR.close()		


### mRNA # do UTR per mRNA
temp_fixgff_added_UTR = open(orig_gff_file_name.replace(".gff", "_temp_fix_1.gff"))
feat_want = set(["mRNA", "5'-UTR", "3'-UTR"])
new_mRNA_coords_dict = {}
seen_mRNA = set()

for line in temp_fixgff_added_UTR:
	line = line.rstrip("\n").split("\t")
	feature = line[2]
	start_c = int(line[3])
	end_c   = int(line[4])


	if feature in feat_want:
		mRNA_ID = None
		if feature == "mRNA":	
			mRNA_ID = line[8].split("ID=")[1].split(";")[0]
		else:
			mRNA_ID = line[8].split("Parent=")[1].split(";")[0]
		
		#print(mRNA_ID)
		
		if mRNA_ID not in seen_mRNA:
			seen_mRNA.add(mRNA_ID)
			new_mRNA_coords_dict[mRNA_ID] = [start_c, end_c]
		else:
			rec = new_mRNA_coords_dict.get(mRNA_ID)
			rec = rec + [start_c, end_c]
			new_mRNA_coords_dict[mRNA_ID] = rec

temp_fixgff_added_UTR.close()		

## output to new temp
temp_fixgff_added_UTR_2 = open(orig_gff_file_name.replace(".gff", "_temp_fix_2.gff"), "w")

#temp_fixgff_added_UTR_2.write("##gff-version 3\n")

temp_fixgff_added_UTR = open(orig_gff_file_name.replace(".gff", "_temp_fix_1.gff"))
for line in temp_fixgff_added_UTR:
	line_o = line
	line = line.rstrip("\n").split("\t")
	feature = line[2]
	
	if feature == "gene":
		gene_ID = line[8].split("ID=")[1].split(";")[0]
		coord_rec = new_gene_coords_dict.get(gene_ID)
		temp_fixgff_added_UTR_2.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" +
									  str(min(coord_rec)) + "\t" + str(max(coord_rec)) + "\t" +
									  line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\n")

	elif feature == "mRNA":
		mRNA_ID = line[8].split("ID=")[1].split(";")[0]
		coord_rec = new_mRNA_coords_dict.get(mRNA_ID)
		temp_fixgff_added_UTR_2.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" +
									  str(min(coord_rec)) + "\t" + str(max(coord_rec)) + "\t" +
									  line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\n")

	else:
		temp_fixgff_added_UTR_2.write(line_o)


temp_fixgff_added_UTR.close()
temp_fixgff_added_UTR_2.close()




####################################################################################################################################
#### merge identical genes


temp_fixgff_added_UTR_2 = open(orig_gff_file_name.replace(".gff", "_temp_fix_2.gff"))

seen_gene_pos = set()
gene_pos_ID_dict = {}
gene_orient_dict = {}

for line in temp_fixgff_added_UTR_2:
	line_o = line
	line = line.rstrip("\n").split("\t")
	feature = line[2]
	if feature == "gene":
		scaf    = line[0]
		start_c = int(line[3])
		end_c   = int(line[4])
		gene_ID = line[8].split("ID=")[1].split(";")[0]
		orient  = line[6] ## if different dir - keep sep.
		pos_id  = (scaf , start_c, end_c, orient)

		if pos_id  not in seen_gene_pos:
			seen_gene_pos.add(pos_id)
			gene_pos_ID_dict[pos_id] = [gene_ID]
		else:
			rec = gene_pos_ID_dict.get(pos_id)
			rec.append(gene_ID)
			gene_pos_ID_dict[pos_id] = rec
		

rename_gene_dict = {}
renamed_genes = set()

for p in gene_pos_ID_dict:
	rec = gene_pos_ID_dict.get(p)
	if len(rec) > 1:
		#print(p)
		#print(rec)
		for i in range(1, len(rec)):
			rename_gene_dict[rec[i]] = rec[0]
			renamed_genes.add(rec[0])
			

rename_mRNA_dict = {}

seen_gene = set()
temp_fixgff_added_UTR_2 = open(orig_gff_file_name.replace(".gff", "_temp_fix_2.gff"))
for line in temp_fixgff_added_UTR_2:
	line_o = line
	line = line.rstrip("\n").split("\t")
	feature = line[2]
	if feature == "mRNA":
		ID      = line[8].split("ID=")[1].split(";")[0]
		gene_ID = line[8].split("gene_id=")[1].split(";")[0]
		new_gene_ID = rename_gene_dict.get(gene_ID)
		if new_gene_ID == None:
			new_gene_ID = gene_ID
		
		if new_gene_ID not in seen_gene:
			seen_gene.add(new_gene_ID)
			rename_mRNA_dict[new_gene_ID] = [ID]
		else:
			rec = rename_mRNA_dict.get(new_gene_ID)
			rec.append(ID)
			rename_mRNA_dict[new_gene_ID] = rec
			
temp_fixgff_added_UTR_2.close()

print("N merged genes: " + str(len(renamed_genes)))

#### function to label transcripts

def N_to_trans_letter(N_trans):
	from string import ascii_lowercase
	import itertools

	def iter_all_strings():
		for size in itertools.count(1):
			for s in itertools.product(ascii_lowercase, repeat=size):
				yield "".join(s)
	last = ""
	for s in itertools.islice(iter_all_strings(), N_trans):
		#print(s)
		last = s
	return(last)


rename_mRNA_dict_2 = {}
for p in rename_mRNA_dict:
	if p in renamed_genes:
		rec = rename_mRNA_dict.get(p)
		# print(p)
		# print(rec)
		for i in range(1,  len(rec) + 1):
			# print(rec[i-1])
			# print(p + "-R" + N_to_trans_letter(i).upper())
			rename_mRNA_dict_2[rec[i-1]] = p + "-R" + N_to_trans_letter(i).upper()
		
		

### merge genes (by renaming)
temp_fixgff_added_UTR_3 = open(orig_gff_file_name.replace(".gff", "_temp_fix_3.gff"), "w")
#temp_fixgff_added_UTR_3.write("##gff-version 3\n")
temp_fixgff_added_UTR_2 = open(orig_gff_file_name.replace(".gff", "_temp_fix_2.gff"))
N_dup = 0

seen_line = set()
for line in temp_fixgff_added_UTR_2:
	line_o = line
	line = line.rstrip("\n").split("\t")
	feature = line[2]
	
	if feature == "gene":
		gene_ID = line[8].split("ID=")[1].split(";")[0]
		new_gene_ID = rename_gene_dict.get(gene_ID)
		if new_gene_ID == None:
			new_gene_ID = gene_ID
		
		line_o = line_o.replace(gene_ID, new_gene_ID)
			
	elif feature == "mRNA":
		gene_ID = line[8].split("gene_id=")[1].split(";")[0]
		new_gene_ID = rename_gene_dict.get(gene_ID)
		if new_gene_ID == None:
			new_gene_ID = gene_ID
		
		mRNA_ID = line[8].split("ID=")[1].split(";")[0]
		new_mRNA_ID = rename_mRNA_dict_2.get(mRNA_ID)
		if new_mRNA_ID == None:
			new_mRNA_ID = mRNA_ID
		
		line_o = line_o.replace(mRNA_ID, new_mRNA_ID)
		line_o = line_o.replace(gene_ID, new_gene_ID)
	
	else:
		gene_ID = line[8].split("gene_id=")[1].split(";")[0]
		new_gene_ID = rename_gene_dict.get(gene_ID)
		if new_gene_ID == None:
			new_gene_ID = gene_ID
		
		mRNA_ID = line[8].split("Parent=")[1].split(";")[0]
		new_mRNA_ID = rename_mRNA_dict_2.get(mRNA_ID)
		if new_mRNA_ID == None:
			new_mRNA_ID = mRNA_ID
			
		
		line_o = line_o.replace(mRNA_ID, new_mRNA_ID)
		line_o = line_o.replace(gene_ID, new_gene_ID)
	
	
	### output all line, but drop dup genes 
	
	if feature == "gene":
		line_start_a = line_o.split(";")[0].strip().split("\t")
		line_start = line_start_a[0] + "XX" + line_start_a[2] +  "XX" + line_start_a[3] +  "XX" + line_start_a[4] +  "XX" + line_start_a[5] +  "XX" + line_start_a[6] +  "XX" + line_start_a[7]  + "XX" + line_start_a[8]  ## drop source of annot as crit
		if line_start not in seen_line:
			seen_line.add(line_start)
			temp_fixgff_added_UTR_3.write(line_o)
		else:
			#print(line_o)
			N_dup = N_dup + 1
	else:
		temp_fixgff_added_UTR_3.write(line_o)

temp_fixgff_added_UTR_3.close()

# print(N_dup )

### ordering of records

###### need to rank within gene


# get_gene_start_coord

gene_start_dict = {}
temp_fixgff_added_UTR_3 = open(orig_gff_file_name.replace(".gff", "_temp_fix_3.gff"))
for line in temp_fixgff_added_UTR_3:
	line_o = line
	line = line.rstrip("\n").split("\t")
	scaf    = line[0]
	start_c = int(line[3])
	feature = line[2]
	ID = line[8].split("ID=")[1].split(";")[0]	
	if feature == "gene":
		gene_ID = line[8].split("ID=")[1].split(";")[0]
		gene_start_dict[gene_ID] = start_c 
temp_fixgff_added_UTR_3.close()

out_fixgff_file         = open(orig_gff_file_name.replace(".gff", "UGF" + str(max_utr_gap) + "M" + str(merge_n) + ".gff"), "w")
out_fixgff_file.write("##gff-version 3\n")

all_record_list = []
all_rec_dict    = {}
temp_fixgff_added_UTR_3 = open(orig_gff_file_name.replace(".gff", "_temp_fix_3.gff"))
for line in temp_fixgff_added_UTR_3:
	line_o = line
	line = line.rstrip("\n").split("\t")
	scaf    = line[0]
	start_c = int(line[3])
	feature = line[2]
	feat_class = 9999
	gene_ID = None
	ID = line[8].split("ID=")[1].split(";")[0]	
	if feature == "gene":
		gene_ID = line[8].split("ID=")[1].split(";")[0]
		mRNA_ID = "0"
	elif feature == "mRNA":
		gene_ID = line[8].split("gene_id=")[1].split(";")[0]
		mRNA_ID = "1" + line[8].split("ID=")[1].split(";")[0]
	else:
		gene_ID = line[8].split("gene_id=")[1].split(";")[0]
		mRNA_ID = "1" + line[8].split("Parent=")[1].split(";")[0]
		
	gene_start = gene_start_dict.get(gene_ID)	
	if feature == "gene":
		feat_class = 1
	if feature == "mRNA":
		feat_class = 2
	if feature == "start_codon":
		feat_class = 3
	if feature == "stop_codon":
		feat_class = 3
	if feature == "CDS":
		feat_class = 4	
	if feature == "exon":
		feat_class = 5
		
	rec_ID  = (scaf, gene_start, mRNA_ID, start_c, feat_class, ID)
	all_record_list.append(rec_ID)
	all_rec_dict[rec_ID] = line_o




### sort

# aa = [("scaf_1",  int(99), int(1)), ("scaf_1",  int(10000000), int(2)),  ("scaf_1",  int(10000000), int(1)), ("scaf_10",  int(1), int(2)), ("scaf_2",  int(1), int(2))]
# print(sorted(aa, key=lambda e: (e[0], e[1], e[2])))


all_record_list_s = sorted(all_record_list, key=lambda e: (e[0], e[1], e[2], e[3], e[4], e[5]))

for r in all_record_list_s:
	out_l = all_rec_dict.get(r)
	#print(out_l)
	out_fixgff_file.write(out_l)
	
	
temp_fixgff_added_UTR_3.close()

# os.remove(orig_gff_file_name.replace(".gff", "_temp_fix_1.gff"))
# os.remove(orig_gff_file_name.replace(".gff", "_temp_fix_2.gff"))
# os.remove(orig_gff_file_name.replace(".gff", "_temp_fix_3.gff"))


print("\n\nFinished Arthur Dent\n\n")

