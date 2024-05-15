import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:w:h')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

in_file_name = None
outprefix    = "testout"
window_size   = 10000

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** bedgraph_cov_to_windows.py | Written by DJP, 14/05/21 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Takes a bedgraph file from bedtools (i.e. made with genomeCoverageBed -ibam mybam.bam -bga) and puts coverage into specified windows\nWindows slide by 1/2 the length of the specified window")	
		print("\n**** Usage****\n")
		print("python3 bedgraph_cov_to_windows.py -i [in file] -o [outprefix] -w [window size]\n\n")
		sys.exit(2)
	elif opt in ('-i'):
		in_file_name  = arg
	elif opt in ('-w'):		
		window_size   = int(arg)
	elif opt in ('-o'):
		outprefix = arg
	else:
		print("i dont know")
		sys.exit(2)

if in_file_name == None:
	print("\n\nERROR. No input file\n\n")

if window_size % 2 == 0:
	print("\nWindow size = " + str(window_size))
else:
	window_size = window_size + 1
	print("\nWindow size specified was not divisable by 2. Window size set to " + str(window_size))


### read in file

### get max scaf

### read in file
max_scaf_len_dict = {}
seen_scaf = set()
in_file = open(in_file_name)
for line in in_file:
	line       = line.rstrip("\n").split("\t")
	scaf_name  = line[0]
	end_c      = int(line[2])
	if scaf_name not in seen_scaf:
		seen_scaf.add(scaf_name)
		max_scaf_len_dict[scaf_name] = end_c 
	else:
		rec = max_scaf_len_dict.get(scaf_name)
		if rec < end_c:
			max_scaf_len_dict[scaf_name] = end_c 

in_file.close()

print("Step 1, done")

### get vals for all positions

window_dict = {}

slide_by    = int(window_size / 2)
curr_total = 0
N_line = 0
line_N_tot = 0

in_file = open(in_file_name)
seen_scaf = set()
for line in in_file:
	line = line.rstrip("\n").split("\t")
	scaf_name  = line[0]
	start_c = int(line[1])
	end_c   = int(line[2])
	cov     = int(line[3])
	max_scaf = max_scaf_len_dict.get(scaf_name)
	if scaf_name not in seen_scaf:
		seen_scaf.add(scaf_name)
		curr_total = 0
		N_line = 0
	#print(line)
	for i in range(start_c, end_c):
		#print(scaf_name + "\t" + str(i) + "\t" + str(i + 1) + "\t" + str(cov))
		N_line = N_line + 1
		line_N_tot = line_N_tot + 1
		if N_line < slide_by :
			curr_total = curr_total + cov
		else:
			
			curr_total = curr_total + cov
			#print(scaf_name + "\t" + str(line_N_tot - N_line) + "\t" + str(line_N_tot) + "\t" + str(curr_total))
			window_dict[scaf_name + ":" + str(line_N_tot - N_line) + "-" + str(line_N_tot)] = curr_total
			curr_total = 0
			N_line = 0
		
		if line_N_tot == max_scaf:
			#print(scaf_name + "\t" + str(line_N_tot - N_line) + "\t" + str(line_N_tot) + "\t" + str(curr_total))
			window_dict[scaf_name + ":" + str(line_N_tot - N_line) + "-" + str(line_N_tot)] = curr_total
			line_N_tot = 0

in_file.close()

print("Step 2, done")

### join windows to make sliding

outfile = open(outprefix + "_cov_window_" + str(window_size) + ".txt", "w")
outfile.write("scaf_name\tstart\tend\ttot_feat\tcov\n")

seen_scaf_sl = sorted(list(seen_scaf))	
for scaf_name in seen_scaf_sl:

	for i in range(0, max_scaf_len_dict.get(scaf_name), slide_by):
		max_scaf = max_scaf_len_dict.get(scaf_name)

		if i + window_size < max_scaf:
			window_1 = window_dict.get(scaf_name + ":" + str(i) + "-" + str(i + slide_by))
			window_2 = window_dict.get(scaf_name + ":" + str(i + slide_by) + "-" + str(i + window_size))
			window_12 = window_1 + window_2
			cov = decimal.Decimal(decimal.Decimal(window_12) / decimal.Decimal(((i + window_size) - i)))
			#print(scaf_name + "\t" + str(i) + "\t" +  str(i + window_size) + "\t" + str(window_12) + "\t" + str(cov))
			outfile.write(scaf_name + "\t" + str(i) + "\t" +  str(i + window_size) + "\t" + str(window_12) + "\t" + str(cov) + "\n")

	
		if i + window_size >= max_scaf:
			window_1 = window_dict.get(scaf_name + ":" + str(i) + "-" + str(i + slide_by))
			window_2 = window_dict.get(scaf_name + ":" + str(i + slide_by) + "-" + str(i + window_size))
			
			if window_1 == None and window_2 == None:
				window_12 = window_dict.get(scaf_name + ":" + str(i) + "-" + str(max_scaf))
				cov = decimal.Decimal(decimal.Decimal(window_12) / decimal.Decimal((max_scaf - i)))
				#print(scaf_name + "\t" + str(i) + "\t" +  str(max_scaf) + "\t" + str(window_12) + "\t" + str(cov))
				outfile.write(scaf_name + "\t" + str(i) + "\t" +  str(max_scaf) + "\t" + str(window_12) + "\t" + str(cov) + "\n")
			elif window_1 != None and window_2 == None:
				window_2 = window_dict.get(scaf_name + ":" + str(i + slide_by) + "-" + str(max_scaf))
				window_12 = window_1 + window_2
				cov = decimal.Decimal(decimal.Decimal(window_12) / decimal.Decimal((max_scaf - i)))
				#print(scaf_name + "\t" + str(i) + "\t" +  str(max_scaf) + "\t" + str(window_12) + "\t" + str(cov))
				outfile.write(scaf_name + "\t" + str(i) + "\t" +  str(max_scaf) + "\t" + str(window_12) + "\t" + str(cov) + "\n")
			else:
				window_12 = window_1 + window_2
				cov = decimal.Decimal(decimal.Decimal(window_12) / decimal.Decimal(((max_scaf - i))))
				#print(scaf_name + "\t" + str(i) + "\t" +  str(max_scaf) + "\t" + str(window_12) + "\t" + str(cov))
				outfile.write(scaf_name + "\t" + str(i) + "\t" +  str(max_scaf) + "\t" + str(window_12) + "\t" + str(cov) + "\n")
			break


print("\n\nAll done, Gully\n\n")
