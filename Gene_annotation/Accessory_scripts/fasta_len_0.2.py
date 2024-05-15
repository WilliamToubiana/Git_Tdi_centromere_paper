## fasta len 

import sys
import os

args = sys.argv
arg_len = len(args)

if arg_len <2:
	print("\n**** Written by DJP, 24/11/15 in Python 3.4 ****\n")
	### update 18/05/16 - also spits out N of seqs and it unwraps the fastafile
	### update 08/06/16 - altered output dormat to make it more useful
	print("This program takes a fasta file as input.") 
	print("It outputs, Number of seqs, min length, max length, mean length and N50 of fasta sequences to screen and output file.")
	print("\n**** USAGE **** \n")
	print("fasta_len_v0.1.py [name of fasta file] \n")

else:
	input_fasta = args[1]


	##### FIRST unwrap fasta - precautionary will be necessary for some files 
	### note making a temp unwrapped fasta file  - removed at end
	output_fasta_name = input_fasta + ".TEMP_extract_fasta_file" 

	output_file = open(output_fasta_name, "w")
	print("Unwrapping fasta file\n")
	count = 0
	in_file = open(input_fasta)
	for line in in_file:
		count = count + 1
		line = line.rstrip("\n")
		if line.startswith(">") and count == 1:
			output_file.write(line + "\n")
		elif line.startswith(">") and count > 1:
			output_file.write("\n" + line + "\n")
		else: 
			output_file.write(line)	
	
	output_file.close()
	
	
	name = []
	seq = []
	all_seq = ""
	
	seq_dict = {}
	
	seq_count = 0
	file_in = open(output_fasta_name)
	for line in file_in:
		line = line.rstrip("\n")
		if line.startswith(">"):
			seq_count = seq_count + 1
			line = line.replace(">", "")
			name.append(line)
		else:
			seq.append(line)
			all_seq = all_seq + line

	#print(all_seq)
	
	out_bases_l = []
	from collections import Counter
	basecounts = Counter(all_seq)
	for el in basecounts:
		out_bases_l.append(el + ":" + str(basecounts.get(el)))

	out_bases_l_s = sorted(out_bases_l, reverse = True)
	out_bases = ""
	for i in out_bases_l_s:
		a1 = i
		out_bases = a1 + ',' + out_bases
	
	out_bases = out_bases.rstrip(",")


	out_bases_l_cap = []
	basecounts = Counter(all_seq.upper())
	for el in basecounts:
		out_bases_l_cap.append(el + ":" + str(basecounts.get(el)))

	out_bases_l_s_cap = sorted(out_bases_l_cap, reverse = True)
	out_bases_cap = ""
	for i in out_bases_l_s_cap:
		a1 = i
		out_bases_cap = a1 + ',' + out_bases_cap
	
	out_bases_cap = out_bases_cap.rstrip(",")	
	
	
	
	for element in range(0,len(name)):
		name1 = name[element]
		seq1 = seq[element]
		seq_dict[name1] = seq1

	#print(seq_dict)	
	#print(seq_dict.get("contig_99"))

	total_Ns = 0 
	len_list = []
	for el in seq_dict:
		seq_len = len(seq_dict.get(el))
		len_list.append(seq_len)
		Ns = seq_dict.get(el).upper().count("N")
		total_Ns = total_Ns + Ns
		#print(Ns)

		
	a1 = len(len_list)
	sum_of_len = sum(len_list)

	mean_len = sum_of_len / a1

	min_len = min(len_list)	
	max_len = max(len_list)

	### calc N50
	sorted_len_list  = sorted(len_list, reverse=True)
	sum_of_len_2 = sum_of_len / 2

	i = 0
	num_tot = 0
	for num in sorted_len_list:
		if sum_of_len_2 > num_tot:
			num_tot = num_tot + num
			i = i + 1
			#print(num)
			#print(i)

	N50 = sorted_len_list[i-1]
	perc_Ns = total_Ns / sum_of_len * 100

	outname = input_fasta + ".stat.txt"
	outname2 = input_fasta + ".base_comp.txt"
	outputfile = open(outname, "w")
	
	outputfile.write("Filename" + "\ttotal_len\ttotal_Ns\t%Ns\t" + "N_seqs" + "\t" + "shortest_seq_len" + "\t" + "longest_seq_len" + "\t" + "mean_len" + "\t" + "N50" + "\n") 
	outputfile.write(input_fasta + "\t" + str(sum_of_len) + "\t" + str(total_Ns) + "\t" + str(perc_Ns) + "\t" + str(seq_count) + "\t" + str(min_len) + "\t" + str(max_len) + "\t" + str(mean_len) + "\t" + str(N50)+ "\n")
	
	outputfile2 = open(outname2, "w")
	outputfile2.write(input_fasta + "," + out_bases + "\n")
	
	#print(len_list)
	
	print(input_fasta + ":")
	print("The number of sequences are: " + str(seq_count))
	print("The total number of base pairs are: " + str(sum_of_len))
	print("The percentage of bases that are Ns: " + str(perc_Ns) + "%")
	print("The total number of Ns are: " + str(total_Ns))
	print("The shortest sequence is: " + str(min_len) + " BP long.")	
	print("The longest sequence is: " + str(max_len) + " BP long.")	
	print("The mean is: " + str(mean_len))	
	print("The N50 is: " + str(N50))
	print("The total number of each base is: " + str(out_bases))
	print("The total number of each base is: " + str(out_bases_cap) + "\n")
	#print(sorted_len_list)

	## tidyup
	file_in.close()
	os.remove(output_fasta_name)




