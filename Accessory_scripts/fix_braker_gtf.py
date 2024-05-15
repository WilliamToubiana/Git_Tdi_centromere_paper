### Maker_gff_to_HTseq_gff.py

import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:p:d:hMG')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_file_name = None
gene_prefix = "Test"
N_geneN_digits = 5
ADD_missing_exons = True
ADD_gene_id = False

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** Maker_gff_to_HTseq_gff.py | Written by DJP, 12/04/22 in Python 3.5 in Bangor, UK. ****\n")
		print("Takes combined gtfs from the braker pipeline and converts them so they can be used properly")
		
		print("\n**** USAGE **** \n")
		print("python3 fix_braker_gtf.py -i [braker gtf file] -p [gene prefix to use] -d [number of digits to use for gene names] \n")
		print("specify -M if you don't want to fix missing exon records")
		print("specify -G if you want to add a gene_id to all non-gene features")
		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-p'):
		gene_prefix = arg
	elif opt in ('-d'):
		N_geneN_digits = int(arg)
	elif opt in ('-M'):
		ADD_missing_exons = False
	elif opt in ('-G'):
		ADD_gene_id = True
	else:
		print("i dont know")
		sys.exit(2)


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


### read gff | get link between transcript name and gene name from CDS

trans_to_gene_dict = {}

gene_names_from_CDS = set()
trans_names_from_CDS = set()
trans_names_all = set()

in_gff = open(in_file_name)
for line in in_gff:
	line = line.rstrip("\n")

	if len(line.split("\t")) >= 5:
		feature = line.split("\t")[2]
		if feature == "CDS":
			descrip_line = line.split("\t")[8]
			#print(descrip_line)
			transcript_id = descrip_line.split("transcript_id")[1].split(";")[0].strip().strip('"')
			#print(transcript_id)
			gene_id = descrip_line.split("gene_id")[1].split(";")[0].strip().strip('"')
			#print(gene_id)
			
			gene_names_from_CDS.add(gene_id)
			trans_names_from_CDS.add(transcript_id)
			trans_to_gene_dict[transcript_id] = gene_id
		if feature == "transcript":
			descrip_line = line.split("\t")[8].strip()
			trans_names_all.add(descrip_line)
in_gff.close()		

print("\nN genes: " + str(len(gene_names_from_CDS)))
print("N transcripts: " + str(len(trans_names_all)))

	

### check we got them all

if len(trans_names_all) != len(trans_names_from_CDS | trans_names_all):
	print("ERROR: some transcripts do not have CDS")
	sys.exit(2)


### set gene numbering length

if N_geneN_digits == None:
	N_geneN_digits == len(str(len(gene_names_from_CDS)))
else:
	if not N_geneN_digits >= len(str(len(gene_names_from_CDS))):
		print("error N gene digits set too low for the number of genes. For this data it should be set to " + str(len(str(len(gene_names_from_CDS)))) + " or higher.")
		sys.exit(2)



###################################################
## correct

seen_gene = set()
seen_gene_trans = set()
seen_exon_trans = set()
trans_orientation_dict = {}
gene_to_newgene_name_dict = {}
newgene_to_oldtrans_name_dict = {}
trans_to_newtrans_dict = {}
new_trans_to_exon_dict = {}

temp_out_file = open(in_file_name + ".temp", "w")

line_N = 0
gene_N = 0
in_gff = open(in_file_name)
for line in in_gff:
	line_N = line_N + 1
	line = line.rstrip("\n").split("\t")
	if len(line) >= 5:
		feature = line[2]

		if feature == "gene":
			gene_N = gene_N + 1
			descrip_line = line[8].strip().split(";")[0].strip()
			if descrip_line not in seen_gene:
				seen_gene.add(descrip_line)
				new_gene_name = gene_prefix + "_" + str(gene_N).zfill(N_geneN_digits)
				#print(new_gene_name )
				gene_to_newgene_name_dict[descrip_line] = new_gene_name 
			else:
				print("Error!")
				sys.exit(2)
			new_description = "ID=" + new_gene_name + ";Alias=" + descrip_line + ";Name=" + new_gene_name
			# print(descrip_line )
			# print(new_description)
			temp_out_file.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + new_description + "\n")
			
		elif feature == "transcript":
			descrip_line = line[8].strip()
			gene_id = trans_to_gene_dict.get(descrip_line)
			new_gene_name = gene_to_newgene_name_dict.get(gene_id)
			
			if new_gene_name not in seen_gene_trans:
				seen_gene_trans.add(new_gene_name)
				new_trans_name = new_gene_name + "-RA"
				newgene_to_oldtrans_name_dict[new_gene_name] = [descrip_line]
				trans_to_newtrans_dict[descrip_line] = new_trans_name
			else:
				rec = newgene_to_oldtrans_name_dict.get(new_gene_name)
				new_trans_name = new_gene_name + "-R" + N_to_trans_letter(len(rec) + 1).upper()
				rec.append(descrip_line)
				newgene_to_oldtrans_name_dict[new_gene_name] = rec
				trans_to_newtrans_dict[descrip_line] = new_trans_name
				
			new_description = "ID=" + new_trans_name + ";Alias=" + descrip_line + ";Parent=" + new_gene_name + ";Name=" + new_trans_name
			trans_orientation_dict[new_trans_name] = line[6]
			temp_out_file.write(line[0] + "\t" + line[1] + "\t" + "mRNA" + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + new_description + "\n")
		
		## number exons uniquely here - after will renumber with order depending on the strand
		elif feature == "exon":
			descrip_line = line[8]
			transcript_id = descrip_line.split("transcript_id")[1].split(";")[0].strip().strip('"')
			new_trans_id = trans_to_newtrans_dict.get(transcript_id)
			gene_id = descrip_line.split("gene_id")[1].split(";")[0].strip().strip('"')
			new_gene_name = gene_to_newgene_name_dict.get(gene_id)
			new_description = "ID=" + new_trans_id + ":" + feature + str(line_N) +  ";" + "Parent=" + new_trans_id 
			
			if new_trans_id not in seen_exon_trans:
				seen_exon_trans.add(new_trans_id)
				new_trans_to_exon_dict[new_trans_id] = [new_trans_id + ":" + feature + str(line_N)]
			else:
				rec = new_trans_to_exon_dict.get(new_trans_id)
				rec.append(new_trans_id + ":" + feature + str(line_N))
				new_trans_to_exon_dict[new_trans_id] = rec
			
			temp_out_file.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + new_description + "\n")	
			
		else:
			descrip_line = line[8]
			transcript_id = descrip_line.split("transcript_id")[1].split(";")[0].strip().strip('"')
			new_trans_id = trans_to_newtrans_dict.get(transcript_id)
			gene_id = descrip_line.split("gene_id")[1].split(";")[0].strip().strip('"')
			new_gene_name = gene_to_newgene_name_dict.get(gene_id)
			new_description = "ID=" + new_trans_id + ":" + feature +  ";" + "Parent=" + new_trans_id 			
			# print(descrip_line )
			# print(new_description)		
		#print(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8])
			temp_out_file.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + new_description + "\n")

		
in_gff.close()		



##########################################################################################3
### renumber exons

old_exon_new_exon_name_dict = {}

for tra in new_trans_to_exon_dict:
	exons = new_trans_to_exon_dict.get(tra)
	orient = trans_orientation_dict.get(tra)
	if orient == "+":
		e_N = 0
		for e in exons:
			e_N = e_N + 1
			old_exon_new_exon_name_dict[e] = tra + ":exon" + str(e_N)
	elif orient == "-":
		eN_M = len(exons) + 1
		for e in exons:
			eN_M = eN_M - 1
			old_exon_new_exon_name_dict[e] = tra + ":exon" + str(eN_M)
	else:
		print("error")
		sys.exit(2)


output_file = open(in_file_name.replace(".gtf", "_FBG.gff"), "w")
output_file.write("##gff-version 3\n")
N_skipped_intron_line = 0

temp_out_file = open(in_file_name + ".temp")
for line in temp_out_file:
	line = line.strip().split("\t")
	feature = line[2]
	descrip_line = line[8]
	if feature == "exon":
		exon_ID   = descrip_line.split("ID=")[1].split(";")[0]
		parent_ID = descrip_line.split("Parent=")[1].split(";")[0]
		if ADD_missing_exons == False:
			new_exon_ID = old_exon_new_exon_name_dict.get(exon_ID)
		else:
			new_exon_ID = exon_ID
		new_description = "ID=" + new_exon_ID +  ";" + "Parent=" + parent_ID
		output_file.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + new_description + "\n")
	elif feature == "intron":
		N_skipped_intron_line = N_skipped_intron_line + 1
	else:
		output_file.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\n")


print("\nNot outputting " + str(N_skipped_intron_line) + " intron lines\n")	

temp_out_file.close()	
output_file.close()
os.remove(in_file_name + ".temp")



#################################################################################################################################
#### When adding UTR and merging - the exons get missed out for some transcripts.
#### if ADD_missing_exons = TRUE I will add these missing exons by duplicating the CDS entries.

if ADD_missing_exons == True:
	
	
	### look through gff file. collect CDS and exons into dict
	mRNA_IDs = set()
	mRNA_exon_dict = {}
	mRNA_CDS_dict = {}

	mRNA_exon_seen = set()
	mRNA_CDS_seen = set()
	print("\nLooking for missing exons. Skip this with -M")
	in_file = open(in_file_name.replace(".gtf", "_FBG.gff"))
	for line in in_file:
		if not line.startswith("#"):
			line = line.strip().split("\t")
			feature = line[2]
			descrip_line = line[8]
			if feature == "mRNA":
				mRNA_ID   = descrip_line.split("ID=")[1].split(";")[0]
				mRNA_IDs.add(mRNA_ID)
			elif feature == "exon":
				parent_ID = descrip_line.split("Parent=")[1].split(";")[0]
				exon_ID   = descrip_line.split("ID=")[1].split(";")[0]
				if parent_ID not in mRNA_exon_seen:
					mRNA_exon_seen.add(parent_ID)
					mRNA_exon_dict[parent_ID] = [exon_ID]
				else:
					rec = mRNA_exon_dict.get(parent_ID)
					rec.append(exon_ID)
					mRNA_exon_dict[parent_ID] = rec
			elif feature == "CDS":
				parent_ID = descrip_line.split("Parent=")[1].split(";")[0]
				CDS_ID   = descrip_line.split("ID=")[1].split(";")[0]
				if parent_ID not in mRNA_CDS_seen:
					mRNA_CDS_seen.add(parent_ID)
					mRNA_CDS_dict[parent_ID] = [CDS_ID]
				else:
					rec = mRNA_CDS_dict.get(parent_ID)
					rec.append(CDS_ID)
					mRNA_CDS_dict[parent_ID] = rec
		
	
	in_file.close()


	mRNA_no_exon = mRNA_IDs.difference(mRNA_exon_seen)
	print("Number of mRNAs with no exons " + str(len(mRNA_no_exon)))
	missing_ex_with_noCDS = mRNA_no_exon.difference(mRNA_CDS_seen) 
	print("Number of mRNAs with no exons or CDS " + str(len(missing_ex_with_noCDS)))
	print("Adding exon records based on CDS for " + str(len(mRNA_no_exon) - len(missing_ex_with_noCDS)) + " mRNAs")
	
	set_list = [mRNA_no_exon, mRNA_CDS_seen]
	add_exon_IDs = set.intersection(*set_list)	
	
	### add exon lines from CDS line when missing - output to tempfile
	
	in_file = open(in_file_name.replace(".gtf", "_FBG.gff"))
	temp_out_file = open(in_file_name + ".temp", "w")
	seen_new_e = set()
	line_N = 0
	for line in in_file:
		line_N = line_N + 1
		if not line.startswith("#"):
			line_o = line.strip()
			line = line.strip().split("\t")
			#print(line)
			feature = line[2]
			descrip_line = line[8]
			
			if feature == "CDS":
				parent_ID = descrip_line.split("Parent=")[1].split(";")[0]
				CDS_ID   = descrip_line.split("ID=")[1].split(";")[0]
				if parent_ID in add_exon_IDs:
					new_description = "ID=" + parent_ID + ":" + "exon" + str(line_N) +  ";" + "Parent=" + parent_ID
					temp_out_file.write(line_o + "\n" + line[0] + "\t" + line[1] + "\t" + "exon" + "\t" + line[3] + "\t" + line[4] + "\t" + "." + "\t" + line[6] + "\t" + "." + "\t" + new_description + "\n" )
					
					if parent_ID not in seen_new_e:
						seen_new_e.add(parent_ID)
						new_trans_to_exon_dict[parent_ID] = [parent_ID + ":" + "exon" + str(line_N)]
					else:
						rec = new_trans_to_exon_dict.get(parent_ID)
						rec.append(parent_ID + ":" + "exon" + str(line_N))
						new_trans_to_exon_dict[parent_ID] = rec
				
				
				else:
					temp_out_file.write(line_o + "\n")
			else:
				temp_out_file.write(line_o + "\n")
	
	in_file.close()
	temp_out_file.close()


	##########################################################################################3
	### renumber exons
	
	old_exon_new_exon_name_dict = {}
	
	for tra in new_trans_to_exon_dict:
		exons = new_trans_to_exon_dict.get(tra)
		orient = trans_orientation_dict.get(tra)
		if orient == "+":
			e_N = 0
			for e in exons:
				e_N = e_N + 1
				old_exon_new_exon_name_dict[e] = tra + ":exon" + str(e_N)
		elif orient == "-":
			eN_M = len(exons) + 1
			for e in exons:
				eN_M = eN_M - 1
				old_exon_new_exon_name_dict[e] = tra + ":exon" + str(eN_M)
		else:
			print("error")
			sys.exit(2)

	output_file = open(in_file_name.replace(".gtf", "_FBG.gff"), "w")
	output_file.write("##gff-version 3\n")
	N_skipped_intron_line = 0
	
	temp_out_file = open(in_file_name + ".temp")
	for line in temp_out_file:
		line = line.strip().split("\t")
		feature = line[2]
		descrip_line = line[8]
		if feature == "exon":
			exon_ID   = descrip_line.split("ID=")[1].split(";")[0]
			parent_ID = descrip_line.split("Parent=")[1].split(";")[0]
			new_exon_ID = old_exon_new_exon_name_dict.get(exon_ID)
			new_description = "ID=" + new_exon_ID +  ";" + "Parent=" + parent_ID
			output_file.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + new_description + "\n")
		elif feature == "intron":
			N_skipped_intron_line = N_skipped_intron_line + 1
		else:
			output_file.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\n")
		
	
	temp_out_file.close()	
	output_file.close()
	os.remove(in_file_name + ".temp")

####### add gene_id flag to all non-gene ids
if ADD_gene_id == True:
	print("\nAdding gene_id to all non-gene records as -G is specified")	
	
	output_file = open(in_file_name.replace(".gtf", "_FBGgi.gff"), "w")
	output_file.write("##gff-version 3\n")
	

	in_file = open(in_file_name.replace(".gtf", "_FBG.gff"))
	for line in in_file:
		if not line.startswith("#"):
			line = line.strip().split("\t")
			feature = line[2]
			descrip_line = line[8]
			if feature == "gene":
				output_file.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\n")
			else:
				gene_ID = descrip_line.split("Parent=")[1].split(";")[0].split("-")[0]
				output_file.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + ";gene_id=" + gene_ID + "\n")
		
	output_file.close()
	in_file.close()
	os.remove(in_file_name.replace(".gtf", "_FBG.gff"))
	
print("\n\nFinished, Kelvin\n\n")






#print(new_trans_to_exon_dict)

#print(trans_orientation_dict)




# awk '{print $3}'  for_jeli/Tps_braker_RNAseq_run_4_prot_run_4_combined_1.gtf | sort | uniq
# CDS
# exon
# gene
# intron
# start_codon
# stop_codon
# transcript


# Tps_LRv5b_scf1  AUGUSTUS        gene    5445    5696    .       +       .       g_64091
# Tps_LRv5b_scf1  AUGUSTUS        transcript      5445    5696    .       +       .       anno2.Tps_LRv5b_scf1+_file_1_file_1_g26508.t1
# Tps_LRv5b_scf1  AUGUSTUS        start_codon     5445    5447    .       +       0       transcript_id "anno2.Tps_LRv5b_scf1+_file_1_file_1_g26508.t1"; gene_id "g_64091";
# Tps_LRv5b_scf1  AUGUSTUS        CDS     5445    5696    0.56    +       0       transcript_id "anno2.Tps_LRv5b_scf1+_file_1_file_1_g26508.t1"; gene_id "g_64091";
# Tps_LRv5b_scf1  AUGUSTUS        exon    5445    5696    .       +       .       transcript_id "anno2.Tps_LRv5b_scf1+_file_1_file_1_g26508.t1"; gene_id "g_64091";
# Tps_LRv5b_scf1  AUGUSTUS        stop_codon      5694    5696    .       +       0       transcript_id "anno2.Tps_LRv5b_scf1+_file_1_file_1_g26508.t1"; gene_id "g_64091";


# ##gff-version 3
# Tdi_LRv5a_scf1  maker   gene    888196  889979  .       -       .       ID=TDI_S1_0001;Alias=maker-Tdi_LRv5a_scf1-snap-gene-3.1;Name=TDI_S1_0001
# Tdi_LRv5a_scf1  maker   mRNA    888196  889979  .       -       .       ID=TDI_S1_0001-RA;Parent=TDI_S1_0001;Alias=maker-Tdi_LRv5a_scf1-snap-gene-3.1-mRNA-1;Name=TDI_S1_0001-RA;_AED=0.21;_QI=0|0.33|0.25|0.75|0|0.5|4|376|161;_eAED=0.21
# Tdi_LRv5a_scf1  maker   exon    888196  888649  .       -       .       ID=TDI_S1_0001-RA:4;Parent=TDI_S1_0001-RA
# Tdi_LRv5a_scf1  maker   exon    888761  888915  .       -       .       ID=TDI_S1_0001-RA:3;Parent=TDI_S1_0001-RA
# Tdi_LRv5a_scf1  maker   exon    888967  889111  .       -       .       ID=TDI_S1_0001-RA:2;Parent=TDI_S1_0001-RA
# Tdi_LRv5a_scf1  maker   exon    889872  889979  .       -       .       ID=TDI_S1_0001-RA:1;Parent=TDI_S1_0001-RA
# Tdi_LRv5a_scf1  maker   CDS     888572  888649  .       -       0       ID=TDI_S1_0001-RA:cds;Parent=TDI_S1_0001-RA
# Tdi_LRv5a_scf1  maker   CDS     888761  888915  .       -       2       ID=TDI_S1_0001-RA:cds;Parent=TDI_S1_0001-RA
# Tdi_LRv5a_scf1  maker   CDS     888967  889111  .       -       0       ID=TDI_S1_0001-RA:cds;Parent=TDI_S1_0001-RA
# Tdi_LRv5a_scf1  maker   CDS     889872  889979  .       -       0       ID=TDI_S1_0001-RA:cds;Parent=TDI_S1_0001-RA
# Tdi_LRv5a_scf1  maker   three_prime_UTR 888196  888571  .       -       .       ID=TDI_S1_0001-RA:three_prime_utr;Parent=TDI_S1_0001-RA
# Tdi_LRv5a_scf1  maker   gene    1086430 1098116 .       +       .       ID=TDI_S1_0002;Alias=snap_masked-Tdi_LRv5a_scf1-processed-gene-4.10;Name=TDI_S1_0002
# Tdi_LRv5a_scf1  maker   mRNA    1086430 1098116 .       +       .       ID=TDI_S1_0002-RA;Parent=TDI_S1_0002;Alias=snap_masked-Tdi_LRv5a_scf1-processed-gene-4.10-mRNA-1;Name=TDI_S1_0002-RA;_AED=0.85;_QI=0|0|0.16|0.16|1|1|6|0|181;_eAED=0.85;_merge_warning=1
# Tdi_LRv5a_scf1  maker   exon    1086430 1086437 .       +       .       ID=TDI_S1_0002-RA:1;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   exon    1087724 1087864 .       +       .       ID=TDI_S1_0002-RA:2;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   exon    1090815 1090950 .       +       .       ID=TDI_S1_0002-RA:3;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   exon    1091479 1091538 .       +       .       ID=TDI_S1_0002-RA:4;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   exon    1095554 1095613 .       +       .       ID=TDI_S1_0002-RA:5;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   exon    1097976 1098116 .       +       .       ID=TDI_S1_0002-RA:6;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   CDS     1086430 1086437 .       +       0       ID=TDI_S1_0002-RA:cds;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   CDS     1087724 1087864 .       +       1       ID=TDI_S1_0002-RA:cds;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   CDS     1090815 1090950 .       +       1       ID=TDI_S1_0002-RA:cds;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   CDS     1091479 1091538 .       +       0       ID=TDI_S1_0002-RA:cds;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   CDS     1095554 1095613 .       +       0       ID=TDI_S1_0002-RA:cds;Parent=TDI_S1_0002-RA
# Tdi_LRv5a_scf1  maker   CDS     1097976 1098116 .       +       0       ID=TDI_S1_0002-RA:cds;Parent=TDI_S1_0002-RA
# 




# print("\n\n\nFinished, Shadow\n\n\n")
