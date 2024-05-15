import sys
import os
import getopt
import decimal
from decimal import *
import re
import collections

try:
	opts, args = getopt.getopt(sys.argv[1:], 'm:g:b:p:o:Ph')
																						
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)

missing_busco_list_filename = None
genome_busco_dir            = None
genome_prefix 				= None 
gff_file_name               = None
out_file_prefix             = "test_BUSCO"
add_phase  					= False

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h'):
		print("\n**** Add_missing_busco.py | Written by DJP, 17/06/22 in Python 3.5 in Porto ****\n")
		print("Takes genome BUSCO output and adds genes it finds missing to the an annotation file")	
		print("\n**** Usage****\n")
		print("python3 Add_missing_busco.py -m [missing busco genes table] -g [gff annotation] -b [genome busco directory] -p [genome prefix] -o [output prefix] \n")
		print("\n**** OPTIONS ****\n")
		print("-P\tSpecify to add phase to added CDS records [default = OFF]\n")
		print("\n**** test example - see 2b_run_annot_T**.sh for real usage ****\n")
		print("python3 Add_missing_busco.py -m BUSCO_play/Tps_braker_prot_run_3e_RNAseq_run_4e_combined_3_FBGgiUGF1000M0_genes_BUSCO/run_insecta_odb10/missing_busco_list.tsv -b BUSCO_play/Tps_LRv5b/ -g /Users/dparker/Desktop/BUSCO_play/Tps_braker_prot_run_3e_RNAseq_run_4e_combined_3_FBGgiUGF1000M0.gff\n\n")
				
		sys.exit(2)
	elif opt in ('-m'):
		missing_busco_list_filename = arg
	elif opt in ('-g'):
		gff_file_name = arg
	elif opt in ('-b'):
		genome_busco_dir = arg
	elif opt in ('-p'):
		genome_prefix = arg
	elif opt in ('-o'):
		out_file_prefix  = arg
	elif opt in ('-P'):
		add_phase  = True	
	else:
		print("i dont know")
		sys.exit(2)


### get missing BUSCO list

missing_busco_list_file = open(missing_busco_list_filename)
missing_BUSCO = set()
for line in missing_busco_list_file :
	line = line.strip()
	if not line.startswith("#"):
		missing_BUSCO.add(line)
print("N missing BUSCO genes:")
print(len(missing_BUSCO))

###### get gff gene names

gene_N_max = -1
gene_prefix = None
if genome_prefix != None: 
	gene_prefix = genome_prefix
else:
	print("Warning, NO genome prefix specified. Guessing it is: ")

N_geneN_digits = None
gff_file = open(gff_file_name)
for line in gff_file:
	line = line.strip()
	if not line.startswith("#"):
		line = line.split("\t")
		feat = line[2]
		#print(line)
		if feat == "gene":
			gene_id = line[8].split("ID=")[1].split(";")[0]
			N_geneN_digits = len(gene_id.split("_")[1])
			if genome_prefix != None:
				gene_N  = int(gene_id.split(gene_prefix)[1].strip("_"))
			else:
				
				gene_N  = int(gene_id.split("_")[1])
				gene_prefix = gene_id.split("_")[0]
			
			# print(gene_id)
			# print(line)
			if gene_N > gene_N_max:
				gene_N_max = gene_N


if genome_prefix == None: 
	print(gene_prefix)
				
print(gene_N_max)
gff_file.close()

### find missing busco in genome BUSCO

to_add_busco_dir_coords = set()

full_genome_table = open(os.path.join(genome_busco_dir, "run_insecta_odb10/full_table.tsv"))
for line in full_genome_table:
	line = line.strip().split("\t")
	if line[0] in missing_BUSCO:
		if line[1] == "Complete":
			#print(line)
			#print([line[2], int(line[3]) + 1])
			to_add_busco_dir_coords.add((line[2], str(int(line[3]) + 1), str(int(line[4]) + 1), line[5])) ### DJP: checked the gff has the correct coords!

print("Number of complete busco found in genome busco run:")
print(len(to_add_busco_dir_coords))
# print(to_add_busco_dir_coords)


################ gffs contain same gene names pointing to different genes! How lovely. Need to filter out the 'genes' i don't want (these look like small genes... )
# grep "TCS_ID=145210at50557|Tps_LRv5b_scf8|+" BUSCO_play/Tps_LRv5b/run_insecta_odb10/metaeuk_output/rerun_results/Tps_LRv5b.fasta.gff 
# Tps_LRv5b_scf8	gene	MetaEuk	8556569	8580652	182	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+
# Tps_LRv5b_scf8	mRNA	MetaEuk	8556569	8580652	182	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_mRNA;Parent=145210at50557|Tps_LRv5b_scf8|+
# Tps_LRv5b_scf8	exon	MetaEuk	8556569	8556685	31	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_exon;Parent=145210at50557|Tps_LRv5b_scf8|+_mRNA
# Tps_LRv5b_scf8	CDS	MetaEuk	8556569	8556685	31	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_CDS;Parent=145210at50557|Tps_LRv5b_scf8|+_exon
# Tps_LRv5b_scf8	exon	MetaEuk	8560412	8560498	44	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_exon;Parent=145210at50557|Tps_LRv5b_scf8|+_mRNA
# Tps_LRv5b_scf8	CDS	MetaEuk	8560412	8560498	44	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_CDS;Parent=145210at50557|Tps_LRv5b_scf8|+_exon
# Tps_LRv5b_scf8	exon	MetaEuk	8576891	8577007	52	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_exon;Parent=145210at50557|Tps_LRv5b_scf8|+_mRNA
# Tps_LRv5b_scf8	CDS	MetaEuk	8576891	8577007	52	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_CDS;Parent=145210at50557|Tps_LRv5b_scf8|+_exon
# Tps_LRv5b_scf8	exon	MetaEuk	8580533	8580652	52	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_exon;Parent=145210at50557|Tps_LRv5b_scf8|+_mRNA
# Tps_LRv5b_scf8	CDS	MetaEuk	8580536	8580652	52	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_CDS;Parent=145210at50557|Tps_LRv5b_scf8|+_exon
# Tps_LRv5b_scf8	gene	MetaEuk	56330405	56366479	46	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+
# Tps_LRv5b_scf8	mRNA	MetaEuk	56330405	56366479	46	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_mRNA;Parent=145210at50557|Tps_LRv5b_scf8|+
# Tps_LRv5b_scf8	exon	MetaEuk	56330405	56330449	21	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_exon;Parent=145210at50557|Tps_LRv5b_scf8|+_mRNA
# Tps_LRv5b_scf8	CDS	MetaEuk	56330405	56330449	21	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_CDS;Parent=145210at50557|Tps_LRv5b_scf8|+_exon
# Tps_LRv5b_scf8	exon	MetaEuk	56366402	56366479	24	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_exon;Parent=145210at50557|Tps_LRv5b_scf8|+_mRNA
# Tps_LRv5b_scf8	CDS	MetaEuk	56366402	56366479	24	+	.	Target_ID=145210at50557;TCS_ID=145210at50557|Tps_LRv5b_scf8|+_CDS;Parent=145210at50557|Tps_LRv5b_scf8|+_exon


gff_1_path = os.path.join(genome_busco_dir, "run_insecta_odb10/metaeuk_output/initial_results/")
for gff_1_path, subdirs, files in os.walk(gff_1_path):
	for name in files:
		if name.endswith(".fasta.gff"):
			print(name)
			gff_1_file = open(os.path.join(gff_1_path, name))
	
gff_2_path = os.path.join(genome_busco_dir, "run_insecta_odb10/metaeuk_output/rerun_results/")
for gff_2_path, subdirs, files in os.walk(gff_2_path):
	for name in files:
		if name.endswith(".fasta.gff"):
			print(name)
			gff_2_file = open(os.path.join(gff_2_path, name))

#gff_1_file = open(os.path.join(genome_busco_dir, "run_insecta_odb10/metaeuk_output/initial_results/*.fasta.gff"))
#gff_2_file = open(os.path.join(genome_busco_dir, "run_insecta_odb10/metaeuk_output/rerun_results/*.fasta.gff"))
temp_file_0 =  open(out_file_prefix + "_AMB_gff0.temp", "w")


for line in gff_1_file:
	line = line.strip()
	feat = line.split("\t")[1]
	if feat == "gene":
		temp_file_0.write("\n" + line)
	else:
		temp_file_0.write("____SPLITDJP_____" + line)

for line in gff_2_file:
	line = line.strip()
	feat = line.split("\t")[1]
	if feat == "gene":
		temp_file_0.write("\n" + line)
	else:
		temp_file_0.write("____SPLITDJP_____" + line)


temp_file_0.close()

temp_file_0 = open(out_file_prefix + "_AMB_gff0.temp")
temp_file_1 = open(out_file_prefix + "_AMB_gff1.temp", "w")
for line in temp_file_0:
	line_o = line
	line = line.strip().split("\t")
	if len(line) > 1:
		coord_c = ((line[0], line[3], line[4],line[6]))
		#print(coord_c)
		if coord_c in to_add_busco_dir_coords:
			line_o = line_o.strip().replace("____SPLITDJP_____", "\n")
			temp_file_1.write(line_o + "\n")

temp_file_1.close()

#####  make gff3 format

temp_file_1 = open(out_file_prefix + "_AMB_gff1.temp")
temp_file_2 = open(out_file_prefix + "_AMB_gff2.temp", "w")
found_gene_ids = set()
seen_exon_trans = set()
new_trans_to_exon_dict = {}
seen_CDS_trans = set()
trans_to_CDS_dict = {}

n_gene = 0
n_mRNA = 0
gene_N = gene_N_max 
gene_to_newgene_name_dict = {}
trans_orientation_dict =  {}
line_N = 0
for line in temp_file_1:
	line = line.strip().split("\t")
	feat = line[1]
	orient = line[6]
	gene_id = None
	line_N = line_N + 1
	new_desc_line = ""
	if feat == "gene":
		gene_id  = line[8].split("TCS_ID=")[1].strip()
		gene_N = gene_N + 1
		new_gene_name = gene_prefix + "_" + str(gene_N).zfill(N_geneN_digits)
		gene_to_newgene_name_dict[gene_id] = new_gene_name
		n_gene = n_gene + 1
		new_desc_line = "ID=" + new_gene_name + ";Alias=" + gene_id + ";Name=" + new_gene_name
	elif feat == "mRNA":
		gene_id  = line[8].split("TCS_ID=")[1].strip().split("Parent=")[1].strip()
		new_gene_name = gene_to_newgene_name_dict.get(gene_id)
		n_mRNA = n_mRNA + 1
		new_desc_line = "ID=" + new_gene_name + "-RA" + ";Alias=" + gene_id + "_mRNA" + ";Parent=" + new_gene_name + ";Name=" + new_gene_name + "-RA"
		trans_orientation_dict[new_gene_name + "-RA"] = orient 
		
	elif feat == "exon":
		gene_id  = line[8].split("TCS_ID=")[1].strip().split("Parent=")[1].strip().split("_exon")[0].split("_mRNA")[0]
		new_gene_name = gene_to_newgene_name_dict.get(gene_id)
		new_desc_line = "ID=" + new_gene_name + "-RA" + ":" + feat + str(line_N) +  ";" + "Parent=" + new_gene_name + "-RA"
		new_trans_id = new_gene_name + "-RA"
		# print(gene_id)
		# print(new_gene_name)
		# print(new_trans_id)
		if new_trans_id not in seen_exon_trans:
			seen_exon_trans.add(new_trans_id)
			new_trans_to_exon_dict[new_trans_id] = [new_trans_id + ":" + feat + str(line_N)]
		else:
			rec = new_trans_to_exon_dict.get(new_trans_id)
			rec.append(new_trans_id + ":" + feat + str(line_N))
			new_trans_to_exon_dict[new_trans_id] = rec
	elif  feat == "CDS":
		gene_id  = line[8].split("TCS_ID=")[1].strip().split("Parent=")[1].strip().split("_exon")[0].split("_mRNA")[0]
		new_gene_name = gene_to_newgene_name_dict.get(gene_id)
		new_desc_line = "ID=" + new_gene_name + "-RA" + ":" + feat +  ";" + "Parent=" + new_gene_name + "-RA"
		if new_trans_id not in seen_CDS_trans:
			seen_CDS_trans.add(new_trans_id)
			trans_to_CDS_dict[new_trans_id] = [(new_trans_id + ":" + feat, int(line[3]), int(line[4]))]
		else:
			rec = trans_to_CDS_dict.get(new_trans_id)
			rec.append((new_trans_id + ":" + feat, int(line[3]), int(line[4])))
			trans_to_CDS_dict[new_trans_id] = rec

	
	else:
		gene_id  = line[8].split("TCS_ID=")[1].strip().split("Parent=")[1].strip().split("_exon")[0].split("_mRNA")[0]
		new_gene_name = gene_to_newgene_name_dict.get(gene_id)
		new_desc_line = "ID=" + new_gene_name + "-RA" + ":" + feat +  ";" + "Parent=" + new_gene_name + "-RA"		
	#print(new_desc_line)
	
	temp_file_2.write(line[0] + "\t" + line[2] + "\t" + line[1] + "\t" + line[3] + "\t" + line[4] + "\t" + "." + "\t" +  line[6]  + "\t" + "." + "\t" + new_desc_line  + "\n")

temp_file_2.close()

if n_mRNA != n_gene:
	print("more than one mRNA per gene. This code can't deal with that. Exiting.")
	sys.exit(2)

print("N exon mRNA")
print(len(seen_exon_trans))

##########################################################################################3
### renumber exons

old_exon_new_exon_name_dict = {}

### gff  is ordered within respect to strand...
for tra in new_trans_to_exon_dict:
	exons = new_trans_to_exon_dict.get(tra)
	orient = trans_orientation_dict.get(tra)
	e_N = 0
	for e in exons:
		e_N = e_N + 1
		old_exon_new_exon_name_dict[e] = tra + ":exon" + str(e_N)


temp_file_2 = open(out_file_prefix + "_AMB_gff2.temp")
temp_file_3 = open(out_file_prefix + "_AMB_gff3.temp", "w")

for line in temp_file_2:
	line = line.strip().split("\t")
	feature = line[2]
	descrip_line = line[8]
	if feature == "exon":
		exon_ID   = descrip_line.split("ID=")[1].split(";")[0]
		parent_ID = descrip_line.split("Parent=")[1].split(";")[0]
		new_exon_ID = old_exon_new_exon_name_dict.get(exon_ID)
		new_description = "ID=" + new_exon_ID +  ";" + "Parent=" + parent_ID
		temp_file_3.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + new_description + "\n")
	else:
		temp_file_3.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\n")

temp_file_3.close()
### workout CDS phase





#####################################################################################################
### option - no need for it as gffread can handle it without
### also all phase == 0 (i think) - pretty sure this is because metaEuk don't bother working out the splice boundries nicely (this is also OK for mapping RNAseq too)

if add_phase == True:
	def CDS_phase_finder(cds_list, ori):
		cds_list_s = None
		out_cds_list = []
		if ori == "+":
			cds_list_s = sorted(cds_list, key=lambda x: x[1])
			#print("++++++")
		elif ori == "-":
			#print("-----")
			cds_list_s = sorted(cds_list, key=lambda x: x[1], reverse = True)
		else:
			cds_list_s = None
	
		remainder_list = [0]
		for i in cds_list_s:
			len_c = i[2] - i[1] + 1
			remainder = len_c % 3
			#print(len_c)
			#print(remainder)
			remainder_list.append(remainder)
		for i in range(0, len(cds_list_s)):
			cds_list_s_i = cds_list_s[i]
			phase = 3 - remainder_list[i]
			if phase == 3:
				phase = 0
			out_cds_list.append((cds_list_s_i[0], cds_list_s_i[1], cds_list_s_i[2], phase))	
		return(out_cds_list)
	
	
	exon_phase_dict = {}
	
	for el in trans_to_CDS_dict:
		rec = trans_to_CDS_dict.get(el)
		orient = trans_orientation_dict.get(el)
		phase_list = CDS_phase_finder(rec, orient)
		for i in phase_list:
			exon_phase_dict[(i[0],i[1],i[2])] = i[3]
	
	temp_file_3 = open(out_file_prefix + "_AMB_gff3.temp")
	temp_file_4 = open(out_file_prefix + "_AMB_gff4.temp", "w")
	
	for line in temp_file_3:
		line = line.strip().split("\t")
		feature = line[2]
		descrip_line = line[8]
		#print(descrip_line)
		if feature == "CDS":
			CDS_id = descrip_line.split("ID=")[1].strip().split(";")[0]
			#print((CDS_id, line[3], line[4]))
			phase = exon_phase_dict.get((CDS_id, int(line[3]), int(line[4])))
			#print(phase)
			temp_file_4.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + str(phase) + "\t" + line[8] + "\n")
		else:
			temp_file_4.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7]    + "\t" + line[8] + "\n")
	temp_file_3.close()
	temp_file_4.close()

####### add gene_id flag to all non-gene ids

temp_file_5 = open(out_file_prefix + "_AMB_gff5.temp", "w")
busco_extra_file = None
if add_phase == True:
	busco_extra_file = open(out_file_prefix + "_AMB_gff4.temp")
else:
	busco_extra_file = open(out_file_prefix + "_AMB_gff3.temp")
	
for line in busco_extra_file:
	if not line.startswith("#"):
		line = line.strip().split("\t")
		feature = line[2]
		descrip_line = line[8]
		if feature == "gene":
			temp_file_5.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\n")
		else:
			gene_ID = descrip_line.split("Parent=")[1].split(";")[0].split("-")[0]
			temp_file_5.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + ";gene_id=" + gene_ID + "\n")
	
temp_file_5.close()
busco_extra_file.close()

	

################################################################################################3
#### add in order to gff


## join files
temp_file_5 = open(out_file_prefix + "_AMB_gff5.temp")
temp_file_6 = open(out_file_prefix + "_AMB_gff6.temp", "w")


for line in temp_file_5:
	temp_file_6.write(line)

gff_file = open(gff_file_name)
for line in gff_file :
	if not line.startswith("#"):
		temp_file_6.write(line)

temp_file_5.close()
temp_file_6.close()

### ordering of records
###### need to rank within gene


# get_gene_start_coord

gene_start_dict = {}
temp_file_6 = open(out_file_prefix + "_AMB_gff6.temp")
for line in temp_file_6:
	line_o = line
	line = line.rstrip("\n").split("\t")
	scaf    = line[0]
	start_c = int(line[3])
	feature = line[2]
	ID = line[8].split("ID=")[1].split(";")[0]	
	if feature == "gene":
		gene_ID = line[8].split("ID=")[1].split(";")[0]
		gene_start_dict[gene_ID] = start_c 
		
temp_file_6.close()



temp_file_6             = open(out_file_prefix + "_AMB_gff6.temp")
out_fixgff_file         = open(gff_file_name.replace(".gff", "_wB.gff"), "w")
out_fixgff_file.write("##gff-version 3\n")

all_record_list = []
all_rec_dict    = {}

for line in temp_file_6:
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
	#print(rec_ID)
	all_record_list.append(rec_ID)
	all_rec_dict[rec_ID] = line_o




### sort

# aa = [("scaf_1",  int(99), int(1)), ("scaf_1",  int(10000000), int(2)),  ("scaf_1",  int(10000000), int(1)), ("scaf_10",  int(1), int(2)), ("scaf_2",  int(1), int(2))]
# print(sorted(aa, key=lambda e: (e[0], e[1], e[2])))


all_record_list_s = sorted(all_record_list, key=lambda e: (e[0], e[1], e[2], e[3], e[4], e[5]))
#print(all_record_list_s )

for r in all_record_list_s:
	out_l = all_rec_dict.get(r)
	#print(out_l)
	out_fixgff_file.write(out_l)
	
print(gff_file_name.replace(".gff", "_wB.gff"))

os.remove(out_file_prefix + "_AMB_gff0.temp")
os.remove(out_file_prefix + "_AMB_gff1.temp")
os.remove(out_file_prefix + "_AMB_gff2.temp")
os.remove(out_file_prefix + "_AMB_gff3.temp")
os.remove(out_file_prefix + "_AMB_gff5.temp")
os.remove(out_file_prefix + "_AMB_gff6.temp")

if add_phase == True:
	os.remove(out_file_prefix + "_AMB_gff4.temp")

print("\nDone Gullan\n\n")
