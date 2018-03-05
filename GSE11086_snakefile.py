import os
from pathlib import Path

rule all:
	input:
		'GSE11086/folder_flag.done',
		'GSE11086/flags/to_fasta.done',
		'GSE11086/flags/Sequential_mapping.done',
		'GSE11086/flags/Count_Sequential_mapping.done',
		'GSE11086/output/Counts/count.txt',
		'GSE11086/flags/Count_to_graph.done'

rule make_folders:
	output:
		touch('GSE11086/folder_flag.done')
	run:
		import os
		from pathlib import Path
		GSE11086_output = Path('GSE11086/output')
		if not GSE11086_output.is_dir():
			os.mkdir('GSE11086/output')
		if not Path('GSE11086/output/Preprocessed').is_dir():
			os.mkdir('GSE11086/output/Preprocessed')
		if not Path('GSE11086/output/Seq_map').is_dir():
			os.mkdir('GSE11086/output/Seq_map')
		if not Path('GSE11086/output/Counts').is_dir():
			os.mkdir('GSE11086/output/Counts')
		GSE11086_flags = Path('GSE11086/flags')
		if not GSE11086_flags.is_dir():
			os.mkdir('GSE11086/flags')

rule download_GEO_and_split_files:
#collapse mapped reads and unmap reads into the same library, change format of csv to common fasta format
	input:
		'GSE11086/folder_flag.done'
	output:
		touch('GSE11086/flags/to_fasta.done')
	run:
		import os
		#pip install GEOparse in conda environment
		#pip install biopython in conda environment
		import GEOparse
		import pandas as pandas

		csv_folder=os.path.dirname('GSE11086/')
		gse = GEOparse.get_GEO(geo="GSE11086", destdir=csv_folder)
		
		for GSM, gsm in gse.gsms.items():	 
			text_file = open('GSE11086/output/Preprocessed/'+GSM+'.fasta',"w")
			for line in range(len(gsm.table)):
				if gsm.table.iloc[line][1]>0:
					fasta_seq,fasta_name=gsm.table.iloc[line][0],str(line+1)+'_'+GSM+'_'+str(gsm.table.iloc[line][1])
					text_file.write((">{}\n{}\n".format(fasta_name,fasta_seq)))
			text_file.close()

rule sequential_mapping:
	input:
		'GSE11086/flags/to_fasta.done'
	output:
		touch('GSE11086/flags/Sequential_mapping.done')
	run:
		import os
		bt_prefix_list=[]
		bt_index_folder='index_generation/output/bt_indexes/'
		for file in os.listdir(os.path.join(os.getcwd(),bt_index_folder)):
			if file.endswith(".rev.1.ebwt"):
				bt_prefix_list.append(file[:-11])
		def sort_index(elem):
			return int(elem[5:])
		bt_prefix_list.sort(key=sort_index)
		files_list=[]
		invalid_GSM=[]
		removal_list=[]
		prepro_fasta_folder=os.path.dirname('GSE11086/output/Preprocessed/')
		seq_map_folder=os.path.dirname('GSE11086/output/Seq_map/')
		for file in os.listdir(prepro_fasta_folder):
			if file.endswith(".fasta"):
				files_list.append(file)
		for fasta in files_list:
			GSM=fasta[:-6]
			for i in range(len(bt_prefix_list)):
				Index=bt_prefix_list[i]
				Previously_unmapped=bt_prefix_list[i-1]
				if i ==0:
					shell('bowtie -f --un {seq_map_folder}/{GSM}_unmapped_{Index}.fasta --al {seq_map_folder}/{GSM}_mapped_{Index}.fasta -v 0 --best -t {bt_index_folder}{Index} {prepro_fasta_folder}/{fasta} > {seq_map_folder}/{GSM}_mapped{Index}.bowtie.txt')
				elif i ==1:
					if os.path.getsize(os.path.join(seq_map_folder,str(GSM)+'_mappedIndex0.bowtie.txt'))==0:
						invalid_GSM.append(GSM)
						files_list.pop(fasta)
					else:
						shell('bowtie -f --un {seq_map_folder}/{GSM}_unmapped_{Index}.fasta -v 0 --best -t {bt_index_folder}{Index} {seq_map_folder}/{GSM}_mapped_Index0.fasta > {seq_map_folder}/{GSM}_mapped{Index}.bowtie.txt')
				else:
					shell('bowtie -f --un {seq_map_folder}/{GSM}_unmapped_{Index}.fasta -v 0 --best -t {bt_index_folder}{Index} {seq_map_folder}/{GSM}_unmapped_{Previously_unmapped}.fasta > {seq_map_folder}/{GSM}_mapped{Index}.bowtie.txt')
			for i in invalid_GSM:
				for file in os.listdir(os.path.join(os.getcwd(),bt_index_folder)):
					if file.startswith(i):
						removal_list.append(file)
			for unwanted_file in removal_list:
				shell('rm {seq_map_folder}/{unwanted_file}')

rule count_by_lib:
	input:
		'GSE11086/flags/Sequential_mapping.done'
	output:
		touch('GSE11086/flags/Count_Sequential_mapping.done'),
		'GSE11086/output/Counts/count.txt'
	run:
		import os
		import pandas as pd
		Full_list=[]
		seq_map_folder=os.path.dirname('GSE11086/output/Seq_map/')
		for file in os.listdir(seq_map_folder):
			if file.endswith(".bowtie.txt"):
				Full_list.append(file)

		out_text_file=open('GSE11086/output/Counts/count.txt',"a")
		for i in Full_list:	
			if os.path.getsize(os.path.join(seq_map_folder,i))<=0:
				out_text_file.write('\n'+i+'\t'+str(0))
			else:
				everything = pd.read_csv(os.path.join(seq_map_folder,i),sep="\t+",header=None,index_col=False, usecols=[0,1,2,3,4,5,6])
				def count_one(x):
					return int(x.split("_",3)[2])
				def sum_all(everything):
					counter=0
					for i in range(len(everything)):
						counter+=count_one(everything[0][i])
					return counter
				out_text_file.write('\n'+i+'\t'+str(sum_all(everything)))
		out_text_file.close()

rule count_unmapped_by_lib:
	input:
		'GSE11086/flags/Sequential_mapping.done'
	output:
		touch('GSE11086/flags/Count_Sequential_mapping_unmapped.done'),
		'GSE11086/output/Counts/unmapped_count.txt'
	run:
		import os
		import pandas as pd
		Full_list=[]
		seq_map_folder=os.path.dirname('GSE11086/output/Seq_map/')
		for file in os.listdir(seq_map_folder):
			if file.endswith("unmapped_Index0.fasta"):
				Full_list.append(file)

		out_text_file=open('GSE11086/output/Counts/unmapped_count.txt',"a")
		for i in Full_list:	
			if os.path.getsize(os.path.join(seq_map_folder,i))<=0:
				continue
			else:
				everything = pd.read_csv(os.path.join(seq_map_folder,i),sep="\t+",header=None)
				def count_one(x):
					return int(x.split("_",3)[2])
				def sum_all(everything):
					counter=0
					for i in range(0,len(everything),2):
						counter+=count_one(everything[0][i])
					return counter
				out_text_file.write('\n'+i+'\t'+str(sum_all(everything)))
		out_text_file.write('\n')
		out_text_file.close()

rule count_to_graph_R_GSE11086:
	input:
		'GSE11086/output/Counts/count.txt',
		'GSE11086/output/Counts/unmapped_count.txt',
		'GSE11086/flags/Count_Sequential_mapping.done',
		'GSE11086/flags/Count_Sequential_mapping_unmapped.done'
	output:
		'GSE11086/output/Counts/count.png',
		'GSE11086/output/Counts/Reads_mapped.png',
		'GSE11086/output/Counts/percentage_per_lib.png',
		touch('GSE11086/flags/Count_to_graph.done')
	script:
		"Input_files/R_scripts/Universal_count_to_graph_ggplot.R"