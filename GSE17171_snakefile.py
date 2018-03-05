import os
from pathlib import Path

rule all:
	input:
		'GSE17171/folder_flag.done',
		'GSE17171/flags/to_fasta.done',
		'GSE17171/flags/Sequential_mapping.done',
		'GSE17171/flags/Count_Sequential_mapping.done',
		'GSE17171/output/Counts/count.txt',
		'GSE17171/flags/Count_to_graph.done'

rule make_folders:
	output:
		touch('GSE17171/folder_flag.done')
	run:
		import os
		from pathlib import Path
		GSE17171_output = Path('GSE17171/output')
		if not GSE17171_output.is_dir():
			os.mkdir('GSE17171/output')
		if not Path('GSE17171/output/Preprocessed').is_dir():
			os.mkdir('GSE17171/output/Preprocessed')
		if not Path('GSE17171/output/Seq_map').is_dir():
			os.mkdir('GSE17171/output/Seq_map')
		if not Path('GSE17171/output/Counts').is_dir():
			os.mkdir('GSE17171/output/Counts')
		GSE17171_flags = Path('GSE17171/flags')
		if not GSE17171_flags.is_dir():
			os.mkdir('GSE17171/flags')


rule csv_to_fasta:
#collapse mapped reads and unmap reads into the same library, change format of csv to common fasta format
	input:
		'GSE17171/folder_flag.done'
	output:
		touch('GSE17171/flags/to_fasta.done')
	run:
		import os
		import pandas as pd

		files_list =[]
		lib_list=[]
		csv_folder=os.path.dirname('GSE17171/')
		for file in os.listdir(os.path.join(os.getcwd(),csv_folder)):
			if file.endswith(".csv"):
				files_list.append(file)
		for i in files_list:
			if not i[:9] in lib_list:
				lib_list.append(i[:9])
			else:
				continue

		for i in lib_list:		
			text_file = open('GSE17171/output/Preprocessed/'+i+'.fasta',"w")
			counter=0
			for j in files_list:
				if j[:9]==i:
					data = pd.read_csv(os.path.join(csv_folder,j), sep='\,')
					for line in range(len(data)):
						if counter==0:
							fasta_name=str(line+1)+'_'+i+'_'+str(data['READS'].iloc[line])
						else:
							fasta_name=str(line+counter+1)+'_'+i+'_'+str(data['READS'].iloc[line])
						fasta_seq=str(data['SEQ'].iloc[line])
						if len(fasta_seq)>=18:
							text_file.write((">{}\n{}\n".format(fasta_name,fasta_seq)))
					counter=len(data)
				else:
					counter=0
			text_file.close()

rule sequential_mapping:
	input:
		'GSE17171/flags/to_fasta.done'
	output:
		touch('GSE17171/flags/Sequential_mapping.done')
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
		prepro_fasta_folder=os.path.dirname('GSE17171/output/Preprocessed/')
		seq_map_folder=os.path.dirname('GSE17171/output/Seq_map/')
		for file in os.listdir(prepro_fasta_folder):
			if file.endswith(".fasta"):
				files_list.append(file)
		for fasta in files_list:
			GSM=fasta[:9]
			for i in range(len(bt_prefix_list)):
				Index=bt_prefix_list[i]
				Previously_unmapped=bt_prefix_list[i-1]
				if i ==0:
					shell('bowtie -f --un {seq_map_folder}/{GSM}_unmapped_{Index}.fasta --al {seq_map_folder}/{GSM}_mapped_{Index}.fasta -v 2 --best -t {bt_index_folder}/{Index} {prepro_fasta_folder}/{fasta} > {seq_map_folder}/{GSM}_mapped{Index}.bowtie.txt')
				elif i ==1:
					shell('bowtie -f --un {seq_map_folder}/{GSM}_unmapped_{Index}.fasta -v 2 --best -t {bt_index_folder}/{Index} {seq_map_folder}/{GSM}_mapped_Index0.fasta > {seq_map_folder}/{GSM}_mapped{Index}.bowtie.txt')
				else:
					shell('bowtie -f --un {seq_map_folder}/{GSM}_unmapped_{Index}.fasta -v 2 --best -t {bt_index_folder}/{Index} {seq_map_folder}/{GSM}_unmapped_{Previously_unmapped}.fasta > {seq_map_folder}/{GSM}_mapped{Index}.bowtie.txt')

rule count_by_lib:
	input:
		'GSE17171/flags/Sequential_mapping.done'
	output:
		touch('GSE17171/flags/Count_Sequential_mapping.done'),
		'GSE17171/output/Counts/count.txt'
	run:
		import os
		import pandas as pd
		Full_list=[]
		seq_map_folder=os.path.dirname('GSE17171/output/Seq_map/')
		for file in os.listdir(seq_map_folder):
			if file.endswith(".bowtie.txt"):
				Full_list.append(file)

		out_text_file=open('GSE17171/output/Counts/count.txt',"a")
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
		'GSE17171/flags/Sequential_mapping.done'
	output:
		touch('GSE17171/flags/Count_Sequential_mapping_unmapped.done'),
		'GSE17171/output/Counts/unmapped_count.txt'
	run:
		import os
		import pandas as pd
		Full_list=[]
		seq_map_folder=os.path.dirname('GSE17171/output/Seq_map/')
		for file in os.listdir(seq_map_folder):
			if file.endswith("unmapped_Index0.fasta"):
				Full_list.append(file)

		out_text_file=open('GSE17171/output/Counts/unmapped_count.txt',"a")
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

rule count_to_graph_R_GSE17171:
	input:
		'GSE17171/output/Counts/count.txt',
		'GSE17171/output/Counts/unmapped_count.txt',
		'GSE17171/flags/Count_Sequential_mapping.done',
		'GSE17171/flags/Count_Sequential_mapping_unmapped.done'
	output:
		'GSE17171/output/Counts/count.png',
		'GSE17171/output/Counts/Reads_mapped.png',
		'GSE17171/output/Counts/percentage_per_lib.png',
		touch('GSE17171/flags/Count_to_graph.done')
	script:
		"Input_files/R_scripts/Universal_count_to_graph_ggplot.R"


#rule count_to_graph_R_2:
#	input:
#		'Seq_map_directory/count.txt',
#		'Spike_in_folder/Output/Spike_in_count.txt',
#		'spike_in_count.done'
#	output:
#		'Seq_map_directory/count2.png'
#	script:
#		"Input_files/R_scripts/Count_to_graph_2.R"

# df1 = data[data.index % 3 != 2]
#fasta_formatter -t -i minitest.txt > minitest.tab
# import pandas as pd
# data = pd.read_table('minitest.tab', sep='\t', header =None)
#letters=lambda x: ''.join(c for c in x if c.isalpha())
#filter_mfe=lambda x: float(x[-6:-1])>-0.50

#out_text_file=open('test_out.txt',"a")
#for i in range(len(data)):
#	if filter_mfe(data[1][i]):
#		out_text_file.write('>'+str(data[0][i])+'\n'+letters(data[1][i])+'\n')
#out_text_file.close()

