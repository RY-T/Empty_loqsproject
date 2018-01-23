import os
from pathlib import Path

rule all:
	input:
		'STRUCT_RNA/folder_flag.done',
		'STRUCT_RNA/flags/fastqc.done',
		'STRUCT_RNA/flags/to_collapsed_tab.done',
		'STRUCT_RNA/flags/tab_to_fasta.done',
		'STRUCT_RNA/flags/artifacts_filter',
		'STRUCT_RNA/flags/Sequential_mapping.done',
		'STRUCT_RNA/flags/RNA_fold_filter',
		'STRUCT_RNA/flags/Count_Sequential_mapping.done',
		'spike_in_count.done',
		'STRUCT_RNA/output/Counts/count.txt',
		'STRUCT_RNA/flags/Count_to_graph1.done',
		'STRUCT_RNA/output/Counts/count1.png',
		'STRUCT_RNA/flags/Count_to_graph2.done',
		'STRUCT_RNA/output/Counts/count2.png',
		'STRUCT_RNA/flags/folders_2.done',
		'STRUCT_RNA/flags/all_mapping.done',
		'STRUCT_RNA/flags/Count_all_mapping.done',
		'STRUCT_RNA/flags/all_Count_to_graph1.done',
		'STRUCT_RNA/flags/all_Count_to_graph2.done'




rule make_folders:
	output:
		touch('STRUCT_RNA/folder_flag.done')
	run:
		import os
		from pathlib import Path
		STRUCT_RNA_output = Path('STRUCT_RNA/output')
		if not STRUCT_RNA_output.is_dir():
			os.mkdir('STRUCT_RNA/output')
		if not Path('STRUCT_RNA/output/Fastqc').is_dir():
			os.mkdir('STRUCT_RNA/output/Fastqc')
		if not Path('STRUCT_RNA/output/Preprocessed').is_dir():
			os.mkdir('STRUCT_RNA/output/Preprocessed')
		if not Path('STRUCT_RNA/output/Seq_map').is_dir():
			os.mkdir('STRUCT_RNA/output/Seq_map')
		if not Path('STRUCT_RNA/output/Counts').is_dir():
			os.mkdir('STRUCT_RNA/output/Counts')
		STRUCT_RNA_flags = Path('STRUCT_RNA/flags')
		if not STRUCT_RNA_flags.is_dir():
			os.mkdir('STRUCT_RNA/flags')

rule running_fastqc:
	input:
		'STRUCT_RNA/folder_flag.done'
	output:
		touch('STRUCT_RNA/flags/fastqc.done')
	run:
		import os
		fastq_folder=os.path.dirname('Input_files/small_RNA_Libs/')
		output_folder='STRUCT_RNA/output/Fastqc/'
		fastqc_path='FastQC/fastqc'
		for file in os.listdir(os.path.join(os.getcwd(),fastq_folder)):
			if file.endswith(".fastq"):
				shell('{fastqc_path} -o {output_folder} {fastq_folder}/{file}')

rule cutting_adapters_fastq_to_fasta_to_collapsed_reads_then_to_tab_format:
	input:
		'STRUCT_RNA/folder_flag.done'
	output:
		touch('STRUCT_RNA/flags/to_collapsed_tab.done')
	run:
		import os
		files_list =[]
		fastq_folder=os.path.dirname('Input_files/small_RNA_Libs/')
		output_folder='STRUCT_RNA/output/'

		adapter = 'TGGAATTCTCGGGTGCCAAGG'

		for file in os.listdir(os.path.join(os.getcwd(),fastq_folder)):
			if file.endswith(".fastq"):
				files_list.append(file)
		for fastq in files_list:
			LIB=fastq[:-6]
			shell('fastx_clipper -a {adapter} -c -v -i {fastq_folder}/{fastq}| fastq_to_fasta -v -r |\
				sed "s/^>/>{LIB}_/"| fastx_collapser -v | fasta_formatter -t | sed "s/-/_/g ; s/_/_{LIB}_/" >\
				{output_folder}/{LIB}.tab')

rule using_pandas_to_filter_reads_and_convert_everything_back_to_fasta:
	input:
		'STRUCT_RNA/flags/to_collapsed_tab.done'
	output:
		touch('STRUCT_RNA/flags/tab_to_fasta.done')
	run:
		import os
		import re
		files_list =[]
		tab_folder=os.path.dirname('STRUCT_RNA/output/')
		output_folder='STRUCT_RNA/output/'
		for file in os.listdir(os.path.join(os.getcwd(),tab_folder)):
			if file.endswith(".tab"):
				files_list.append(file)

		import pandas as pd
		for file in files_list:
			LIB=file[:-4]
			#intermediate_file=os.path.join(tab_folder,file+'len')
			text_file = open(os.path.join(output_folder,LIB+'len.fasta'),"w")
			data = pd.read_csv(os.path.join(tab_folder,file), sep='\t', names = ["Identifier", "Read"])
			data['length'] = data['Read'].str.len()
			filtered_data = data.query('18 <= length <= 30')
			#filtered_data.to_csv(intermediate_file, sep='\t', header=None, index=False)
			
			#new_data = pd.read_csv(intermediate_file, sep='\t', names = ["Identifier", "Read", "length"])
			fasta_data = filtered_data[['Read', 'Identifier']]
			fasta_data['Read'] = '-' + fasta_data['Read'].astype(str)
			fasta_data['Identifier'] = '!' + fasta_data['Identifier'].astype(str)
			fasta_data['fasta'] = fasta_data['Identifier'] + fasta_data['Read']
			fasta = fasta_data['fasta']
			string = fasta.to_string(header=False)

			for line in string.split('\n'):
				fields = re.split('!|-',line)
				text_file.write((">{}\n{}\n".format(fields[1],fields[2])))
			text_file.close()

rule artifacts_filter:
	input:
		'STRUCT_RNA/flags/tab_to_fasta.done'
	output:
		touch('STRUCT_RNA/flags/artifacts_filter')
	run:
		import os
		files_list =[]
		fasta_folder=os.path.dirname('STRUCT_RNA/output/')
		for file in os.listdir(os.path.join(os.getcwd(),fasta_folder)):
			if file.endswith("len.fasta"):
				files_list.append(file)
		for fasta in files_list:
			LIB=fasta[:-9]
			shell('fastx_artifacts_filter -v -i {fasta_folder}/{fasta} > {fasta_folder}/Preprocessed/{LIB}_final_col_filter.fasta')
			#shell('echo "Final read identiier count:"')
			#shell('grep -c {fasta_folder}/Preprocessed/{LIB}_final_col_filter.fasta')

rule RNA_fold_filter:
	input:'STRUCT_RNA/flags/artifacts_filter'
	output:
		touch('STRUCT_RNA/flags/RNA_fold_filter')
	run:
		import os
		import pandas as pd
		letters=lambda x: ''.join(c for c in x if c.isalpha())
		def filter_mfe(fasta_line):
			value=fasta_line[-6:-1]
			if value[0]=='-':
				return float(value)<=-0.50
			else:
				return float(value[1:])<=-0.50
#		filter_mfe=lambda x: float(x[-6:-1])<=-0.50

		files_list=[]
		prepro_fasta_folder=os.path.dirname('STRUCT_RNA/output/Preprocessed/')
		for file in os.listdir(prepro_fasta_folder):
			if file.endswith(".fasta"):
				files_list.append(file)
		for fasta in files_list:
			shell('RNAfold -v --noPS {prepro_fasta_folder}/{fasta} | fasta_formatter -t > {prepro_fasta_folder}/{fasta}.temp')
			data = pd.read_table(os.path.join(prepro_fasta_folder,fasta+'.temp'), sep='\t', header =None)
			out_text_file=open(os.path.join(prepro_fasta_folder,fasta+'.structure'),"a")
			for i in range(len(data)):
				if filter_mfe(data[1][i]):
					out_text_file.write('>'+str(data[0][i])+'\n'+letters(data[1][i])+'\n')
			out_text_file.close()

rule sequential_mapping:
	input:
		'STRUCT_RNA/flags/RNA_fold_filter'
	output:
		touch('STRUCT_RNA/flags/Sequential_mapping.done')
	run:
		import os
		bt_prefix_list=[]
		bt_index_folder='Index_gen/Repeat_Masker/modified_bed/fasta/bt_indexes/'
		for file in os.listdir(os.path.join(os.getcwd(),bt_index_folder)):
			if file.endswith(".rev.1.ebwt"):
				bt_prefix_list.append(file[:-11])
		def sort_index(elem):
			return int(elem[5:])
		bt_prefix_list.sort(key=sort_index)
		files_list=[]
		prepro_fasta_folder=os.path.dirname('STRUCT_RNA/output/Preprocessed/')
		seq_map_folder=os.path.dirname('STRUCT_RNA/output/Seq_map/')
		for file in os.listdir(prepro_fasta_folder):
			if file.endswith(".fasta.structure"):
				files_list.append(file)
		for fasta in files_list:
			LIB=fasta[:-33]
			for i in range(len(bt_prefix_list)):
				Index=bt_prefix_list[i]
				Previously_unmapped=bt_prefix_list[i-1]
				if i ==0:
					shell('bowtie -f --un {seq_map_folder}/{LIB}_unmapped_{Index}.fasta --al {seq_map_folder}/{LIB}_mapped_{Index}.fasta -v 0 -p 1 -t {bt_index_folder}/{Index} {prepro_fasta_folder}/{fasta} > {seq_map_folder}/{LIB}_mapped{Index}.bowtie.txt')
				elif i ==1:
					shell('bowtie -f --un {seq_map_folder}/{LIB}_unmapped_{Index}.fasta -v 0 -p 1 -t {bt_index_folder}/{Index} {seq_map_folder}/{LIB}_mapped_Index0.fasta > {seq_map_folder}/{LIB}_mapped{Index}.bowtie.txt')
				else:
					shell('bowtie -f --un {seq_map_folder}/{LIB}_unmapped_{Index}.fasta -v 0 -p 1 -t {bt_index_folder}/{Index} {seq_map_folder}/{LIB}_unmapped_{Previously_unmapped}.fasta > {seq_map_folder}/{LIB}_mapped{Index}.bowtie.txt')

rule count_by_lib:
	input:
		'STRUCT_RNA/flags/Sequential_mapping.done'
	output:
		touch('STRUCT_RNA/flags/Count_Sequential_mapping.done'),
		'STRUCT_RNA/output/Counts/count.txt'
	run:
		import os
		import pandas as pd
		Full_list=[]
		seq_map_folder=os.path.dirname('STRUCT_RNA/output/Seq_map/')
		for file in os.listdir(seq_map_folder):
			if file.endswith(".bowtie.txt"):
				Full_list.append(file)

		out_text_file=open('STRUCT_RNA/output/Counts/count.txt',"a")
		for i in Full_list:	
			if os.path.getsize(os.path.join(seq_map_folder,i))<=0:
				continue
			else:
				everything = pd.read_csv(os.path.join(seq_map_folder,i),sep="\t+",header=None)
				def count_one(x):
					return int(x.split("_",3)[2])
				def sum_all(everything):
					counter=0
					for i in range(len(everything)):
						counter+=count_one(everything[0][i])
					return counter
				out_text_file.write('\n'+i+'\t'+str(sum_all(everything)))
		out_text_file.close()

rule run_spike_in_script:
	output:touch('spike_in_count.done')
	run:
		shell('snakemake -s spikein_normalisation_snakefile')

rule count_to_graph_R:
	input:
		'STRUCT_RNA/output/Counts/count.txt',
		'Spike_in_folder/Output/Spike_in_count.txt',
		'spike_in_count.done'
	output:
		'STRUCT_RNA/output/Counts/count1.png',
		touch('STRUCT_RNA/flags/Count_to_graph1.done')
	script:
		"Input_files/R_scripts/Count_to_graph_1.R"

rule count_to_graph_R_2:
	input:
		'STRUCT_RNA/output/Counts/count.txt',
		'Spike_in_folder/Output/Spike_in_count.txt',
		'spike_in_count.done'
	output:
		'STRUCT_RNA/output/Counts/count2.png',
		touch('STRUCT_RNA/flags/Count_to_graph2.done')
	script:
		"Input_files/R_scripts/Count_to_graph_2.R"

rule phase_2_make_folders:
	input:
		'STRUCT_RNA/flags/artifacts_filter',
		'spike_in_count.done',

	output:
		touch('STRUCT_RNA/flags/folders_2.done')
	run:
		import os
		from pathlib import Path
		all_map_output = Path('STRUCT_RNA/output/All_map')
		if not all_map_output.is_dir():
			os.mkdir('STRUCT_RNA/output/All_map')
			os.mkdir('STRUCT_RNA/output/All_map/Counts')

rule all_mapping:
	input:
		'STRUCT_RNA/flags/folders_2.done'
	output:
		touch('STRUCT_RNA/flags/all_mapping.done')
	run:
		import os
		bt_prefix_list=[]
		bt_index_folder='Index_gen/Repeat_Masker/modified_bed/fasta/bt_indexes/'
		for file in os.listdir(os.path.join(os.getcwd(),bt_index_folder)):
			if file.endswith(".rev.1.ebwt"):
				bt_prefix_list.append(file[:-11])
		def sort_index(elem):
			return int(elem[5:])
		bt_prefix_list.sort(key=sort_index)
		files_list=[]
		prepro_fasta_folder=os.path.dirname('STRUCT_RNA/output/Preprocessed/')
		all_map_folder=os.path.dirname('STRUCT_RNA/output/All_map/')
		for file in os.listdir(prepro_fasta_folder):
			if file.endswith(".fasta"):
				files_list.append(file)
		for fasta in files_list:
			LIB=fasta[:-23]
			for i in range(len(bt_prefix_list)):
				Index=bt_prefix_list[i]
				if i ==0:
					shell('bowtie -f --un {all_map_folder}/{LIB}_unmapped_{Index}.fasta --al {all_map_folder}/{LIB}_mapped_{Index}.fasta -v 0 -p 1 -t {bt_index_folder}/{Index} {prepro_fasta_folder}/{fasta}')
				else:
					shell('bowtie -f -v 0 -p 1 -t {bt_index_folder}/{Index} {all_map_folder}/{LIB}_mapped_Index0.fasta > {all_map_folder}/{LIB}_mapped{Index}.bowtie.txt')

rule count_by_lib_all_map:
	input:
		'STRUCT_RNA/flags/all_mapping.done'
	output:
		touch('STRUCT_RNA/flags/Count_all_mapping.done'),
		'STRUCT_RNA/output/All_map/Counts/count.txt'
	run:
		import os
		import pandas as pd
		Full_list=[]
		all_map_folder=os.path.dirname('STRUCT_RNA/output/All_map/')
		for file in os.listdir(all_map_folder):
			if file.endswith(".bowtie.txt"):
				Full_list.append(file)

		out_text_file=open('STRUCT_RNA/output/All_map/Counts/count.txt',"a")
		for i in Full_list:	
			if os.path.getsize(os.path.join(all_map_folder,i))<=0:
				continue
			else:
				everything = pd.read_csv(os.path.join(all_map_folder,i),sep="\t+",header=None)
				def count_one(x):
					return int(x.split("_",3)[2])
				def sum_all(everything):
					counter=0
					for i in range(len(everything)):
						counter+=count_one(everything[0][i])
					return counter
				out_text_file.write('\n'+i+'\t'+str(sum_all(everything)))
		out_text_file.close()

rule all_count_to_graph_R:
	input:
		'STRUCT_RNA/output/All_map/Counts/count.txt',
		'Spike_in_folder/Output/Spike_in_count.txt',
		'spike_in_count.done'
	output:
		'STRUCT_RNA/output/All_map/Counts/count1.png',
		touch('STRUCT_RNA/flags/all_Count_to_graph1.done')
	script:
		"Input_files/R_scripts/Count_to_graph_1.R"

rule all_count_to_graph_R_2:
	input:
		'STRUCT_RNA/output/All_map/Counts/count.txt',
		'Spike_in_folder/Output/Spike_in_count.txt',
		'spike_in_count.done'
	output:
		'STRUCT_RNA/output/All_map/Counts/count2.png',
		touch('STRUCT_RNA/flags/all_Count_to_graph2.done')
	script:
		"Input_files/R_scripts/Count_to_graph_2.R"
