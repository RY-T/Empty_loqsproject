import os
from pathlib import Path

rule all:
	input:
		'RNA_LIB/folder_flag.done',
		'RNA_LIB/flags/renaming.done',
		'RNA_LIB/flags/fastqc.done',
		'RNA_LIB/flags/to_collapsed_tab.done',
		'RNA_LIB/flags/tab_to_fasta.done',
		'RNA_LIB/flags/artifacts_filter',
		'RNA_LIB/flags/Sequential_mapping.done',
		'RNA_LIB/flags/Count_Sequential_mapping.done',
		'spike_in_count.done',
		'RNA_LIB/output/Counts/count.txt',
		'RNA_LIB/flags/Count_to_graph1.done',
		'RNA_LIB/output/Counts/count1.png',
		'RNA_LIB/flags/Count_to_graph2.done',
		'RNA_LIB/output/Counts/count2.png'


rule make_folders:
	output:
		touch('RNA_LIB/folder_flag.done')
	run:
		import os
		from pathlib import Path
		RNA_LIB_output = Path('RNA_LIB/output')
		if not RNA_LIB_output.is_dir():
			os.mkdir('RNA_LIB/output')
		if not Path('RNA_LIB/output/Fastqc').is_dir():
			os.mkdir('RNA_LIB/output/Fastqc')
		if not Path('RNA_LIB/output/Preprocessed').is_dir():
			os.mkdir('RNA_LIB/output/Preprocessed')
		if not Path('RNA_LIB/output/Seq_map').is_dir():
			os.mkdir('RNA_LIB/output/Seq_map')
		if not Path('RNA_LIB/output/Counts').is_dir():
			os.mkdir('RNA_LIB/output/Counts')
		RNA_LIB_flags = Path('RNA_LIB/flags')
		if not RNA_LIB_flags.is_dir():
			os.mkdir('RNA_LIB/flags')

rule running_fastqc:
	input:
		'RNA_LIB/folder_flag.done'
	output:
		touch('RNA_LIB/flags/fastqc.done')
	run:
		import os
		fastq_folder=os.path.dirname('Input_files/small_RNA_Libs/')
		output_folder='RNA_LIB/output/Fastqc/'
		fastqc_path='RNA_LIB/FastQC/fastqc'
		for file in os.listdir(os.path.join(os.getcwd(),fastq_folder)):
			if file.endswith(".fastq"):
				shell('{fastqc_path} -o {output_folder} {fastq_folder}/{file}')

rule cutting_adapters_fastq_to_fasta_to_collapsed_reads_then_to_tab_format:
	input:
		'RNA_LIB/folder_flag.done'
	output:
		touch('RNA_LIB/flags/to_collapsed_tab.done')
	run:
		import os
		files_list =[]
		fastq_folder=os.path.dirname('Input_files/small_RNA_Libs/')
		output_folder='RNA_LIB/output/'

		adapter = 'TGGAATTCTCGGGTGCCAAGG'

		for file in os.listdir(os.path.join(os.getcwd(),fastq_folder)):
			if file.endswith(".fastq"):
				files_list.append(file)
		for fastq in files_list:
			LIB=fastq[:6]
			shell('fastx_clipper -a {adapter} -c -v -i {fastq_folder}/{fastq}| fastq_to_fasta -v -r |\
				sed "s/^>/>{LIB}_/"| fastx_collapser -v | fasta_formatter -t | sed "s/-/_/g ; s/_/_{LIB}_/" >\
				{output_folder}/{LIB}.tab')

rule using_pandas_to_filter_reads_and_convert_everything_back_to_fasta:
	input:
		'RNA_LIB/flags/to_collapsed_tab.done'
	output:
		touch('RNA_LIB/flags/tab_to_fasta.done')
	run:
		import os
		import re
		files_list =[]
		tab_folder=os.path.dirname('RNA_LIB/output/')
		output_folder='RNA_LIB/output/'
		for file in os.listdir(os.path.join(os.getcwd(),tab_folder)):
			if file.endswith(".tab"):
				files_list.append(file)

		import pandas as pd
		for file in files_list:
			LIB=file[:4]
			text_file = open(output_folder/LIB+'len.fasta',"w")
			data = pd.read_csv(tab_folder/file, sep='\t', names = ["Identifier", "Read"])
			data['length'] = data['Read'].str.len()
			filtered_data = data.query('18 <= length <= 30')
			
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
		'RNA_LIB/flags/tab_to_fasta.done'
	output:
		touch('RNA_LIB/flags/artifacts_filter')
	run:
		import os
		files_list =[]
		fasta_folder=os.path.dirname('RNA_LIB/output/')
		for file in os.listdir(os.path.join(os.getcwd(),fasta_folder)):
			if file.endswith("len.fasta"):
				files_list.append(file)
		for fasta in files_list:
			LIB=fasta[:-9]
			shell('fastx_artifacts_filter -v -i {fasta_folder}/{fasta} > {fasta_folder}/Preprocessed/{LIB}_final_col_filter.fasta\
				echo "Final read identiier count:"\
				grep -c {fasta_folder}/Preprocessed/{LIB}_final_col_filter.fasta')

rule sequential_mapping:
	input:
		'RNA_LIB/flags/artifacts_filter'
	output:
		touch('RNA_LIB/flags/Sequential_mapping.done')
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
		prepro_fasta_folder=os.path.dirname('RNA_LIB/output/Preprocessed/')
		seq_map_folder=os.path.dirname('RNA_LIB/output/Seq_map/')
		for file in os.listdir(prepro_fasta_folder):
			if file.endswith(".fasta"):
				files_list.append(file)
		for fasta in files_list:
			LIB=fasta[:-23]
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
		'RNA_LIB/flags/Sequential_mapping.done'
	output:
		touch('RNA_LIB/flags/Count_Sequential_mapping.done'),
		'RNA_LIB/output/Counts/count.txt'
	run:
		import os
		import pandas as pd
		Full_list=[]
		seq_map_folder=os.path.dirname('RNA_LIB/output/Seq_map/')
		for file in os.listdir(seq_map_folder):
			if file.endswith(".bowtie.txt"):
				Full_list.append(file)

		out_text_file=open('RNA_LIB/output/Counts/count.txt',"a")
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
		'RNA_LIB/output/Counts/count.txt',
		'Spike_in_folder/Output/Spike_in_count.txt',
		'spike_in_count.done'
	output:
		'RNA_LIB/output/Counts/count1.png',
		touch('RNA_LIB/flags/Count_to_graph1.done')
	script:
		"Input_files/R_scripts/Count_to_graph_1.R"

rule count_to_graph_R_2:
	input:
		'RNA_LIB/output/Counts/count.txt',
		'Spike_in_folder/Output/Spike_in_count.txt',
		'spike_in_count.done'
	output:
		'RNA_LIB/output/Counts/count2.png',
		touch('RNA_LIB/flags/Count_to_graph2.done')
	script:
		"Input_files/R_scripts/Count_to_graph_2.R"