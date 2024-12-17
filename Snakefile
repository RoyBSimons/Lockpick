configfile:     config['path_to_config_file']
workdir:	config['path_to_work_dir']
rule all:
        input:
               outconfig = config['output_directory'] + "config.json",
               multiqc_1 = config['output_directory'] + "multiqc_report.html",
               target_cpgs_cov = expand(config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bismark_target_cpgs.cov", sample=config["samples"]),
               off_target_summary = expand(config['output_directory'] + "{sample}/" + "{sample}_off_target_summary.txt", sample = config["samples"]),
               target_cpgs_cov_deduplicated = expand(config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe_barcode.deduplicated.bismark_target_cpgs.cov", sample=config["samples"])
#---------------------------------------------------------------------------------------------------

rule copy_config_to_output_directory:
        input:
                config['path_to_config_file']
        output:
                config['output_directory'] + "config.json"
        shell:
                """
		cp {input} {output}
		cp {config[path_to_scripts]}../Snakefile {config[output_directory]}Snakefile
		"""
#STEP 1: Prepare reads: A. fastqc & Multiqc, C. cut adapter and probe sequences. ### Add a filter-on-read-quality step! (B)
rule fastqc:
        input:
                config['fastq_folder'] + "{sample}_{read}.fq"
        output:
                html = config['output_directory'] + "{sample}/" + "{sample}_{read}_fastqc.html",
                zip = config['output_directory'] + "{sample}/" + "{sample}_{read}_fastqc.zip"
	shell:
               "fastqc {input} -o {config[output_directory]}"+ "{wildcards.sample}/"

rule multiqc:
        input:
                fastqcs = expand(config['output_directory'] + "{sample}/" + "{sample}_{read}_fastqc.zip", sample=config["samples"], read = ["1","2"])
        output:
                config['output_directory'] + "multiqc_report.html"
	shell:
               "multiqc {input.fastqcs} -o {config[output_directory]}"

rule create_primer_files_for_filtering:
        input:
                chosen_panel = config['chosen_panel_csv']
        output:
                primers_A_2 = config['output_directory'] + "primers_A_2.fasta",
		primers_a_1 = config['output_directory'] + "primers_a_1.fasta",
		primers_G_2 = config['output_directory'] + "primers_G_2.fasta",
		primers_g_1 = config['output_directory'] + "primers_g_1.fasta"
        shell:
                "python {config[path_to_scripts]}chosen_panel_to_trim_fastas.py -c {input.chosen_panel} -o {config[output_directory]}"

rule perform_cutadapt:
        input:
                fq_in_1 = config['fastq_folder'] + "{sample}_1.fq",
                fq_in_2 = config['fastq_folder'] + "{sample}_2.fq",
                primers_A_2 = config['output_directory'] + "primers_A_2.fasta",
                primers_a_1 = config['output_directory'] + "primers_a_1.fasta",
		primers_g_1 = config['output_directory'] + "primers_g_1.fasta",
		primers_G_2 = config['output_directory'] + "primers_G_2.fasta"
	output:
		fq_adapter_1 = config['output_directory'] + "{sample}/" + "{sample}_adapter_1.fq",
		fq_adapter_2 = config['output_directory'] + "{sample}/" + "{sample}_adapter_2.fq",
		fq_a_1 = config['output_directory'] + "{sample}/" + "{sample}_a_1.fq",
		info_a_1 = config['output_directory'] + "{sample}/" + "{sample}_info_a_1.txt",
		fq_A_2 = config['output_directory'] + "{sample}/" + "{sample}_A_2.fq",
		info_A_2 = config['output_directory'] + "{sample}/" + "{sample}_info_A_2.txt",
                fq_out_1 = config['output_directory'] + "{sample}/" + "{sample}_filtered_1.fq",
		info_g_1 = config['output_directory'] + "{sample}/" + "{sample}_info_g_1.txt",
                fq_out_2 = config['output_directory'] + "{sample}/" + "{sample}_filtered_2.fq",
		info_G_2 = config['output_directory'] + "{sample}/" + "{sample}_info_G_2.txt"
	resources:
		load=10
	shell:
		"""
		cutadapt {config[cutadapt_flags]} -o {output.fq_adapter_1} -p {output.fq_adapter_2} {input.fq_in_1} {input.fq_in_2} --cores 10
		cutadapt -a file:{input.primers_a_1} -o {output.fq_a_1}  {output.fq_adapter_1} --info-file {output.info_a_1} --action retain --cores 10
		cutadapt --cut 9 -g ^file:{input.primers_g_1} -o {output.fq_out_1} {output.fq_a_1} --info-file {output.info_g_1} --action retain --cores 10
		cutadapt -a file:{input.primers_A_2} -o {output.fq_A_2} {output.fq_adapter_2} --info-file {output.info_A_2} --action retain --cores 10
		cutadapt -g ^file:{input.primers_G_2} -o {output.fq_out_2} {output.fq_A_2} --info-file {output.info_G_2} --action retain --cores 10
		"""

#STEP 2: Align to target references
rule prepare_reference_genome_from_chosen_panel:
	input:
                chosen_panel = config['chosen_panel_csv'],
                reference_genome = config['reference_genome_path']
	output:
                targets_fasta = config['output_directory'] + "targets.fasta",
                targets_bed = config['output_directory'] + "targets.bed",
                amplicons  = config['output_directory'] + "reference_amplicons/reference_amplicons.fasta"
	shell:
               """
               python {config[path_to_scripts]}prepare_reference_genome_from_chosen_panel.py -c {input.chosen_panel} -t {output.targets_fasta} -b {output.targets_bed} -r {input.reference_genome} -a {output.amplicons}
               """

rule bismark_prepare_genome:
	input:
               amplicons  = config['output_directory'] + "reference_amplicons/reference_amplicons.fasta"
	output:
               CT_converted_genome = config['output_directory'] + "reference_amplicons/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
               GA_converted_genome = config['output_directory'] + "reference_amplicons/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
	shell:
              "bismark_genome_preparation {config[output_directory]}reference_amplicons/"
	
rule bismark:
	input:
              CT_converted_genome = config['output_directory'] + "reference_amplicons/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
              GA_converted_genome = config['output_directory'] + "reference_amplicons/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
              fq_in_1 = config['output_directory'] + "{sample}/" + "{sample}_filtered_1.fq",
              fq_in_2 = config['output_directory'] + "{sample}/" + "{sample}_filtered_2.fq"
	output:
              bam = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bam"
	resources:
	      load=10
	shell:
		"""
		bismark --genome {config[output_directory]}reference_amplicons/ {config[bismark_flags]} -1 {input.fq_in_1} -2 {input.fq_in_2} -o {config[output_directory]}{wildcards.sample}/ --temp_dir {config[output_directory]}
		mv {config[output_directory]}{wildcards.sample}/{wildcards.sample}_filtered_1_bismark_bt2_pe.bam {config[output_directory]}{wildcards.sample}/{wildcards.sample}_1_bismark_bt2_pe.bam
		"""



rule bismark_report_and_summary: #summary will combine reports per sample
	input:
             bam = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bam"
	output:
             report = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_PE_report.html",
             summary = config['output_directory'] + "{sample}/" + "{sample}_bismark_summary.html"
	shell:
             """
             bismark2report --dir {config[output_directory]}{wildcards.sample}/
             bismark2summary {input.bam} -o {config[output_directory]}{wildcards.sample}/bismark_summary
             """
rule bismark_methylation_extractor:
	input:
             bam = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bam"
	output:
             cov = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bismark.cov.gz"
	resources:
	     load=10
	shell:
             "bismark_methylation_extractor -p --cutoff 1 --gzip --multicore 3 --comprehensive --bedGraph --zero_based --cytosine_report --no_overlap --genome_folder {config[output_directory]}reference_amplicons/ -o {config[output_directory]}{wildcards.sample}/ {input.bam}"

#STEP 3 Probe panel analysis (probe distribution)
rule probe_panel_analysis:
	input:
             targeted_cov = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bismark.cov.gz"
	output:
             genomic_cov = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bismark_genomic.cov"
	shell:
             "python {config[path_to_scripts]}targeted_cov_to_genome_cov.py -t {input.targeted_cov} -g {output.genomic_cov} -o {config[output_directory]}{wildcards.sample}/"

rule keep_target_cpgs:
	input:
            genomic_cov = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bismark_genomic.cov",
            target_cpg_bed = config['target_cpg_bed']
	output:
            target_cpgs_cov = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bismark_target_cpgs.cov"
	shell:
            "bedtools intersect -a {input.target_cpg_bed} -b {input.genomic_cov} -wb > {output.target_cpgs_cov}"

rule off_target_analysis:
	input:
		info_a_1 = config['output_directory'] + "{sample}/" + "{sample}_info_a_1.txt",
		info_A_2 = config['output_directory'] + "{sample}/" + "{sample}_info_A_2.txt",
		info_g_1 = config['output_directory'] + "{sample}/" + "{sample}_info_g_1.txt",
		info_G_2 = config['output_directory'] + "{sample}/" + "{sample}_info_G_2.txt",
		bam = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bam"
	output:
		origin_a_1 = config['output_directory'] + "{sample}/" + "{sample}_origin_a_1.tsv",
		origin_A_2 = config['output_directory'] + "{sample}/" + "{sample}_origin_g_1.tsv",
		origin_g_1 = config['output_directory'] + "{sample}/" + "{sample}_origin_A_2.tsv",
		origin_G_2 = config['output_directory'] + "{sample}/" + "{sample}_origin_G_2.tsv",
		origin = config['output_directory'] + "{sample}/" + "{sample}_origin.tsv",
		insert_lengths = config['output_directory'] + "{sample}/" + "{sample}_insert_lengths.tsv",
		read_id = config['output_directory'] + "{sample}/" + "{sample}_read_id.tsv",
		read_id_mapped = config['output_directory'] + "{sample}/" + "{sample}_read_id_mapped.tsv",
		off_target_summary = config['output_directory'] + "{sample}/" + "{sample}_off_target_summary.txt",
		read_groups =  config['output_directory'] + "{sample}/" + "{sample}_read_groups.csv"

	shell:
		"""
		cat {input.info_a_1} | awk -F '\t' '{{print $8}}' > {output.origin_a_1}
		cat {input.info_g_1} | awk -F '\t' '{{print $8}}' > {output.origin_g_1}
		cat {input.info_A_2} | awk -F '\t' '{{print $8}}' > {output.origin_A_2}
		cat {input.info_G_2} | awk -F '\t' '{{print $8}}' > {output.origin_G_2}
		cat {input.info_G_2} | awk -F '\t' '{{print $1}}' > {output.read_id}
		python {config[path_to_scripts]}info_to_insert.py -a {input.info_a_1} -A {input.info_A_2} -g {input.info_g_1} -G {input.info_G_2} -o {output.insert_lengths}
		paste -d '\t' {output.read_id} {output.origin_a_1} {output.origin_g_1} {output.origin_A_2} {output.origin_G_2} {output.insert_lengths}  > {output.origin}
		samtools view {input.bam} | awk '{{sub("_"," ")}}{{print $1" "$2}}' > {output.read_id_mapped}
		python {config[path_to_scripts]}group_and_count_reads.py -m {output.read_id_mapped} -r {output.origin} -o {output.read_groups}> {output.off_target_summary}
		"""

rule deduplication_UMIs:
	input:
		bam = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.bam",
		target_cpg_bed = config['target_cpg_bed']
	output:
		header = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.sam_header",
		alignment = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe.sam_alignments",
		sam = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe_barcode.sam",
		deduplicated_sam = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe_barcode.deduplicated.sam",
		cov = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe_barcode.deduplicated.bismark.cov.gz",
		cov_genomic = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe_barcode.deduplicated.bismark_genomic.gz",
		cov_target = config['output_directory'] + "{sample}/" + "{sample}_1_bismark_bt2_pe_barcode.deduplicated.bismark_target_cpgs.cov"
	shell:
		"""
		samtools view -H {input.bam} > {output.header} ;
		samtools view {input.bam} > {output.alignment} ;
		awk -v OFS="\t" '{split($1,a,"_")}{split(a[1],umi,":")}{split($1,r,":")}{split(a[2],read,":")}{$1=""; print r[1]":"r[2]":"r[3]":"r[4]":"r[5]":"r[6]":"r[7]":"r[11]"_"read[1]":"r[9]":"r[10]":"umi[8]$0;}' {output.alignment} > {output.sam} ;
		cat {output.sam} >> {output.header} ;
		mv {output.header} {output.sam} ;
		deduplicate_bismark -p --barcode {output.sam} --output_dir {config[output_directory]}{sample}/;
		bismark_methylation_extractor -p --cutoff 1 --gzip --multicore 3 --comprehensive --bedGraph --zero_based --cytosine_report --no_overlap --genome_folder {config[output_directory]}reference_amplicons/ -o {config[output_directory]}{sample}/ {output.deduplicated_sam} ;
		python {config[path_to_scripts]}targeted_cov_to_genome_cov.py -t {output.cov} -g {output.cov_genomic} -o {config[output_directory]}{sample}/ ;
		bedtools intersect -a {input.target_cpg_bed} -b {output.cov_genomic} -wb > {output.cov_target}
		"""
