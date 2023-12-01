configfile:     config['path_to_config_file']
rule all:
        input:
                config = config['path_to_config_file'],
                multiqc_1 = config['output_directory'] + "multiqc_report.html",
                target_cpgs_cov = expand(config['output_directory'] + "{sample}_filtered_1_bismark_bt2_pe.bismark_target_cpgs.cov", sample=config["samples"])
#---------------------------------------------------------------------------------------------------
rule copy_config_to_output_directory:
        input:
                config['path_to_config_file']
        output:
                config['output_directory'] + "/config.json"
        shell:
                "cp {input} {output}"
#STEP 1: Prepare reads: A. fastqc & Multiqc, C. cut adapter and probe sequences. ### Add a filter-on-read-quality step! (B)
rule fastqc:
        input:
                config['fastq_folder'] + "{sample}_{read}.fq"
        output:
                html = config['output_directory'] + "{sample}_{read}_fastqc.html",
                zip = config['output_directory'] + "{sample}_{read}_fastqc.zip"
	shell:
               "fastqc {input} -o {config[output_directory]}"

rule multiqc:
        input:
                fastqcs = expand(config['output_directory'] + "{sample}_{read}_fastqc.zip", sample=config["samples"], read = ["1","2"])
        output:
                config['output_directory'] + "multiqc_report.html"
	shell:
               "multiqc {input.fastqcs} -o {config[output_directory]}"

rule create_primer_files_for_filtering:
        input:
                chosen_panel = config['chosen_panel_csv']
        output:
                primers_A = config['output_directory'] + "primers_A.fasta",
		primers_a = config['output_directory'] + "primers_a.fasta"
        shell:
                "python {config[path_to_scripts]}chosen_panel_to_trim_fastas.py -A {output.primers_A} -a {output.primers_a} -c {input.chosen_panel}"
rule perform_cutadapt:
        input:
                fq_in_1 = config['fastq_folder'] + "{sample}_1.fq",
                fq_in_2 = config['fastq_folder'] + "{sample}_2.fq",
                primers_A = config['output_directory'] + "primers_A.fasta",
                primers_a = config['output_directory'] + "primers_a.fasta"
        output:
                fq_out_1 = config['output_directory'] + "{sample}_filtered_1.fq",
                fq_out_2 = config['output_directory'] + "{sample}_filtered_2.fq"
        shell:
                "cutadapt -A file:{input.primers_A} -a file:{input.primers_a} -o {output.fq_out_1} -p {output.fq_out_2} {input.fq_in_1} {input.fq_in_2}"

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
              fq_in_1 = config['output_directory'] + "{sample}_filtered_1.fq",
              fq_in_2 = config['output_directory'] + "{sample}_filtered_2.fq"
	output:
              bam = config['output_directory'] + "{sample}_filtered_1_bismark_bt2_pe.bam"
	shell:
             "bismark --genome {config[output_directory]}reference_amplicons/ --non_directional --multicore 10 -D 15 -R 2 -1 {input.fq_in_1} -2 {input.fq_in_2} -o {config[output_directory]}"
rule bismark_report_and_summary: #summary will combine reports per sample
	input:
             bam = config['output_directory'] + "{sample}_filtered_1_bismark_bt2_pe.bam"
	output:
             report = config['output_directory'] + "{sample}_filtered_1_bismark_bt2_PE_report.html",
             summary = config['output_directory'] + "{sample}_bismark_summary.html"
	shell:
             """
             bismark2report --dir {config[output_directory]}
             bismark2summary {input.bam} -o {config[output_directory]}/bismark_summary
             """
rule bismark_methylation_extractor:
	input:
             bam = config['output_directory'] + "{sample}_filtered_1_bismark_bt2_pe.bam"
	output:
             cov = config['output_directory'] + "{sample}_filtered_1_bismark_bt2_pe.bismark.cov.gz"
	shell:
             "bismark_methylation_extractor -p --cutoff 1 --gzip --multicore 10 --comprehensive --bedGraph --zero_based --cytosine_report --no_overlap --genome_folder {config[output_directory]}reference_amplicons/ -o {config[output_directory]} {input.bam}"

#STEP 3 Probe panel analysis (probe distribution)
rule probe_panel_analysis:
	input:
             targeted_cov = config['output_directory'] + "{sample}_filtered_1_bismark_bt2_pe.bismark.cov.gz"
	output:
             genomic_cov = config['output_directory'] + "{sample}_filtered_1_bismark_bt2_pe.bismark_genomic.cov"
	shell:
             "python {config[path_to_scripts]}targeted_cov_to_genome_cov.py -t {input.targeted_cov} -g {output.genomic_cov} -o {config[output_directory]}"

rule keep_target_cpgs:
	input:
            genomic_cov = config['output_directory'] + "{sample}_filtered_1_bismark_bt2_pe.bismark_genomic.cov",
            target_cpg_bed = config['target_cpg_bed']
	output:
            target_cpgs_cov = config['output_directory'] + "{sample}_filtered_1_bismark_bt2_pe.bismark_target_cpgs.cov"
	shell:
            "bedtools intersect -a {input.genomic_cov} -b {input.target_cpg_bed} > {output.target_cpgs_cov}"
