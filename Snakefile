# java >= 11

configfile: "config.yaml"

samples, = glob_wildcards("input/{id}.fastq")

rule all:
    input:
        expand("preqc/{id}_fastqc.html", id=samples),
        #expand("postqc/{id}_nanostat.txt", id=samples),
        expand("postqc/{id}_fastqc.html", id=samples),
        expand("rel-abundances/{id}_qc_rel-abundance.tsv", id=samples),
        "output/formated_emu_table_combined_counts.tsv",
        "output/krona_composition_samples.html"
    threads: 
        workflow.cores

rule fastqc_pre:
    input: 
        "input/{id}.fastq"
    output: 
        "preqc/{id}_fastqc.html",
        temp("preqc/{id}_fastqc.zip")
    log:
        "logs/preqc/{id}_fastqc.log"
    threads: 
        workflow.cores
    shell:
        "scripts/fast_qc/fastqc -o preqc -t {workflow.cores} {input} &> {log}"

rule porechop:
    input: 
        "input/{id}.fastq"
    output:
        temp("postqc/{id}_porechop.fastq")
    log:
        "logs/postqc/{id}_porechop.log"
    threads: 
        workflow.cores
    shell:
        "porechop -t {workflow.cores} -i {input} -o {output} &> {log}"

rule nanofilt:
    input:
        "postqc/{id}_porechop.fastq"
    output:
        temp("postqc/{id}_porechop_nanofilt.fastq")
    params:
        min_len = config["nanofilt"]["minlen"],
        max_len = config["nanofilt"]["maxlen"],
        quality = config["nanofilt"]["quality"]
    log:
        "logs/postqc/{id}_nanofilt.log"
    threads: 
        workflow.cores
    shell:
        "NanoFilt --logfile {log} -l {params.min_len} --maxlength {params.max_len} --q {params.quality} {input} > {output} "

rule minimap:
    input:
        "postqc/{id}_porechop_nanofilt.fastq"
    output:
        temp("postqc/{id}_porechop_nanofilt.paf")
    log:
        "logs/postqc/{id}_minimap.log"
    threads: 
        workflow.cores
    shell: 
        "minimap2 -x ava-ont -g 500 -t {workflow.cores} {input} {input} > {output} 2> {log}"

rule yacrd:
    input:
        "postqc/{id}_porechop_nanofilt.paf",
        "postqc/{id}_porechop_nanofilt.fastq"
    output:
        temp("postqc/{id}.yacrd"),
        "postqc/{id}.scrubb.fastq"
    log:
        "logs/postqc/{id}_yacrd"
    threads: 
        workflow.cores
    shell:
        "yacrd -t {workflow.cores} -i {input[0]} -o {output[0]} -c 4 -n 0.4 scrubb -i {input[1]} -o {output[1]} &> {log}"

rule fastqc_post:
    input: 
        "postqc/{id}.scrubb.fastq"
    output: 
        "postqc/{id}.scrubb_fastqc.html",
        temp("postqc/{id}.scrubb_fastqc.zip")
    log:
        "logs/postqc/{id}.scrubb_fastqc.log"
    threads: 
        workflow.cores
    shell:
        "scripts/fast_qc/fastqc -o postqc -t {workflow.cores} {input} &> {log}"

rule rename:
    input:
        "postqc/{id}.scrubb.fastq",
        "postqc/{id}.scrubb_fastqc.html"
    output:
        "postqc/{id}_qc.fastq",
        "postqc/{id}_fastqc.html"
    shell:
        """
        mv {input[0]} {output[0]}
        mv {input[1]} {output[1]}
        """

""" rule nanostat:
    input:
        "postqc/{id}_qc.fastq"
    output:
        "postqc/{id}_nanostat.txt"
    log:
        "logs/postqc/{id}_nanostat"
    threads:  
        workflow.cores
    shell:
        "NanoStat -t {workflow.cores} --fastq {input} > {output} 2> {log}"
 """
rule emu:
    input:
        "postqc/{id}_qc.fastq"
    output:
        "rel-abundances/{id}_qc_rel-abundance.tsv"
    log:
        "logs/emu/{id}_emu"
    params:
        db = config["emu"]["db"]
    threads: 
        workflow.cores
    shell:
        "emu abundance --keep-counts --type map-ont --db databases/{params.db} --output-dir rel-abundances --threads {workflow.cores} {input} &> {log}"

rule combine_outputs:
    input:
        "rel-abundances"
    output:
        "output/emu-combined-tax_id-counts.tsv"
    log:
        "logs/combine_outputs.log"
    shell:
        """
        emu combine-outputs rel-abundances tax_id --counts > {log}
        mv rel-abundances/emu-combined-tax_id-counts.tsv output/emu-combined-tax_id-counts.tsv
        """

rule counts_krona:
    input:
        "output/emu-combined-tax_id-counts.tsv"
    output:
        "output/formated_emu_table_combined_counts.tsv",
        "output/krona_composition_samples.html"
    log:
        "logs/count_krona.log"
    shell:  
        "python scripts/emu_table_2_krona_v2.py -t {input} -o output"

rule barplot_and_heatmap:
    input:
        "output/formated_emu_table_combined_counts.tsv"