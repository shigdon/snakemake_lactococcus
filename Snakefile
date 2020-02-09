shell.executable('/bin/bash')

localrules: all, copy_assembly, gather_gff, compute_sig, copy_sigs, compare_sigs

import glob, os

outdir = config['outdir'] + 'data/'
logdir = config['outdir'] + 'log/'
analysis_dir = config['outdir'] + 'analysis/'
gbdir = config['gb_asm']

#--------Lactococcus Isolate Dictionary--------#

r1_dict = {}
r1_list = glob.glob('reads/*_R1.fastq.gz')

for name in r1_list:
    i = os.path.basename(name).split('_')[0]
    r1_dict[i] = (name, name.replace('_R1', '_R2'))

def lookup_isolate(wildcards):
   return r1_dict[wildcards.isolate]

print(r1_dict)

#--------GenBank Assemblies--------#

GB_IDS, = glob_wildcards("genbank_asm/{id}.fna.gz")

#--------Workflow--------#

rule all:
    input:
        expand(outdir + "{isolate}/sourmash/{isolate}.sig", isolate=r1_dict.keys()),
        expand(outdir + "{isolate}/prokka", isolate=r1_dict.keys()),
        expand(outdir + "{id}/prokka", id=GB_IDS),
        expand(outdir + "{isolate}/quast", isolate=r1_dict.keys()),
        analysis_dir + "roary/gff-inputs",
        analysis_dir + "roary/roary_out",
        analysis_dir + "sourmash/smash_k31-cmp.csv",
        analysis_dir + "sourmash/smash.k31.lca-classify.csv"

rule trim_isolate:
    input:
        lookup_isolate

    output:
        outdir + "{isolate}/{isolate}.trim_R1.fq.gz",
        outdir + "{isolate}/{isolate}.orphan_R1.fq.gz",
        outdir + "{isolate}/{isolate}.trim_R2.fq.gz",
        outdir + "{isolate}/{isolate}.orphan_R2.fq.gz"

    log:
        logdir + "trimmomatic/{isolate}.log"

    conda:
        "trimmomatic.yml"

    threads: 8

    shell:
        '''
        trimmomatic PE \
            -threads {threads} \
            {input[0]} \
            {input[1]} \
            {output} \
            ILLUMINACLIP:TruSeq3-PE-2.fa:2:40:15:8:TRUE \
            LEADING:2 \
            TRAILING:2 \
            SLIDINGWINDOW:4:15 \
            MINLEN:50 \
            > /dev/null 2> {log}
        '''

rule assemble_isolate:
    input:
        outdir + "{isolate}/{isolate}.trim_R1.fq.gz",
        outdir + "{isolate}/{isolate}.trim_R2.fq.gz"

    output:
        outdir + "{isolate}/megahit"

    conda:
        "megahit.yml"

    threads: 8

    shell:
        '''
        megahit \
            -1 {input[0]} \
            -2 {input[1]} \
            --out-dir {output} \
            -t {threads}
        '''

rule copy_assembly:
    input:
        outdir + "{isolate}/megahit"

    output:
        outdir + "{isolate}/{isolate}.contigs.fa"

    shell:
        '''
        cp {input}/final.contigs.fa {output}
        '''

rule quast:
    input:
        outdir + "{isolate}/{isolate}.contigs.fa"

    output:
        outdir + "{isolate}/quast"

    conda:
        "quast.yml"

    threads: 4

    shell:
        '''
        quast.py \
            {input} \
            --output-dir {output} \
            --labels {wildcards.isolate}
        '''

rule prokka_annotation:
    input:
        outdir + "{isolate}/{isolate}.contigs.fa"

    output:
        outdir + "{isolate}/prokka"

    conda:
        "prokka.yml"

    threads: 8

    shell:
        '''
        prokka \
            --outdir {output} \
            --prefix {wildcards.isolate} \
            --locustag {wildcards.isolate} \
            --cpus {threads} \
            {input}
        '''

rule prokka_gb_asm:
    input:
        expand(gbdir + "{id}.fna.gz", id=GB_IDS)

    output:
        outdir + "{id}/prokka"

    conda:
        "prokka.yml"

    threads: 8

    shell:
        '''
        prokka \
            --outdir {output} \
            --prefix {id} \
            --locustag {id} \
            --cpus {threads} \
            {input}
        '''

rule gather_gff:
    input:
        expand(outdir + "{isolate}/prokka/{isolate}.gff", isolate=r1_dict.keys()),
        expand(outdir + "{id}/prokka/{id}.gff", id=GB_IDS)

    output:
        expand(analysis_dir + "roary/gff-inputs/{isolate}.gff", isolate=r1_dict.keys()),
        expand(analysis_dir + "roary/gff-inputs/{id}.gff", id=GB_IDS)

    shell:
        '''
        cp {input[0]} {output}[0] && \
        cp {input[1]} {output}[1]
        '''

ROARY_INPUTS = glob_wildcards(analysis_dir + "roary/{id}.gff")

rule roary:
    input:
        expand(analysis_dir + "roary/{id}.gff", id=ROARY_INPUTS)

    output:
        analysis_dir + "roary/roary_out"

    conda:
        "roary.yml"

    threads: 24

    shell:
        '''
        roary \
            -e \
            -p {threads} \
            -i 95 \
            -f {output} \
            -o clustered_proteins \
            {input}
        '''

rule compute_sig:
    input:
        outdir + "{isolate}/{isolate}.contigs.fa"

    output:
        outdir + "{isolate}/sourmash/{isolate}.sig"

    conda:
        "sourmash.yml"

    shell:
        '''
        sourmash compute \
            -k 31 \
            --scaled 2000 \
            {input} \
            -o {output}
        '''

rule copy_sigs:
    input:
        outdir + "{isolate}/sourmash/{isolate}.sig"

    output:
        analysis_dir + "sourmash/{isolate}.sig"

    shell:
        '''
        cp {input} {output}
        '''

rule compare_sigs:
    input:
        expand(analysis_dir + "sourmash/{isolate}.sig", isolate=r1_dict.keys())
    output:
        analysis_dir + "sourmash/smash.k31.cmp",
        analysis_dir + "sourmash/smash.k31.cmp.labels.txt",
        analysis_dir + "sourmash/smash_k31-cmp.csv"
    conda:
        "sourmash.yml"
    shell:
        '''
        sourmash compare \
            -k 31 {input} \
            -o {output[0]} \
            --csv {output[2]}
        '''

rule smash_lca_classify:
    input:
        sigs = expand(analysis_dir + "sourmash/{isolate}.sig", isolate=r1_dict.keys()),
        db = config['smash_lca_db']

    output:
        analysis_dir + "sourmash/smash.k31.lca-classify.csv",

    conda:
        "sourmash.yml"

    shell:
        '''
        sourmash lca classify \
            --db {input.db} \
            --query {input.sigs} \
            -o {output} \
        '''

