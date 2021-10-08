samples = ['A','B']
rule all:
    input:
        expand("")

rule fastq:
    input: #SRRname
    output: #SRR.fastq
    container:
        "docker://evolbioinf/sratoolkit:v2.5.7"
        # or is there no fasterqdump insisde ?
    shell:
        #fasterq-dump $(SRAID)


rule chrm:
    input: #URL
    output: #.ref.fa
    shell:
        " arr_chr = 1,2,3,4,5,6...
        for i in "${arr_chr[@]}"
        do
            wget -o <chromosome>.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.!{chr}.fa.gz
        done
        gunzip -c *.fa.gz > ref.fa"

rule index:
    input: ref.fa
    output:
    container:
        "docker://evolbioinfo/star:v2.7.6a"
    shell:
        "STAR --runThreadN <nb cpus> --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ref.fa"

rule annotation:
    input:
    output:
    container:
    shell:

rule mapping:
    input:
    output:
    container:
    shell:

rule counting:
    input:
    output:
    container:
    shell:
