samples = ['SRA062369','SRA062359']

rule all:
    input:
        expand("")

rule fastq:
    input:
    output: #SRR.fastq
    container:
        "docker://pegi3s/sratoolkit"
        # or is there no fasterqdump insisde ?
    shell:
    """
        fasterq-dump {sample}
    """

rule chrm:
    input: #URL
    output: #.ref.fa
    shell:
        """
         arr_chr = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 Mt
        for chr in "${arr_chr[@]}"
        do
            wget -o "$chr.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.!{chr}.fa.gz"
        done
        gunzip -c *.fa.gz > ref.fa
        """

rule index:
    input: ref.fa
    output:
    container:
        "docker://evolbioinfo/star:v2.7.6a"
    shell:
    """
        STAR --runThreadN <nb cpus> --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ref.fa
    """

rule annotation:
    input:
    output:
    container:
    shell:
    """
    """

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
