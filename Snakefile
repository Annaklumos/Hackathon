samples = ["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]

rule all:
    input:
        expand("{SAMPLE}_output.count", SAMPLE=samples)

rule fastq: #permet d'obtenir les deux fichiers du paired-end (sample_1.fastq + sample_2.fastq)
    input :
        "{SAMPLE}"
    output:
        "{SAMPLE}_1.fastq", "{SAMPLE}_2.fastq"
    container:
        "docker://pegi3s/sratoolkit"
    shell:
    """
        fasterq-dump {SAMPLE}
    """

rule chrm: #je ne sais pas quoi faire pour le input
    input:
        #???
    output:
        "ref/human_genome.fa"
    shell:
        """
            mkdir ref
            arr_chr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "MT")
            for chr in "${arr_chr[@]}";
            do
                wget -O "$chr".fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."$chr".fa.gz
            done
            gunzip -c *.fa.gz > {output}
        """

rule index:
    input:
        "human_genome.fa"
    output:
        "ref/Log.out"  "ref/chrLength.txt"  "ref/chrName.txt"  "ref/chrNameLength.txt"  "ref/chrStart.txt" "ref/genomeParameters.txt"

    container:
        "docker://evolbioinfo/star:v2.7.6a"
    shell: #peut donner max de CPU de la machine ?
    """
        STAR --runThreadN <nb cpus> --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles {input}
    """

rule mapping_star: #changer nb CPUs
    input:
        "human_genome.fa","{SAMPLE}"
    output:
        "{SAMPLE}.bam"
    container:
        "docker://evolbioinfo/star:v2.7.6a"
    shell:
    """
        STAR --outSAMstrandField intronMotif \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --genomeDir ref \
        --readFilesIn {SAMPLE}_1.fastq {SAMPLE}_2.fastq \
        --runThreadN <Nb CPUS> \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM <Memory in Bytes> \
        > {output}
    """


rule mapping_samtools: # hmm pourquoi *.bam ?
    input:
        "{SAMPLE}.bam"
    output:
        "{SAMPLE}.bai"
    container:
        "docker://evolbioinfo/samtools:v1.11"
    shell:
    """
        samtools index *.bam
    """

rule annotation:
    input:
        "ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz"
    output:
        "human_genome.gtf"
    shell:
        """
            wget -O human_genome.gtf {input}
        """

rule counting:
    input:
        "{SAMPLE}.bam"
    output:
        "{SAMPLE}.count"
    container:
        "docker://evolbioinfo/subread:v2.0.1"
    shell:
        """
            featureCounts -T <CPUS=2> -t gene -g gene_id -s 0 -a human_genome.gtf -o {SAMPLE}.counts {SAMPLE}.bam
        """
