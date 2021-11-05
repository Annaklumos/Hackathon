# # TODO:
# Fichier config ?



# pour lancer : snakemake --use-singularity --resources load=100
samples = ["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]
type_counts = ["gene", "exon"]
strand = ["0","1","2"]

rule all: #recupere fastqc + counts
 input:
  expand(["{SAMPLE}_fastqc.html","{SAMPLE}_1_fastqc.html","{SAMPLE}_2_fastqc.html",  "{TYPE}_strand_{STRAND}/{SAMPLE}_{TYPE}_{STRAND}.counts"], TYPE=type_counts, STRAND=strand, SAMPLE=samples)


rule fastq: #permet d'obtenir les deux fichiers du paired-end (sample_1.fastq + sample_2.fastq)
 output:
  "{sample}_1.fastq","{sample}_2.fastq"
 resources:
  load=20
 singularity:
   "docker://pegi3s/sratoolkit:2.10.0"
 shell:
  "fasterq-dump {wildcards.sample}"


rule chrm: #télécharge tous les chromosomes humains et regroupe en un seul fasta
 output:
  "ref/human_genome.fa"
 resources:
  load=1
 shell:
  """
   for chr in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "MT" "X" "Y";
   do
   	wget -O "$chr".fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."$chr".fa.gz
   done
   gunzip -c *.fa.gz > {output}
  """

rule index: #index genome
 input:
  "ref/human_genome.fa"
 output:
  "ref/SAindex" #dernier fichier généré théoriquement
 resources:
  load=1
 threads:
  16
 singularity:
  "docker://evolbioinfo/star:v2.7.6a"
 shell: "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles {input}"

rule mapping_star: #aligne reads sur genome
 input:
  index="ref/SAindex", #index fait
  genome="ref/human_genome.fa", #chrm fait
  sample1="{SAMPLE}_1.fastq",
  sample2="{SAMPLE}_2.fastq"
 output:
  "{SAMPLE}.bam"
 resources:
  load=1
 singularity:
  "docker://evolbioinfo/star:v2.7.6a"
 threads:
  8
 shell:
  """
   STAR --outSAMstrandField intronMotif \
   --outFilterMismatchNmax 4 \
   --outFilterMultimapNmax 10 \
   --genomeDir ref \
   --readFilesIn {input.sample1} {input.sample2} \
   --runThreadN {threads} \
   --outSAMunmapped None \
   --outSAMtype BAM SortedByCoordinate \
   --outStd BAM_SortedByCoordinate \
   --genomeLoad NoSharedMemory \
   --limitBAMsortRAM 31000000000 \
   --outFileNamePrefix {wildcards.SAMPLE} > {output}
  """


rule mapping_samtools: #indexe bam
 input:
  "{SAMPLE}.bam"
 output:
  "{SAMPLE}.bam.bai"
 resources:
  load=1
 singularity:
  "docker://evolbioinfo/samtools:v1.11"
 shell:
  """
   samtools index {input}
  """

rule annotation: # télécharge gtf du génome
 output:
  "human_genome.gtf"
 resources:
  load=1
 shell:
  """
   wget -O {output}.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
   gunzip {output}
  """

rule counting: #compte nombre de reads par gène|exon et chaque strand
 input:
  bam="{SAMPLE}.bam",
  bai="{SAMPLE}.bam.bai",
  genome="human_genome.gtf"
 output:
  "{TYPE}_strand_{STRAND}/{SAMPLE}_{TYPE}_{STRAND}.counts"
 resources:
  load=1
 singularity:
  "docker://evolbioinfo/subread:v2.0.1"
 params:
  dir="{TYPE}_strand_{STRAND}"
 threads:
  4
 shell:
  """
   featureCounts -T {threads} -t {wildcards.TYPE} -g gene_id -s {wildcards.STRAND} -a {input.genome} -o {output} {input.bam}
  """

rule fastqc: #verifie qualité des fastq et bam
 input:
  bam="{SAMPLE}.bam",
  sample1="{SAMPLE}_1.fastq",
  sample2="{SAMPLE}_2.fastq"
 output:
  "{SAMPLE}_fastqc.html", "{SAMPLE}_1_fastqc.html", "{SAMPLE}_2_fastqc.html"
 resources:
  load=1
 singularity:
  "docker://pegi3s/fastqc"
 shell:
  """
   fastqc {input.bam}
   fastqc {input.sample1}
   fastqc {input.sample2}
  """
