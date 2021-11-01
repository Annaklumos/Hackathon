# # TODO:
# remplacer les ressources load par threads et préciser le nombre de core dans la ligne de commande ?
# Fichier config ?



# pour lancer : snakemake --use-singularity --resources load=100
samples = ["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]
type_counts = ["gene", "exon"]
strand = ["0","1","2"]

rule all:
 input:
  expand(["{SAMPLE}_fastqc.html","{SAMPLE}_1_fastqc.html","{SAMPLE}_2_fastqc.html"], SAMPLE=samples)
<<<<<<< HEAD
=======

>>>>>>> cf176a0 (Rule COUNTS analyse tous les strands et position ? + modifie fastqc pour avoir docker (fonctionnait pas sinon ?) + ajout thread (=nb de cpu))

rule fastq: #permet d'obtenir les deux fichiers du paired-end (sample_1.fastq + sample_2.fastq)
 output:
  "{sample}_1.fastq","{sample}_2.fastq"
 resources:
  load=20
 singularity:
   "docker://pegi3s/sratoolkit:2.10.0"
 shell:
  "fasterq-dump {wildcards.sample}"

# Remplacer par ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ? fichier plus
# verifier [!-d]
rule chrm: #télécharge tous les chromosomes humains et regroupe en un seul fasta
 output:
  "ref/human_genome.fa"
 resources:
  load=1
 shell:
  """
   if [ ! -d ref ];then
   	mkdir ref
   fi
   for chr in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "MT" "X" "Y";
   do
   	wget -O "$chr".fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."$chr".fa.gz
   done
   gunzip -c *.fa.gz > {output}
  """
# SA.index
rule index:
 input:
  "ref/human_genome.fa"
 output:
  "ref/SAindex"
 resources: #précedememnt load = 1
  load=1
 threads:
  16
 singularity:
  "docker://evolbioinfo/star:v2.7.6a"
 shell: "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles {input}"

rule mapping_star: #CPU : minimum 16 et ram minimum autant que le génome en théorie
 input:
  index="ref/SAindex",
  genome="ref/human_genome.fa",
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
   --outFileNamePrefix {wildcards.SAMPLE}> {output}
  """


rule mapping_samtools:
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

rule annotation:
 output:
  "human_genome.gtf"
 resources:
  load=1
 shell:
  """
   wget -O {output} ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
  """

rule counting:
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
   if [ ! -d {params.dir} ];then
    mkdir {params.dir}
   fi
   featureCounts -T {threads} -t {wildcards.TYPE} -g gene_id -s {wildcards.STRAND} -a {input.genome} -o {output} {input.bam}
  """

rule fastqc:
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

rule fastqc:
 input:
  bam="{SAMPLE}.bam",
  sample1="{SAMPLE}_1.fastq",
  sample2="{SAMPLE}_2.fastq"
 output:
  "{SAMPLE}_fastqc.html", "{SAMPLE}_1_fastqc.html", "{SAMPLE}_2_fastqc.html"
 resources:
  load=1
 conda:
  "../envs/fastqc.yaml"
 shell:
  """
   fastqc {input.bam}
   fastqc {input.sample1}
   fastqc {input.sample2}
  """

