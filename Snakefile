# pour lancer : snakemake --use-singularity --resources load=100
samples = ["SRR628582","SRR628583","SRR628584","SRR628585","SRR628586","SRR628587","SRR628588","SRR628589"]

rule all:
 input:
  expand("{SAMPLE}.counts", SAMPLE=samples)

rule fastq: #permet d'obtenir les deux fichiers du paired-end (sample_1.fastq + sample_2.fastq)
 output:
  "{sample}_1.fastq","{sample}_2.fastq"
 resources:
  load=25
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
   if [! -d ref];then
   	mkdir ref
   fi
   for chr in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "MT";
   do
   	wget -O "$chr".fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."$chr".fa.gz
   done
   gunzip -c *.fa.gz > {output}
  """

rule index:
 input:
  "ref/human_genome.fa"
 output:
  "ref/Log.out","ref/chrLength.txt","ref/chrName.txt","ref/chrNameLength.txt","ref/chrStart.txt"
 resources:
  load=1
 singularity:
  "docker://evolbioinfo/star:v2.7.6a"
 shell: "STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles {input}"

rule mapping_star:
 input:
  index="ref/Log.out",
  genome="ref/human_genome.fa",
  sample1="{SAMPLE}_1.fastq",
  sample2="{SAMPLE}_2.fastq"
 output:
  "{SAMPLE}.bam"
 resources:
  load=100
 singularity:
  "docker://evolbioinfo/star:v2.7.6a"
 shell:
  """
   STAR --outSAMstrandField intronMotif \
   --outFilterMismatchNmax 4 \
   --outFilterMultimapNmax 10 \
   --genomeDir ref \
   --readFilesIn {input.sample1} {input.sample2} \
   --runThreadN 20 \
   --outSAMunmapped None \
   --outSAMtype BAM SortedByCoordinate \
   --outStd BAM_SortedByCoordinate \
   --genomeLoad NoSharedMemory \
   --limitBAMsortRAM 42000000000 \
   > {output}
  """


rule mapping_samtools:
 input:
  "{SAMPLE}.bam"
 output:
  "{SAMPLE}.bai"
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
  genome="human_genome.gtf"
 output:
  "{SAMPLE}.counts"
 resources:
  load=1
 singularity:
  "docker://evolbioinfo/subread:v2.0.1"
 shell:
  """
   featureCounts -T 16 -t gene -g gene_id -s 0 -a {input.genome} -o {output} {input.bam}
  """
