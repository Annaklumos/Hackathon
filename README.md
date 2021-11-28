# Projet fait par Pierre BERGERET, Camille RABIER, Lydie TRAN et Dalyan VENTURA

## Installations pré-requises et commande pour lancer le script
La commande à lancer est :
$ snakemake --use-singularity --resources load=100 --cores all

Afin de lancer le workflow, Snakemake et singularity doivent être installés.

## Résumé du workflow
Le script télécharge l'ensemble des FASTQ issus des expériences de RNAseq présent dans la liste "samples", puis télécharge les chromosomes humains et les regroupe en un fichier du génome.Il télécharge aussi le fichier d'annotation du génome humain.

Le génome humain est ensuite indexé à l'aide de STAR, puis les reads issus des FASTQ sont alignés sur le génome à l'aide de STAR.

L'alignement par STAR donne des fichiers bam qui sont ensuite indexés par samtools.

L'ensemble des FASTQ et bam ont été analysés par FASTQC. Les pages résumants ces analyses sont dans le dossier FASTQC

Le workflow utilise ensuite featureCounts pour compter le nombre de reads par gènes et par exons et pour toutes les possbilités de brins (0: unstranded, 1: stranded, 2: reversely stranded).

Finalement, nous avons analysés à l'aide de DESEQ2 et de PCA, les tables de comptages pour les gènes unstranded. Avec nos samples, les stranded et reversely stranded donnant des tables de comptages avec uniquement des 0.

## Ressources demandées par le script
La machine virtuelle utilisée pour les tests correspond à celle de ifb-core avec les caractéristiques suivantes:
- 16 cores
- 64GB de RAM
- 920GB de stockage

Il est cependant normalement possible de se limiter à une machine avec 32GB de RAM au lieu de 64GB. 
