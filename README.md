Projet fait par Pierre BERGERET, Camille RABIER, Lydie TRAN et Dalyan VENTURA
La machine virtuelle utilisée pour les tests correspond à celle de ifb-core avec les caractéristiques suivantes:
- 16 cores
- 64GB de RAM
- 920GB de stockage
La commande à lancer est :
$ snakemake --use-singularity --resources load=100 --cores all
