# R + DESeq2 + factomineR + factoextra
# https://cran.r-project.org/
# v4.0.2 / v1.39.0

# base image: r-base:4.0.3
FROM r-base:4.0.3

# File Author / Maintainer
MAINTAINER Camille RABIER <camrabier@gmail.com>

RUN apt-get update && apt-get install -y procps libssl-dev libcurl4-gnutls-dev curl git libopenmpi-dev openmpi-bin openmpi-doc libxml2-dev \
    && R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}' \
    && R -e 'BiocManager::install("DESeq2")' \
    && R -e 'BiocManager::install("FactoMineR")' \
    && R -e 'BiocManager::install("factoextra")'

ENTRYPOINT ["/usr/bin/R"]
