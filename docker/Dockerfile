FROM continuumio/miniconda3:4.10.3

RUN apt-get update && apt-get install -y procps
RUN conda create -n main -c bioconda -c conda-forge strling==0.5.0
RUN echo "source activate main" > ~/.bashrc
ENV PATH /opt/conda/envs/main/bin:$PATH
