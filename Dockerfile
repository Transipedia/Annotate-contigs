FROM continuumio/miniconda3

ARG conda_env=annot-gen

WORKDIR /annot

COPY . .

RUN conda env create -f environment.yml

ENV PATH /opt/conda/envs/$conda_env/bin:$PATH
ENV CONDA_DEFAULT_ENV $conda_env

RUN /bin/bash -c "source activate annot-gen"

ENTRYPOINT [ "snakemake", "-s", "/annot/Snakefile" ]
