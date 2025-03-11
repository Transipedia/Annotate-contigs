FROM mambaorg/micromamba:latest

ARG conda_env=AnnotateContig

WORKDIR /pipeline

COPY . .

RUN micromamba create --yes --file annotatecontig.yml --name $conda_env

ENV PATH=/opt/conda/envs/$conda_env/bin:$PATH
ENV MAMBA_DEFAULT_ENV=$conda_env

ENTRYPOINT [ "snakemake", "-s", "/pipeline/Snakefile" ]
