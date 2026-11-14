FROM condaforge/miniforge3:latest

WORKDIR /analysis
RUN conda create -c conda-forge --name cern-hep python=3.11 root uproot awkward vector mplhep pyhf pytest pyyaml -y && \
    /opt/conda/envs/cern-hep/bin/pip install snakemake fsspec-xrootd

COPY . /analysis
ENV PATH="/opt/conda/envs/cern-hep/bin:$PATH"

CMD ["sh", "-c", "PYTHONPATH=. pytest -v && snakemake --cores 1"]
