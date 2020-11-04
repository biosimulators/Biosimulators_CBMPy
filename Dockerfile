# Base OS
FROM continuumio/miniconda3:4.8.2

# metadata
LABEL base_image="continuumio/miniconda3:4.8.2"
LABEL version="0.0.1"
LABEL software="CBMPy"
LABEL software.version="0.7.25"
LABEL about.summary="CBMPy is a platform for constraint based modelling and analysis. CBMPy implements popular analyses such as FBA, FVA, element/charge balancing, network analysis and model editing as well as advanced methods developed specifically for the ecosystem modelling."
LABEL about.home="http://cbmpy.sourceforge.net/"
LABEL about.documentation="http://cbmpy.sourceforge.net/reference/cbmpy.html"
LABEL about.license_file="https://github.com/SystemsBioinformatics/cbmpy/blob/master/LICENCE_GPLv3.txt"
LABEL about.license="GPL-3.0-only"
LABEL about.tags="BioSimulators,mathematical model,constraint-based model,flux balance analysis,simulation,systems biology,computational biology,SBML,SED-ML,COMBINE,OMEX"
LABEL maintainer="BioSimulators Team <info@biosimulators.org>"

# Install requirements
ENV CONDA_ENV=py37 \
    PATH=/opt/conda/envs/${CONDA_ENV}/bin:${PATH}
RUN conda update -y -n base -c defaults conda \
    && conda create -y -n ${CONDA_ENV} python=3.7 \
    && conda install --name ${CONDA_ENV} -y -c bgoli -c SBMLTeam cbmpy \
    && /bin/bash -c "source activate ${CONDA_ENV}" \
    && apt-get update -y \
    && apt-get install -y --no-install-recommends \
        gcc \
        libglpk-dev \
    && pip install glpk \
    && apt-get remove -y \
        gcc \
        libglpk-dev \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_cbmpy
RUN pip install /root/Biosimulators_cbmpy \
    && rm -rf /root/Biosimulators_cbmpy

# Entrypoint
ENTRYPOINT ["cbmpy"]
CMD []
