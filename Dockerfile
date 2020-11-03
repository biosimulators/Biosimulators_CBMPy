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
RUN conda create -y -n py37 python=3.7 \
    && conda activate py37
    && conda install -y -c bgoli -c sbmlteam cbmpy

# Copy code for command-line interface into image and install it
COPY . /root/biosimulators_cbmpy
RUN pip install /root/biosimulators_cbmpy \
    && rm -rf /root/biosimulators_cbmpy

# Entrypoint
ENTRYPOINT ["cbmpy"]
CMD []
