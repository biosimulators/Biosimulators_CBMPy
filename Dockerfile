# Base OS
FROM ubuntu:20.04

# metadata
LABEL base_image="ubuntu:20.04"
LABEL version="1.0.0"
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
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
    && pip3 install -U pip \
    && pip3 install -U setuptools \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Copy code for command-line interface into image and install it
COPY . /root/biosimulators_cbmpy
RUN pip3 install /root/biosimulators_cbmpy

# Entrypoint
ENTRYPOINT ["cbmpy"]
CMD []
