# Base OS
FROM python:3.7.9-slim-buster

# metadata
LABEL \
    org.opencontainers.image.title="CBMPy" \
    org.opencontainers.image.version="0.7.25" \
    org.opencontainers.image.description="Platform for constraint based modelling and analysis. CBMPy implements popular analyses such as FBA, FVA, element/charge balancing, network analysis and model editing as well as advanced methods developed specifically for the ecosystem modelling." \
    org.opencontainers.image.url="http://cbmpy.sourceforge.net/" \
    org.opencontainers.image.documentation="http://cbmpy.sourceforge.net/reference/cbmpy.html" \
    org.opencontainers.image.source="https://github.com/biosimulators/Biosimulators_CBMPy" \
    org.opencontainers.image.authors="BioSimulators Team <info@biosimulators.org>" \
    org.opencontainers.image.vendor="BioSimulators Team" \
    org.opencontainers.image.licenses="GPL-3.0-only" \
    \
    base_image="python:3.7.9-slim-buster" \
    version="0.0.1" \
    software="CBMPy" \
    software.version="0.7.25" \
    about.summary="Platform for constraint based modelling and analysis. CBMPy implements popular analyses such as FBA, FVA, element/charge balancing, network analysis and model editing as well as advanced methods developed specifically for the ecosystem modelling." \
    about.home="http://cbmpy.sourceforge.net/" \
    about.documentation="http://cbmpy.sourceforge.net/reference/cbmpy.html" \
    about.license_file="https://github.com/SystemsBioinformatics/cbmpy/blob/master/LICENCE_GPLv3.txt" \
    about.license="SPDX:GPL-3.0-only" \
    about.tags="BioSimulators,mathematical model,constraint-based model,flux balance analysis,simulation,systems biology,computational biology,SBML,SED-ML,COMBINE,OMEX" \
    maintainer="BioSimulators Team <info@biosimulators.org>"

# Install requirements
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        gcc \
        libglpk-dev \
    && pip install glpk numpy python_libsbml scipy sympy \
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
