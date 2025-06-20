FROM python:3.11-slim

# Avoid interactive prompts during package install
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    zip \
    ssh \
    mono-devel \
 && rm -rf /var/lib/apt/lists/*

# set root directory
ENV HOME=/root
WORKDIR /root

RUN pip install poetry==2.1.3
# poetry useses virtualenvs by default -> we want global installation
RUN poetry config virtualenvs.create false
ADD pyproject.toml /root/pyproject.toml
ADD poetry.lock /root/poetry.lock
RUN poetry install --no-root

# install percolator
RUN DEB=percolator-v3-06-linux-amd64.deb && \
    wget https://github.com/percolator/percolator/releases/download/rel-3-06-01/$DEB -O /tmp/$DEB && \
    dpkg -i /tmp/$DEB && \
    rm /tmp/$DEB

# install ThermoRawFileParser
RUN ZIP=ThermoRawFileParser1.4.3.zip && \
    wget https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.3/$ZIP -O /tmp/$ZIP && \
    unzip /tmp/$ZIP -d /opt/compomics/ && \
    rm /tmp/$ZIP

# Copy source folder
ADD oktoberfest/ /root/oktoberfest

# Copy make files into root folder to allow bootstrapping
ADD Makefile* /root/

# Used by ProteomicsDB runs to describe the oktoberfest version
ADD hash.file /root/hash.file
