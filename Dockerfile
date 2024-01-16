FROM python:3.8.12

# Tell docker that we don't want to be bothered with questions
ARG DEBIAN_FRONTEND=noninteractive

# for mono installation
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb https://download.mono-project.com/repo/debian stable-buster main" | tee /etc/apt/sources.list.d/mono-official-stable.list

RUN apt-get update && apt-get install -y \
        mono-devel \
        ssh \
        zip \
    && rm -rf /var/lib/apt/lists/*

# set root directory
ENV HOME /root
WORKDIR /root

RUN pip install poetry==1.6.1
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
RUN ZIP=ThermoRawFileParser1.4.2.zip && \
    wget https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.2/$ZIP -O /tmp/$ZIP && \
    unzip /tmp/$ZIP -d /root/ && \
    rm /tmp/$ZIP

# Copy source folder
ADD oktoberfest/ /root/oktoberfest

# Copy make files into root folder to allow bootstrapping
ADD Makefile* /root/

# Used by ProteomicsDB runs to describe the oktoberfest version
ADD hash.file /root/hash.file
