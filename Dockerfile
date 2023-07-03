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

RUN pip install poetry==1.3.2
# poetry useses virtualenvs by default -> we want global installation
RUN poetry config virtualenvs.create false
ADD pyproject.toml /root/pyproject.toml
ADD poetry.lock /root/poetry.lock
RUN poetry install --no-root

# install percolator
RUN ZIP=ubuntu.tar.gz && \
    wget https://github.com/percolator/percolator/releases/download/rel-3-05/$ZIP -O /tmp/$ZIP && \
    tar xvzf /tmp/$ZIP && \
    chmod -R 755 /tmp/* && \
    dpkg -i percolator-v3-05-linux-amd64.deb && \
    rm /tmp/$ZIP

# Delete ssh keys
ADD oktoberfest/ /root/oktoberfest

# Used by ProteomicsDB runs to describe the oktoberfest version
ADD hash.file /root/hash.file
