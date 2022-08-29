FROM python:3.8.12

# Tell docker that we don't want to be bothered with questions
ARG DEBIAN_FRONTEND=noninteractive

# for mono installation
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb https://download.mono-project.com/repo/debian stable-buster main" | tee /etc/apt/sources.list.d/mono-official-stable.list

