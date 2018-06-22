FROM continuumio/miniconda

RUN apt-get update && apt-get upgrade -y
#RUN apt-get install -y g++ libz-dev cmake scons libgsl0-dev libncurses5-dev libxml2-dev libxslt1-dev mafft hmmer

WORKDIR /ecgtheow
COPY . /ecgtheow
#RUN ./INSTALL_docker

