FROM ubuntu:latest

RUN apt-get update
RUN apt-get -y install gcc
RUN apt-get -y install g++
RUN apt-get -y install gfortran
RUN apt-get -y install gdb
RUN apt-get -y install valgrind
RUN apt-get -y install libblas-dev
RUN apt-get -y install liblapack-dev
RUN apt-get -y install make
RUN apt-get -y install wget
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
RUN apt-get -y install pipx
RUN apt-get -y install vim
RUN apt-get -y install git
RUN apt-get -y install openmpi-bin openmpi-common libopenmpi-dev
RUN useradd -m -r -g users -s /bin/bash modem

ENV LD_LIBRARY_PATH=/usr/local/lib/
# For Python Installs
ENV PATH="$PATH:/root/.local/bin"

COPY . /home/modem/ModEM
RUN chown -R modem /home/modem/ModEM
RUN chgrp -R users /home/modem/ModEM

USER modem
WORKDIR /home/modem/ModEM/f90

# Make SP2 MPI
RUN ./CONFIG/configure -m mpi -l MF -g release Makefile gfortran && make && mv Mod3DMT Mod3DMT_MF && make clean
RUN ./CONFIG/configure -m mpi -l SP -g release Makefile gfortran && make && mv Mod3DMT Mod3DMT_SP && make clean
RUN ./CONFIG/configure -m mpi -l SP2 -g release Makefile gfortran && make && mv Mod3DMT Mod3DMT_SP2 && make clean

WORKDIR /home/modem/