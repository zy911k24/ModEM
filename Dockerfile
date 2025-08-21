FROM ubuntu:latest

ARG ncpus=1

WORKDIR /root/

RUN apt-get update
RUN echo "Building with ${ncpus} processes"

RUN apt-get -y install gcc-14
RUN apt-get -y install g++-14
RUN apt-get -y install gfortran-14
RUN apt-get -y install libblas-dev
RUN apt-get -y install liblapack-dev
RUN apt-get -y install make
RUN apt-get -y install wget
RUN apt-get -y install python3
RUN apt-get -y install vim
RUN apt-get -y install git

RUN wget https://www.mpich.org/static/downloads/4.3.1/mpich-4.3.1.tar.gz
RUN tar -xzvf mpich-4.3.1.tar.gz

# Install MPICH Manually as the apt-get install uses pmix and 
# that will mess with things
WORKDIR mpich-4.3.1
RUN CC=gcc-14 CXX=g++-14 FC=gfortran-14 ./configure --with-pm=hydra
RUN make -j${ncpus} && make -j${ncpus} install
ENV LD_LIBRARY_PATH=/usr/local/lib/
WORKDIR /root/
RUN rm -rf mpich-4.3.1.tar.gz mpich-4.3.1

