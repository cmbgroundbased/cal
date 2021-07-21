FROM debian:buster

LABEL Description="Pycal Docker images"

ENV DEBIAN_FRONTEND noninteractive

RUN bash -c ' \
	set -euxo pipefail && \
	apt-get update && \
	apt-get upgrade -y && \
	apt-get install -y --no-install-recommends \
		vim \ 
		build-essential \ 
		wget \
		gfortran \ 
		bison \ 
		libibverbs-dev \ 
		libibmad-dev \ 
		libibumad-dev \ 
		librdmacm-dev \ 
		apt-utils \ 
		gcc \ 
		g++ \ 
		kmod \ 
		openssh-client \ 
		python3 \ 
		python3-pip \ 
		python3-mpi4py \ 
		python3-dev \ 
		git \ 
		autoconf \ 
		locales \ 
		libtool \ 
		libcfitsio-dev \
		libfftw3-3 \
		python3-ephem \
		cmake \
		gfortran \
		libopenblas-dev \
		liblapack-dev \
		libfftw3-dev \
		libsuitesparse-dev \
		libpython3-dev \
		python3-scipy \
		python3-matplotlib \
		python3-healpy \
		python3-astropy \
		python3-yaml \
		curl \
		libhwloc-dev \
		pkg-config && \
	wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.0.tar.gz && \
	tar xfz openmpi-4.0.0.tar.gz && \
	cd openmpi-4.0.0 && \
	./configure --with-slurm --with-pmix --enable-static --enable-openib-udcm --enable-openib-rdmacm && \
	make -j12 && \
	make install && \
	ldconfig && \
	mpicc examples/hello_c.c -o hello_c.x'

# LibAATM

RUN bash -c ' \	
	cd /usr/local/src/ && \
	git clone https://github.com/hpc4cmb/libaatm && \
	cd libaatm && \
	mkdir build && \
	cd build && \
	cmake ../ && \
	make && make install'

CMD ["/bin/bash"]

