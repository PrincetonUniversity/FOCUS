FROM public.ecr.aws/bitcompat/python:3.11 AS build

SHELL [ "/bin/bash", "-o", "errexit", "-o", "nounset", "-o", "pipefail", "-c" ]

# Setup the build environment with necessary build packages
RUN <<EOT
  sed -i 's/main/main non-free/' /etc/apt/sources.list
  install_packages git gfortran openmpi-bin libopenmpi-dev g++ libhdf5-openmpi-dev hdf5-tools

  /sbin/ldconfig

EOT

# Copy the app directory to the build
COPY . /app
WORKDIR /app

# Set some environment variables
ARG HOSTNAME="docker"
ARG FC="mpif90.openmpi"
ARG CC="gfortran"

# Build the code
WORKDIR /app/sources
RUN make clean
RUN make

# Clean up the directory
#RUN find . -name *.o -exec rm {} \;

#SHELL [ "/bin/bash", "-o", "errexit", "-o", "nounset", "-o", "pipefail", "-c" ]

RUN <<EOT
  install_packages ca-certificates openmpi-bin hdf5-tools libhdf5-openmpi-fortran-102 libhdf5-openmpi-dev ssh rsync zip
  mkdir -p /data
EOT

WORKDIR /data
