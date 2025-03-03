FROM public.ecr.aws/bitcompat/python:3.11 AS build

SHELL [ "/bin/bash", "-o", "errexit", "-o", "nounset", "-o", "pipefail", "-c" ]

# Setup the build environment with necessary build packages
RUN <<EOT
  sed -i 's/main/main non-free/' /etc/apt/sources.list
  install_packages git gfortran openmpi-bin libopenmpi-dev g++ libhdf5-openmpi-dev hdf5-tools ca-certificates libhdf5-openmpi-fortran-102 ssh rsync zip
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
RUN <<EOT
  mkdir -p /data

  cd sources
  make clean
  make
  find . -name '*.F90' -exec rm {} \;
  find . -name '*.mod' -exec rm {} \;
  find . -name '*.o' -exec rm {} \;
EOT

WORKDIR /data
