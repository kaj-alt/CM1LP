#!/bin/bash
# Build cm1lp_kj on Derecho — NVHPC + OpenACC (GPU) + MPI + NetCDF.
set -e

module load ncarenv/23.06 nvhpc/23.1 cuda/11.7.1 craype/2.7.20 cray-mpich/8.1.25 \
            ncarcompilers/1.0.0 hdf5/1.12.2 netcdf/4.9.2 cdo/2.1.1 cutensor/1.7.0.1

SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPTDIR"

# Derecho requires mpif90; patch Makefile back if a CRC build changed it
sed -i 's/override FC  = mpiifort/override FC  = mpif90/' src/Makefile

cd src
make clean
make -j 32 USE_MPI=true USE_OPENACC=true GPU_TYPE=a100 USE_NETCDF=true

echo "Build complete: ../run/cm1.exe"
