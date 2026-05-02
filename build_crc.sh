#!/bin/bash
# Build cm1lp_kj on CRC (epycfe) — Intel + MPI + NetCDF, CPU-only.
set -e

module load intel intelmpi netcdf cuda

SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPTDIR"

# CRC requires mpiifort; patch Makefile if needed
sed -i 's/override FC = mpif90/override FC = mpiifort/' src/Makefile

cd src
make clean
make -j 32 FC=ifort USE_MPI=true USE_NETCDF=true

echo "Build complete: ../run/cm1.exe"
