#!/bin/bash
# Build cm1lp_kj on CRC (epycfe) — Intel + MPI + NetCDF, CPU-only.
set -e

module load intel intelmpi netcdf cuda

SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPTDIR"

# CRC intel/25+ uses ifx/mpiifx (ifort was removed); patch Makefile regardless of prior state
sed -i 's/override FC  = mp[^ ]*/override FC  = mpiifx/' src/Makefile

cd src
make clean
make -j 32 FC=ifx USE_MPI=true USE_NETCDF=true

echo "Build complete: ../run/cm1.exe"
