#! /bin/bash

# This script downloads, configures, builds, and installs FFTW on the totient
# cluster. Simply run `scripts/fftw-install.sh` from the root `ppm` directory.

set -euo pipefail # http://goo.gl/jMkzRN

main() {
    if [[ "$(basename $PWD)" != "ppm" ]]; then
        echo "ERROR: please run this script in the root ppm directory."
        exit -1
    fi

    readonly fftwname="fftw-3.3.4"
    readonly ffturl="http://www.fftw.org/$fftwname.tar.gz"
    readonly installdir="$PWD/fftw"

    # Download and unzip source code
    if [[ ! -f "$fftwname.tar.gz" ]]; then
        wget "$ffturl"
    fi
    if [[ ! -d "$fftwname" ]]; then
        tar -xzvf "$fftwname.tar.gz"
    fi

    # Configure and build code
    mkdir -p "$installdir"
    cd "$fftwname"
    ./configure --prefix="$installdir" CC=icc --enable-openmp --enable-avx
    make
    make install
}

main
