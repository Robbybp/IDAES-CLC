#!/bin/sh
# This is not very secure, but someone needs to complain to couenne about their
# certificate.
set -e
cd /usr/local/src
svn co -q https://projects.coin-or.org/svn/Couenne/stable/0.5 couenne --non-interactive --trust-server-cert
cd couenne/ThirdParty/Blas && ./get.Blas
cd ../Lapack && ./get.Lapack
cd ../ASL && ./get.ASL
cd ../Mumps && ./get.Mumps
#cd ../Metis && ./get.Metis
cd ../../ && mkdir build && cd build
../configure -C --prefix=/usr/local
make && make install
