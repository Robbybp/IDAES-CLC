#!/bin/sh
set -e
cd /usr/local/src
wget -q http://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.5.tgz
gunzip Ipopt-3.12.5.tgz 
tar xf Ipopt-3.12.5.tar
cd Ipopt-3.12.5/ThirdParty/Blas && ./get.Blas
cd ../Lapack && echo `pwd` && ./get.Lapack
cd ../ASL
sed 's/coinasl=1.3.0/coinasl=3.1.0/g' get.ASL > get.ASL.idaes
chmod +x get.ASL.idaes && ./get.ASL.idaes
cd ../Mumps && ./get.Mumps
#cd ../Metis && ./get.Metis
cd ../../ && ./configure --prefix=/usr/local && make && make install
