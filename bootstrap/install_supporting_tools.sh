#!/bin/sh

# Installing ADOL-C and pyadolc for wrapper
cd /usr/local/src/Ipopt-3.12.5
wget -q http://www.coin-or.org/download/source/ADOL-C/ADOL-C-2.6.2.tgz
tar -zxf ADOL-C-2.6.2.tgz
cd ADOL-C-2.6.2 && ./configure --prefix=/usr/local && make && make install
ldconfig

cd /opt/conda/envs/python2/lib/python2.7/site-packages
git clone -q https://github.com/b45ch1/pyadolc.git
cd /opt/conda/envs/python2/lib/python2.7/site-packages/pyadolc
sed -i 's/.*raw_input.*/#raw_input/' setup.py
./bootstrap.sh && source activate python2 && python setup.py install
cd /home/jovyan/models

# Installing casadi: 
cd /home/jovyan
wget -q http://files.casadi.org/3.1.1/linux/casadi-py27-np1.9.1-v3.1.1.tar.gz
mkdir casadi-py27-np1.9.1-v3.1.1 && cd casadi-py27-np1.9.1-v3.1.1 && \
	tar -zxvf ../casadi-py27-np1.9.1-v3.1.1.tar.gz

echo "/home/jovyan/casadi-py27-np1.9.1-v3.1.1" >> /opt/conda/envs/python2/lib/python2.7/site-packages/casadi.pth

mkdir -p /opt/conda/envs/python2/etc/conda/activate.d;  
mkdir -p /opt/conda/envs/python2/etc/conda/deactivate.d;
echo '#!/bin/sh' > /opt/conda/envs/python2/etc/conda/activate.d/env_vars.sh;
echo 'export CPPFLAGS="-I/usr/local/src/Ipopt-3.12.5/ThirdParty/ASL -I/opt/conda/pkgs/python-2.7.12-0/include/python2.7"' >> /opt/conda/envs/python2/etc/conda/activate.d/env_vars.sh;
echo 'export MODELS_HOME=/home/jovyan/models' >> /opt/conda/envs/python2/etc/conda/activate.d/env_vars.sh;
echo 'export ADOLC_LIBS=/usr/local/lib64' >> /opt/conda/envs/python2/etc/conda/activate.d/env_vars.sh;
echo '#!/bin/sh' > /opt/conda/envs/python2/etc/conda/deactivate.d/env_vars.sh;
echo 'unset CPPFLAGS' >> /opt/conda/envs/python2/etc/conda/deactivate.d/env_vars.sh;
echo 'unset MODELS_HOME' >> /opt/conda/envs/python2/etc/conda/deactivate.d/env_vars.sh;
echo 'unset ADOLC_LIBS' >> /opt/conda/envs/python2/etc/conda/deactivate.d/env_vars.sh