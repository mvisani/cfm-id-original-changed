FROM ubuntu:20.04
LABEL author="Fei Wang" 
LABEL name="CFM-ID 4 Dev ENV Image"

###########################################################################################
# Set ENV
###########################################################################################
ENV CFM_ROOT /opt/cfm
ENV PATH /opt/cfm/bin:$PATH
ARG BUILD_CFM_TRAIN="ON"
ARG BUILD_CFM_TEST="ON"
ARG MODEL_DIR="/trained_models_cfmid4.0"

## Suppress prompting from apt by setting the DEBIAN_FRONTEND variable,
ARG DEBIAN_FRONTEND=noninteractive

############################################################################################
# DevTools for debug
############################################################################################
RUN apt-get update && apt-get install -y gdb valgrind g++ gcc make wget unzip tar libboost-all-dev cmake build-essential libopenmpi-dev;

############################################################################################
# Build RdKit
############################################################################################
# download rdkit
ARG RDKit_VERSION="2017_09_3"
RUN cd /tmp;\
    wget -O rdkit.tar.gz https://github.com/rdkit/rdkit/archive/Release_${RDKit_VERSION}.tar.gz ;\
    tar xvzf rdkit.tar.gz;\
    rm rdkit.tar.gz;

#WORKDIR /tmp/rdkit-Release_${RDKit_VERSION}
# build and install rdkit
# DRDK_OPTIMIZE_NATIVE set to OFF since we are building for docker
RUN cd /tmp/rdkit-Release_${RDKit_VERSION};\
    mkdir build;\
    cd build;\
    cmake .. -DRDK_PGSQL_STATIC=OFF -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_CPP_TESTS=OFF -DRDK_BUILD_DESCRIPTORS3D=OFF -DRDK_INSTALL_INTREE=OFF -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_OPTIMIZE_NATIVE=OFF -DCMAKE_CXX_STANDARD=14 -DCMAKE_BUILD_TYPE=Release;\
    make install -j 12;\
    mkdir -p /usr/local/include/rdkit/External/INCHI-API/;\
    cp ../External/INCHI-API/*.h /usr/local/include/rdkit/External/INCHI-API/;\
    cd /tmp;\
    rm -rf ./rdkit-Release_${RDKit_VERSION};

RUN echo "/usr/local/lib" > /etc/ld.so.conf.d/local.conf;
RUN ldconfig;

############################################################################################
# Installs an LP Solver 5.5.2.11
############################################################################################
RUN cd /tmp; \
    wget -O lp_solve_dev_ux64.tar.gz 'https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_dev_ux64.tar.gz/download' ;\
    mkdir lp_solve_dev_ux64;\
    tar xvzf lp_solve_dev_ux64.tar.gz -C lp_solve_dev_ux64;\
    rm lp_solve_dev_ux64.tar.gz;\
    mkdir -p /usr/local/include/lp_solve;\
    cp /tmp/lp_solve_dev_ux64/lp_Hash.h /usr/local/include/lp_solve/lp_Hash.h;\
    cp /tmp/lp_solve_dev_ux64/lp_SOS.h /usr/local/include/lp_solve/lp_SOS.h;\
    cp /tmp/lp_solve_dev_ux64/lp_lib.h /usr/local/include/lp_solve/lp_lib.h;\
    cp /tmp/lp_solve_dev_ux64/lp_matrix.h /usr/local/include/lp_solve/lp_matrix.h;\
    cp /tmp/lp_solve_dev_ux64/lp_mipbb.h /usr/local/include/lp_solve/lp_mipbb.h;\
    cp /tmp/lp_solve_dev_ux64/lp_types.h /usr/local/include/lp_solve/lp_types.h;\
    cp /tmp/lp_solve_dev_ux64/lp_utils.h /usr/local/include/lp_solve/lp_utils.h;\
    cp /tmp/lp_solve_dev_ux64/liblpsolve55.so /usr/local/lib/liblpsolve55.so;\
    cp /tmp/lp_solve_dev_ux64/liblpsolve55.a /usr/local/lib/liblpsolve55.a;\
    rm -rf /tmp/lp_solve_dev_ux64;

############################################################################################
# Installs perf
############################################################################################
RUN apt-get install -y flex bison;\
    cd /tmp; \
    wget -O perf-5.10.0.tar.gz https://mirrors.edge.kernel.org/pub/linux/kernel/tools/perf/v5.10.0/perf-5.10.0.tar.gz;\
    tar xvzf perf-5.10.0.tar.gz;\
    cd ./perf-5.10.0/tools/perf;\
    make;\
    cp perf /usr/bin;\
    cd /tmp;\
    rm -rf ./perf-5.10.0;\
    rm -rf ./perf-5.10.0.tar.gz;
