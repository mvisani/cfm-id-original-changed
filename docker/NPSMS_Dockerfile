FROM alpine:3.12
LABEL author="Fei Wang" 
LABEL name="NPS-MS Docker Image"

###########################################################################################
# Set ENV
###########################################################################################
ENV CFM_ROOT /opt/cfm
ENV PATH /opt/cfm/bin:$PATH
ARG BUILD_CFM_TRAIN="ON"
ARG BUILD_CFM_TEST="OFF"
ARG MODEL_DIR="/trained_models_npsms"
ARG RDKit_VERSION="2017_09_3"
ARG NPSMS_VERSION="1.0.0"

############################################################################################
# Installs an LP Solver
############################################################################################
RUN	apk update && apk add --no-cache wget tar;\
	cd /tmp;\
	wget -O lp_solve_src.tar.gz 'https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.5/lp_solve_5.5.2.5_source.tar.gz/download' ;\
	tar xvzf lp_solve_src.tar.gz;\
	rm lp_solve_src.tar.gz;\
	apk del wget tar;

WORKDIR /tmp/lp_solve_5.5
#Patch build
# COPY lpsolve55_patched_ccc ./lpsolve55
# Copy headers
RUN	apk update && apk add --no-cache g++ gcc make;\
	mkdir -p /usr/local/include/lp_solve;\
	cp ./lp_Hash.h /usr/local/include/lp_solve/lp_Hash.h;\
	cp ./lp_SOS.h /usr/local/include/lp_solve/lp_SOS.h;\
	cp ./lp_lib.h /usr/local/include/lp_solve/lp_lib.h;\
	cp ./lp_matrix.h /usr/local/include/lp_solve/lp_matrix.h;\
	cp ./lp_mipbb.h /usr/local/include/lp_solve/lp_mipbb.h;\
	cp ./lp_types.h /usr/local/include/lp_solve/lp_types.h;\
	cp ./lp_utils.h /usr/local/include/lp_solve/lp_utils.h;\
	cd ./lpsolve55 && sh ccc;\
	mkdir -p /usr/local/lib/lp_solve;\
	cp ./bin/ux64/liblpsolve55.so /usr/local/lib/liblpsolve55.so;\
	cp ./bin/ux64/liblpsolve55.a /usr/local/lib/liblpsolve55.a;\
 	rm -rf /tmp/lp_solve_5.5;\
	apk del g++ gcc make;

############################################################################################
# Build RdKit
############################################################################################
# download rdkit
RUN apk update && apk add --no-cache wget tar;\
	cd /tmp;\
	wget -O rdkit.tar.gz https://github.com/rdkit/rdkit/archive/Release_${RDKit_VERSION}.tar.gz ;\
	tar xvzf rdkit.tar.gz;\
	rm rdkit.tar.gz;\
	apk del wget tar;

WORKDIR /tmp/rdkit-Release_${RDKit_VERSION}
# build and install rdkit
# DRDK_OPTIMIZE_NATIVE set to OFF since we are building for docker
RUN apk add --no-cache g++ gcc make cmake boost-dev;\
	mkdir build;\
	cd build;\
	cmake .. -DRDK_PGSQL_STATIC=OFF -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_CPP_TESTS=OFF -DRDK_BUILD_DESCRIPTORS3D=OFF -DRDK_INSTALL_STATIC_LIBS=OFF -DRDK_INSTALL_INTREE=OFF -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_OPTIMIZE_NATIVE=OFF -DCMAKE_CXX_STANDARD=11 -DCMAKE_BUILD_TYPE=Release;\
	make install -j 12;\
	mkdir -p /usr/local/include/rdkit/External/INCHI-API/;\
	cp ../External/INCHI-API/*.h /usr/local/include/rdkit/External/INCHI-API/;\
	rm -rf /tmp/rdkit-Release_${RDKit_VERSION};\
	apk del g++ gcc make cmake boost-dev;\
	apk add --no-cache boost-serialization boost-regex;

############################################################################################
# Build MPICH
############################################################################################
# Source is available at http://www.mpich.org/static/downloads/

# Build Options:
# See installation guide of target MPICH version
# Ex: http://www.mpich.org/static/downloads/3.2/mpich-3.2-installguide.pdf
# These options are passed to the steps below
ARG MPICH_VERSION="3.2"
ARG MPICH_CONFIGURE_OPTIONS="--disable-fortran"
ARG MPICH_MAKE_OPTIONS

# Download, build, and install MPICH
RUN mkdir /tmp/mpich-src
WORKDIR /tmp/mpich-src
RUN apk update && apk add --no-cache wget tar g++ gcc make;\
	wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
      && tar xfz mpich-${MPICH_VERSION}.tar.gz  \
      && cd mpich-${MPICH_VERSION}  \
      && ./configure ${MPICH_CONFIGURE_OPTIONS}  \
      && make ${MPICH_MAKE_OPTIONS} && make install -j 6\
      && rm -rf /tmp/mpich-src\
	  && apk del wget tar g++ gcc make;

############################################################################################
# Build NSP-MS/CFM-ID 4 MSML
############################################################################################
# build and install cfm

RUN apk update && apk add --no-cache wget tar;\
	cd /opt;\
	wget -O cfmid.tar.gz https://bitbucket.org/wishartlab/cfm-id-code/get/NPS-MS_${NPSMS_VERSION}.tar.gz ;\
	mkdir cfmid;\
	tar xf cfmid.tar.gz -C cfmid;\
	rm cfmid.tar.gz;\
	mv cfmid/wishartlab* cfmid/v4;\
	mv cfmid/v4/cfm cfm;\
	mkdir -p ${MODEL_DIR}/[M+H]+;\
	mv cfmid/v4/cfm-pretrained-models/npsms/npsms/* ${MODEL_DIR}/[M+H]+;\
	rm -rf cfmid;\
	apk del wget tar;
	
RUN apk add --no-cache g++ gcc make cmake boost-dev;\
	cd /opt/cfm;\
	mkdir build;\
	cd build;\
	cmake .. -DINCLUDE_TESTS=${BUILD_CFM_TEST}\
			 -DINCLUDE_TRAIN=${BUILD_CFM_TRAIN}\
			 -DLPSOLVE_INCLUDE_DIR=/usr/local/include/lp_solve\
			 -DLPSOLVE_LIBRARY_DIR=/usr/local/lib\
			 -DRDKIT_INCLUDE_DIR=/usr/local/include/rdkit\
			 -DRDKIT_INCLUDE_EXT_DIR=/usr/local/include/rdkit/External\
			 -DRDKIT_LIBRARY_DIR=/usr/local/lib\
			 -DCMAKE_CXX_STANDARD=11;\
	make install -j 12;\
	cd .. && rm -rf build;\
	apk del g++ gcc make cmake boost-dev;\
	apk add --no-cache boost-filesystem boost-serialization boost-system boost-thread;\
	if [ ${BUILD_CFM_TRAIN} = "ON" ] ; then apk add --no-cache boost-program_options; fi;\
	if [ ${BUILD_CFM_TEST} = "ON" ] ; then apk add --no-cache boost-unit_test_framework; fi;

############################################################################################
# Create required folders, symlink
############################################################################################
RUN  mkdir /root/output;\
     mkdir -p /cfmid/public/system;\
     mkdir /cfmid/public/spectra;

# Set WORKDIR to /cfmid so that cfm-id-precomputed will be able to use relative paths
# starting with "public".
WORKDIR /cfmid

############################################################################################
# move version label to bottom for faster build
############################################################################################
LABEL npsms.version = NPSMS_VERSION