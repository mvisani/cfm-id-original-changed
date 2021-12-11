FROM alpine:3.11
LABEL author="Fei Wang" 
LABEL name="Latest CFM-ID 4 Docker MSML Image"
LABEL cfmid.msml.version="latest"

###########################################################################################
# Set ENV
###########################################################################################
ENV CFM_ROOT /opt/cfm
ENV PATH /opt/cfm/bin:$PATH
ARG BUILD_CFM_TRAIN="ON"
ARG BUILD_CFM_TEST="OFF"
ARG MODEL_DIR="/trained_models_cfmid4.0"
ARG MAKE_PROCESSOR_COUNT=4

############################################################################################
# Installs Boost and other apks
############################################################################################
RUN apk update;\ 
	apk add --no-cache wget tar g++ gcc make cmake boost-dev;

############################################################################################
# Installs an LP Solver
############################################################################################
RUN	cd /tmp;\
	wget -O lp_solve_src.tar.gz https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.5/lp_solve_5.5.2.5_source.tar.gz/download;\
	tar xvzf lp_solve_src.tar.gz;\
	rm lp_solve_src.tar.gz;

WORKDIR /tmp/lp_solve_5.5
#Patch build
COPY lpsolve55_patched_ccc ./lpsolve55
# Copy headers
RUN	mkdir -p /usr/local/include/lp_solve;\
	cp ./lp_Hash.h /usr/local/include/lp_solve/lp_Hash.h;\
	cp ./lp_SOS.h /usr/local/include/lp_solve/lp_SOS.h;\
	cp ./lp_lib.h /usr/local/include/lp_solve/lp_lib.h;\
	cp ./lp_matrix.h /usr/local/include/lp_solve/lp_matrix.h;\
	cp ./lp_mipbb.h /usr/local/include/lp_solve/lp_mipbb.h;\
	cp ./lp_types.h /usr/local/include/lp_solve/lp_types.h;\
	cp ./lp_utils.h /usr/local/include/lp_solve/lp_utils.h;\
	cd ./lpsolve55 && sh lpsolve55_patched_ccc;\
	mkdir -p /usr/local/lib/lp_solve;\
	cp ./bin/ux64/liblpsolve55.so /usr/local/lib/liblpsolve55.so;\
	cp ./bin/ux64/liblpsolve55.a /usr/local/lib/liblpsolve55.a;
RUN rm -rf /tmp/lp_solve_5.5

############################################################################################
# Build RDkit
############################################################################################
# download rdkit

RUN cd /tmp;\
	wget -O rdkit.tar.gz https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz;\
	tar xvzf rdkit.tar.gz;\
	rm rdkit.tar.gz

WORKDIR /tmp/rdkit-Release_2017_09_3
# build and install rdkit
RUN mkdir build;\
	cd build;\
	# DRDK_OPTIMIZE_NATIVE set to OFF since we are building for docker
	cmake .. -DRDK_PGSQL_STATIC=OFF\
		 	  -DRDK_BUILD_PYTHON_WRAPPERS=OFF\
			  -DRDK_BUILD_CPP_TESTS=OFF  \
			  -DRDK_BUILD_DESCRIPTORS3D=OFF\
			  -DRDK_INSTALL_STATIC_LIBS=OFF\
			  -DRDK_INSTALL_INTREE=OFF\
			  -DRDK_BUILD_INCHI_SUPPORT=ON\ 
			  -DCMAKE_CXX_STANDARD=11\
			  -DCMAKE_BUILD_TYPE=Release;\
	make install -j ${MAKE_PROCESSOR_COUNT};
# Install External and Clean Up
RUN mkdir -p /usr/local/include/rdkit/External/INCHI-API/;\
	cp ./External/INCHI-API/*.h /usr/local/include/rdkit/External/INCHI-API/;\
	rm -rf /tmp/rdkit-Release_2017_09_3;

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
RUN wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
      && tar xfz mpich-${MPICH_VERSION}.tar.gz  \
      && cd mpich-${MPICH_VERSION}  \
      && ./configure ${MPICH_CONFIGURE_OPTIONS}  \
      && make ${MPICH_MAKE_OPTIONS} && make install -j ${MAKE_PROCESSOR_COUNT}\
      && rm -rf /tmp/mpich-src

############################################################################################
# Build CFM-ID 4 MSML
############################################################################################
# build and install cfm
RUN cd /opt;\
	wget -O cfmid.tar.gz https://bitbucket.org/wishartlab/cfm-id-code/get/CFM-ID_4.0.5.tar.gz;\
	mkdir cfmid;\
	ls /opt;\
	tar xf cfmid.tar.gz -C cfmid;\
	rm cfmid.tar.gz;\
	mv cfmid/wishartlab* cfmid/v4;\
	mv cfmid/v4/cfm cfm;\
	mv cfmid/v4/cfm-pretrained-models ${MODEL_DIR};\
	rm -rf cfmid;\
	ls /opt/cfm;

WORKDIR /opt/cfm
RUN	mkdir build;\
	cd build;\
	cmake .. -DINCLUDE_TESTS=${BUILD_CFM_TEST}\
			 -DINCLUDE_TRAIN=${BUILD_CFM_TRAIN}\
			 -DLPSOLVE_INCLUDE_DIR=/usr/local/include/lp_solve\
			 -DLPSOLVE_LIBRARY_DIR=/usr/local/lib\
			 -DRDKIT_INCLUDE_DIR=/usr/local/include/rdkit\
			 -DRDKIT_INCLUDE_EXT_DIR=/usr/local/include/rdkit/External\
			 -DRDKIT_LIBRARY_DIR=/usr/local/lib\
			 -DCMAKE_CXX_STANDARD=11;\
	make install -j ${MAKE_PROCESSOR_COUNT} && make clean;\
	cd ..;\
	rm -rf build;

########################################################################
# clean Up
########################################################################
RUN apk del wget tar g++ gcc make cmake boost-dev;
RUN apk add --no-cache boost-filesystem boost-regex boost-serialization boost-system boost-thread;

RUN if [ ${BUILD_CFM_TRAIN} = "ON" ] ; then apk add --no-cache boost-program_options; fi;
RUN if [ ${BUILD_CFM_TEST} = "ON" ] ; then apk add --no-cache boost-unit_test_framework; fi;

RUN rm -rf /tmp/* /var/tmp/*
############################################################################################
# output folder
############################################################################################
RUN  mkdir /root/output;
RUN ls /tmp