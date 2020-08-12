FROM alpine:3.11
LABEL author="Fei Wang" 
LABEL name="CFM-ID 4 Docker Image"
LABEL image_version = "4.1.1.6"
LABEL cfm_msml_version="4.0.5"
LABEL msrb_version="1.1.0"

###########################################################################################
# Set ENV
###########################################################################################
ENV CFM_ROOT /opt/cfm
ENV PATH /opt/cfm/bin:$PATH

############################################################################################
# Installs Boost and other apks
############################################################################################
RUN apk update;\ 
	apk add --no-cache openmpi-dev --repository=http://dl-cdn.alpinelinux.org/alpine/edge/testing;\
	apk add --no-cache wget tar g++ gcc make cmake boost-dev;

############################################################################################
# Installs an LP Solver
############################################################################################
RUN	cd /tmp;\
	wget -O lp_solve_src.tar.gz https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.5/lp_solve_5.5.2.5_source.tar.gz/download;\
	tar xvzf lp_solve_src.tar.gz;\
	rm lp_solve_src.tar.gz;

# Copy headers
RUN	mkdir -p /usr/local/include/lp_solve;\
	cp /tmp/lp_solve_5.5/lp_Hash.h /usr/local/include/lp_solve/lp_Hash.h;\
	cp /tmp/lp_solve_5.5/lp_SOS.h /usr/local/include/lp_solve/lp_SOS.h;\
	cp /tmp/lp_solve_5.5/lp_lib.h /usr/local/include/lp_solve/lp_lib.h;\
	cp /tmp/lp_solve_5.5/lp_matrix.h /usr/local/include/lp_solve/lp_matrix.h;\
	cp /tmp/lp_solve_5.5/lp_mipbb.h /usr/local/include/lp_solve/lp_mipbb.h;\
	cp /tmp/lp_solve_5.5/lp_types.h /usr/local/include/lp_solve/lp_types.h;\
	cp /tmp/lp_solve_5.5/lp_utils.h /usr/local/include/lp_solve/lp_utils.h;

#Patch build
COPY lpsolve55_patched_ccc /tmp/lp_solve_5.5/lpsolve55
#Build lib
RUN cd /tmp/lp_solve_5.5/lpsolve55 && sh lpsolve55_patched_ccc;\
	mkdir -p /usr/local/lib/lp_solve;\
	cp /tmp/lp_solve_5.5/lpsolve55/bin/ux64/liblpsolve55.so /usr/local/lib/liblpsolve55.so;\
	cp /tmp/lp_solve_5.5/lpsolve55/bin/ux64/liblpsolve55.a /usr/local/lib/liblpsolve55.a;
RUN rm -rf /tmp/lp_solve_5.5

############################################################################################
# Build RDkit
############################################################################################
# download rdkit
RUN cd /opt;\
	wget -O rdkit.tar.gz https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz;\
	tar xvzf rdkit.tar.gz;\
	rm rdkit.tar.gz

# build and install rdkit
RUN cd /opt/rdkit-Release_2017_09_3;\
	mkdir build;\
	cd build;\
	# DRDK_OPTIMIZE_NATIVE set to OFF since we are building for docker
	cmake .. -DRDK_PGSQL_STATIC=OFF\
		 	  -DRDK_BUILD_PYTHON_WRAPPERS=OFF\
			  -DRDK_BUILD_CPP_TESTS=OFF  \
			  -DRDK_BUILD_DESCRIPTORS3D=OFF\
			  -DRDK_INSTALL_STATIC_LIBS=OFF\
			  -DRDK_INSTALL_INTREE=OFF\
			  -DRDK_BUILD_INCHI_SUPPORT=ON\ 
			  -DRDK_OPTIMIZE_NATIVE=OFF\
			  -DCMAKE_CXX_STANDARD=11\
			  -DCMAKE_BUILD_TYPE=Release;\
	make install -j 4;\
	mkdir -p /usr/local/include/rdkit/External/INCHI-API/;\
	cp /opt/rdkit-Release_2017_09_3/External/INCHI-API/*.h /usr/local/include/rdkit/External/INCHI-API/;\
	cd /opt && rm -rf ./rdkit-Release_2017_09_3;

############################################################################################
# Build CFM-ID 4 MSML
############################################################################################
ARG BUILD_CFM_TRAIN="OFF"
ARG BUILD_CFM_TEST="OFF"
ARG MODEL_DIR="/trained_models_cfmid4.0"

# build and install cfm
RUN cd /opt;\
	wget -O cfmid.tar.gz https://bitbucket.org/wishartlab/cfm-id-code/get/CFM-ID_4.0.5.tar.gz;\
	mkdir cfmid;\
	tar xf cfmid.tar.gz -C cfmid;\
	rm cfmid.tar.gz;\
	mv cfmid/wishartlab* cfmid/v4;\
	mv cfmid/v4/cfm cfm;\
	mv cfmid/v4/cfm-pretrained-models ${MODEL_DIR};\
	rm -rf cfmid;
	
RUN cd /opt/cfm;\
	mkdir build;\
	cd build;\
	cmake .. -DINCLUDE_TESTS=${BUILD_CFM_TEST}\
			 -DINCLUDE_TRAIN=${BUILD_CFM_TRAIN}\
			 -DLPSOLVE_INCLUDE_DIR=/usr/local/include/lp_solve\
			 -DLPSOLVE_LIBRARY_DIR=/usr/local/lib\
			 -DCMAKE_CXX_STANDARD=11;\
	make install -j 4 && make clean;

############################################################################################
# Build CFM-ID 4 MSRB
############################################################################################
ENV PATH /opt/msrb/bin:$PATH
#RUN  apt-get update;\
#     apt-get install -y openjdk-11-jre-headless;
RUN apk update && apk add --no-cache openjdk11-jre-headless;
RUN  mkdir /opt/msrb;\
     cd /opt/msrb;\
	 wget -O msrb.tar.gz https://bitbucket.org/wishartlab/msrb-fragmenter/downloads/msrb-fragmenter-1.1.0.jar;

########################################################################
# clean Up
########################################################################
RUN apk del wget tar g++ gcc make cmake boost-dev;
RUN apk add --no-cache boost-filesystem boost-regex boost-serialization boost-system boost-thread;
RUN rm -rf /tmp/* /var/tmp/*
#RUN apk add --no-cache boost-program-options boost-test;

############################################################################################
# output folder
############################################################################################
RUN  mkdir /root/output;