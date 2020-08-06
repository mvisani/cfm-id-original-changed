FROM ubuntu:18.04
LABEL author="Fei Wang" 
LABEL name="CFM-ID 3 Docker Image"
LABEL image_version = "3.0.0.0"
LABEL cfm_msml_version="2.4"
LABEL msrb_version="1.0.6"
LABEL org.label-schema.build-date=$BUILD_DATE

###########################################################################################
# INSTALL CFM-ID 2 MSML
###########################################################################################
ENV CFM_ROOT /opt/cfm
ENV RDBASE /opt/rdkit-Release_2017_09_3
ENV LD_LIBRARY_PATH /opt/rdkit-Release_2017_09_3/lib:/opt/lp_solve/lpsolve55/bin/ux64:/usr/lib64:$LD_LIBRARY_PATH
ENV PATH /opt/cfm/bin:$PATH

############################################################################################
# Installs Boost and other apks
############################################################################################
RUN apt-get update;\
	apt-get install -y wget tar g++ gcc mpich cmake;\
	apt-get install -y 	libboost-dev\
						libboost-filesystem-dev \
						libboost-regex-dev\
						libboost-serialization-dev\
						libboost-system-dev\
						libboost-test-dev\
						libboost-thread-dev\
						libboost-program-options-dev;
						
############################################################################################
# Installs an LP Solver
############################################################################################
RUN mkdir /usr/lib64;
RUN	cd /usr/lib64;\
	wget -O lp_solve.tar.gz https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.5/lp_solve_5.5.2.5_dev_ux64.tar.gz/download;\
	tar xvzf lp_solve.tar.gz;\
	rm lp_solve.tar.gz

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
	cmake .. -DRDK_PGSQL_STATIC=OFF\
		 	  -DRDK_BUILD_PYTHON_WRAPPERS=OFF\
			  -DRDK_BUILD_CPP_TESTS=OFF  \
			  -DRDK_BUILD_DESCRIPTORS3D=OFF\
			  -DRDK_INSTALL_STATIC_LIBS=OFF\
			  -DRDK_INSTALL_INTREE=ON\
			  -DRDK_BUILD_INCHI_SUPPORT=ON\
			  -DRDK_OPTIMIZE_NATIVE=ON\
			  -DCMAKE_CXX_STANDARD=11;\
	make install -j 4

############################################################################################
# Build CFM-ID 3 MSML
############################################################################################
ARG BUILD_CFM_TRAIN="OFF"
ARG BUILD_CFM_TEST="OFF"

# build and install cfm
RUN cd /opt;\
	wget -O cfmid.tar.gz https://bitbucket.org/wishartlab/cfm-id-code/get/CFM-ID_2.4.tar.gz;\
	mkdir cfmid;\
	tar xf cfmid.tar.gz -C cfmid;\
	rm cfmid.tar.gz;\
	mv cfmid/wishartlab* cfmid/v2;\
	mv cfmid/v2/cfm cfm;\
	rm -rf cfmid;
	
RUN cd /opt/cfm;\
	mkdir build;\
	cd build;\
	cmake .. -DINCLUDE_TESTS=${BUILD_CFM_TEST}\
			 -DINCLUDE_TRAIN=${BUILD_CFM_TRAIN}\
			 -DLPSOLVE_INCLUDE_DIR=/usr/lib64 -DLPSOLVE_LIBRARY_DIR=/usr/lib64 \
			 -DCMAKE_CXX_STANDARD=11;\
	make install -j 4 && make clean;

RUN cd /;\
	wget -O cfmid2_trained_models.tar.gz  https://bitbucket.org/wishartlab/cfm-id-code/downloads/cfmid2_trained_models.tar.gz;\
	tar xf cfmid2_trained_models.tar.gz;

############################################################################################
# Build CFM-ID 3 MSRB
############################################################################################
ENV PATH /opt/msrb/bin:$PATH
RUN  apt-get update;\
     apt-get install -y openjdk-11-jre-headless;\
     mkdir /opt/msrb;\
     cd /opt/msrb;\
	 wget -O msrb.tar.gz https://bitbucket.org/wishartlab/msrb-fragmenter/downloads/msrb-fragmenter-1.0.8.jar;

########################################################################
# clean Up
########################################################################
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /tmp/* /var/tmp/*

############################################################################################
# output folder
############################################################################################
RUN  mkdir /root/output;