FROM ubuntu:18.04
MAINTAINER Fei Wang

ENV RDBASE /opt/rdkit-Release_2017_09_3
ENV CFM_ROOT /opt/cfm
ENV LD_LIBRARY_PATH /opt/rdkit-Release_2017_09_3/lib:/opt/lp_solve/lpsolve55/bin/ux64:/usr/lib64:$LD_LIBRARY_PATH
ENV PATH /opt/cfm/bin:/opt/miniconda/bin:$PATH

RUN apt-get update && apt-get install -y wget tar build-essential mpich libboost-all-dev cmake

# Installs an LP Solver
RUN mkdir /usr/lib64;\
	cd /usr/lib64;\
	wget -O lp_solve.tar.gz https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.0/lp_solve_5.5.2.0_dev_ux64.tar.gz/download;\
	tar xvzf lp_solve.tar.gz

# download rdkit
RUN cd /opt;\
	wget -O rdkit.tar.gz https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz;\
	tar xvzf rdkit.tar.gz;\
	rm rdkit.tar.gz

# build and install rdkit
RUN cd /opt/rdkit-Release_2017_09_3;\
	mkdir build;\
	cd build;\
	cmake .. -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_INCHI_SUPPORT=ON -DCMAKE_CXX_STANDARD=11;\
	make install -j4

# build and install cfm
RUN cd /opt;\
	wget -O cfmid.tar.gz https://bitbucket.org/7FeiW/cfm-id-code/get/CFM-ID_4.0.1.tar.gz;\
	mkdir cfmid;\
	tar xf cfmid.tar.gz -C cfmid;\
	mv cfmid/7FeiW* cfmid/fei;\
	mv cfmid/fei/cfm cfm;\
	rm -rf cfmid;\
	cd cfm;\
	mkdir build;\
	cd build;\
	cmake .. -DINCLUDE_TESTS=OFF -DINCLUDE_TRAIN=ON -DLPSOLVE_INCLUDE_DIR=/usr/lib64 -DLPSOLVE_LIBRARY_DIR=/usr/lib64 \
	-DCMAKE_CXX_STANDARD=11;\
	make install -j4

RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

