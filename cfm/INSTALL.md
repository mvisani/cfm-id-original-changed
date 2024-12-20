# INSTALL GUIDE

* [COMPILING FOR LINUX](#compling-for-linux)
* [COMPILING FOR MACOS](#compling-for-mac)
* [COMPILING FOR WINDOWS](#compling-for-windows)
* [Build On ComputeCanada Cedar Cluster](#build-on-computecanada-cedar-cluster)

## COMPILING FOR LINUX or WSL on Windows
For WSL check out [https://learn.microsoft.com/en-us/windows/wsl/install]

### Install Cmake and GCC

Install CMake VERSION 3.7.0 or higher, and GCC 7 or higher. Eariler GCC version may still work.

### Get Boost library

Install Boost VERSION 1.62 or higher. At minimum, include the regex,serialization,filesystem,system,thread,program_options,test modules.
You can install this from prebuild repo via ppa or your favorite package managment tool, here is an example on ubuntu:
** Install Boost VERSION 1.73 would not work RDKit_2017_09_3 **

```bash
sudo apt-get install libboost-all-dev
```

If you wish to install this from source code
Download boost_1_71_0.tar.gz from <http://www.boost.org/users/history/version_1_71_0.html>

```bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.71.0/source/boost_1_71_0.tar.gz
tar -xvzf boost_1_71_0.tar.gz
cd boost_1_71_0

BOOST_ROOT=/opt/boost_1_71_0
./bootstrap.sh --prefix=${BOOST_ROOT} --includedir=headers --libdir=dist
sudo ./b2 --prefix=${BOOST_ROOT} address-model=64 install
```

### Get RDKit library

Install RDKit (see <http://rdkit.org/>), including the InChI Extensions (python extensions are not required). Set the environment variable RDBASE to the RDKit install location (the directory with Code, lib etc..).
Download RDKit_2017_09_3.tgz from <https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz>;

*NOTE:Newer RDKit may work but we have not test it yet*

#### Install to insource location
```bash
wget https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz;
tar -zxvf Release_2017_09_3.tar.gz
cd ./rdkit-Release_2017_09_3
mkdir build
cd build
cmake .. -DRDK_PGSQL_STATIC=OFF -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_CPP_TESTS=OFF -DRDK_BUILD_DESCRIPTORS3D=OFF -DRDK_INSTALL_STATIC_LIBS=OFF-DRDK_INSTALL_INTREE=ON -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_OPTIMIZE_NATIVE=ON -DCMAKE_CXX_STANDARD=14 -DCMAKE_BUILD_TYPE=Release
sudo make install -j8
```
#### Install to Custom location
If You have a custom boost location such as ``` BOOST_ROOT=/opt/boost_1_71_0 ``` and rdkit install location such as ``` RDBASE=/opt/rdkit-Release_2017_09_3 ```:
```bash
BOOST_ROOT=/opt/boost_1_71_0
RDBASE=/opt/rdkit-Release_2017_09_3
cmake .. -DRDK_PGSQL_STATIC=OFF -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_CPP_TESTS=OFF -DRDK_BUILD_DESCRIPTORS3D=OFF -DRDK_INSTALL_STATIC_LIBS=OFF-DRDK_INSTALL_INTREE=ON -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_OPTIMIZE_NATIVE=ON -DCMAKE_CXX_STANDARD=14 -DCMAKE_BUILD_TYPE=Release  -DRDK_INSTALL_INTREE=OFF -DBoost_NO_BOOST_CMAKE=ON -DBOOST_ROOT=${BOOST_ROOT} -DCMAKE_INSTALL_PREFIX=${RDBASE}
sudo make install -j8
```

**Note** that ```-DRDK_INSTALL_INTREE=ON``` will install RDKit lib within its source file, while ```-DRDK_INSTALL_INTREE=OFF``` will install RDKit in the ```/usr/local/``` by default or to ```DCMAKE_INSTALL_PREFIX``` if set. However, RDKit will not automaticlly install  InChI Extension in the  ```/usr/local/```. You can move InChI Extension with (assume you are inside ```build ``` directory):

```bash
sudo mkdir  -p  /usr/local/include/rdkit/External/INCHI-API/;\
sudo cp  ../External/INCHI-API/*.h  /usr/local/include/rdkit/External/INCHI-API/;\
```

In case you have a different install location such as ``` RDBASE=/opt/rdkit-Release_2017_09_3 ```
```bash
sudo mkdir  -p  ${RDBASE}/include/rdkit/External/INCHI-API/;\
sudo cp  ../External/INCHI-API/*.h  ${RDBASE}/include/rdkit/External/INCHI-API/;\
```

### Get LPSolve library

You may be able to use one of the pre-compiled dev versions: <https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_dev_ux64.tar.gz/download> .

* NOTE : Prebuild library does not work on Apline linux because Apline does not use glibc (if you have no idea what this is, this most likey does not apply to you)

```bash
cd /tmp;
wget -O lp_solve_dev_ux64.tar.gz 'https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_dev_ux64.tar.gz/download' ;
mkdir lp_solve_dev_ux64
tar xvzf lp_solve_dev_ux64.tar.gz -C lp_solve_dev_ux64;
rm lp_solve_dev_ux64.tar.gz;

LPSOLVE_INCLUDE_DIR=/opt/lp_solve/
LPSOLVE_LIB_DIR=/opt/lp_solve/

sudo mkdir -p ${LPSOLVE_INCLUDE_DIR};
sudo mkdir -p ${LPSOLVE_LIB_DIR};
sudo cp ./lp_solve_dev_ux64/lp_Hash.h ${LPSOLVE_INCLUDE_DIR}/lp_Hash.h;
sudo cp ./lp_solve_dev_ux64/lp_SOS.h ${LPSOLVE_INCLUDE_DIR}/lp_SOS.h;
sudo cp ./lp_solve_dev_ux64/lp_lib.h ${LPSOLVE_INCLUDE_DIR}/lp_lib.h;
sudo cp ./lp_solve_dev_ux64/lp_matrix.h ${LPSOLVE_INCLUDE_DIR}/lp_matrix.h;
sudo cp ./lp_solve_dev_ux64/lp_mipbb.h ${LPSOLVE_INCLUDE_DIR}/lp_mipbb.h;
sudo cp ./lp_solve_dev_ux64/lp_types.h ${LPSOLVE_INCLUDE_DIR}/lp_types.h;
sudo cp ./lp_solve_dev_ux64/lp_utils.h ${LPSOLVE_INCLUDE_DIR}/lp_utils.h;

sudo cp ./lp_solve_dev_ux64/liblpsolve55.so ${LPSOLVE_LIB_DIR}/liblpsolve55.so;
sudo cp ./lp_solve_dev_ux64/liblpsolve55.a ${LPSOLVE_LIB_DIR}/liblpsolve55.a;
rm -rf ./lp_solve_dev_ux64;
```

If you wish to build your own, download compile the source code for LPSolve. Download lp_solve_5.5.2.11_source.tar.gz from <https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11>

```bash
cd /tmp;
wget -O lp_solve_src.tar.gz 'https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz/download' ;
tar xvzf lp_solve_src.tar.gz;
rm lp_solve_src.tar.gz;


cd /tmp/lp_solve_5.5;
LPSOLVE_INCLUDE_DIR=/opt/lp_solve/
LPSOLVE_LIB_DIR=/opt/lp_solve/
sudo mkdir -p ${LPSOLVE_INCLUDE_DIR};
sudo mkdir -p ${LPSOLVE_LIB_DIR};

sudo cp /tmp/lp_solve_dev_ux64/lp_Hash.h ${LPSOLVE_INCLUDE_DIR}/lp_Hash.h;
sudo cp /tmp/lp_solve_dev_ux64/lp_SOS.h ${LPSOLVE_INCLUDE_DIR}/lp_SOS.h;
sudo cp /tmp/lp_solve_dev_ux64/lp_lib.h ${LPSOLVE_INCLUDE_DIR}/lp_lib.h;
sudo cp /tmp/lp_solve_dev_ux64/lp_matrix.h ${LPSOLVE_INCLUDE_DIR}/lp_matrix.h;
sudo cp /tmp/lp_solve_dev_ux64/lp_mipbb.h ${LPSOLVE_INCLUDE_DIR}/lp_mipbb.h;
sudo cp /tmp/lp_solve_dev_ux64/lp_types.h ${LPSOLVE_INCLUDE_DIR}/lp_types.h;
sudo cp /tmp/lp_solve_dev_ux64/lp_utils.h ${LPSOLVE_INCLUDE_DIR}/lp_utils.h;

cd ./lpsolve55 && sh ccc;
cp ./bin/ux64/liblpsolve55.so ${LPSOLVE_LIB_DIR}/liblpsolve55.so;
cp ./bin/ux64/liblpsolve55.a ${LPSOLVE_LIB_DIR}/liblpsolve55.a;
rm -rf /tmp/lp_solve_5.5;
```

### Get MPI library (optional)

If Compiling the cfm-train and cfm-test executables, Install a version of MPI. Current CFM-ID are comptable with MPI-3.2

If library is not inistalled, on ubuntu you could also use

```Bash
sudo apt-get install libopenmpi-dev 
```

### Setup Libraries

if your libaray installation is not under the stanard ```/usr/``` directory, you will need to set ```LD_LIBRARY_PATH``` to include Boost, RDKit and LPSolve library locations. This can be done in **one of following method**,

1. Add following command in ```~/.bashrc```

    ```bash
   export BOOST_ROOT=/opt/boost_1_71_0
   export RDBASE=/opt/rdkit-Release_2017_09_3
   export LP_SOLVE=/opt/lp_solve
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BOOST_ROOT/lib:$RDBASE/lib:$LP_SOLVE
   ```

   Then reload by  ```source ~/.bashrc```

2. go to ```/etc/ld.so.conf.d``` add ```*.conf``` for each library Boost, RDKit and LPSolve library locations

```
> boost.conf
/opt/boost_1_71_0/lib
rdkit.conf
/opt/rdkit-Release_2017_09_3/lib
lp_solve.conf
/opt/lp_solve
```
Then reload ld by  ```sudo ldconfig```
### Build CFM-ID

Download CFM-ID Release version from <https://bitbucket.org/wishartlab/cfm-id-code/downloads/?tab=tags>, or clone git repo. Run cmake ```CFM_ROOT``` where ```CFM_ROOT``` is the location of the cfm directory
e.g. if you are in cfm/build, you can use cmake .. , setting the ```LPSOLVE_INCLUDE_DIR``` and ```LPSOLVE_LIBRARY_DIR``` values appropriately. Use  ```INCLUDE_TESTS``` and ```INCLUDE_TRAIN``` to enable or disbale```cfm-tests``` and ```cfm-train```

```bash
    cd cfm
    mkdir  build;\
    cd  build;\
    cmake  .. \
        -DLPSOLVE_INCLUDE_DIR=${LP_SOLVE}\
        -DLPSOLVE_LIBRARY_DIR=${LP_SOLVE}\
        -DBOOST_ROOT=${BOOST_ROOT}\
        -DBoost_NO_BOOST_CMAKE=ON \
        -DRDKIT_INCLUDE_DIR=${RDBASE}/include/rdkit\
        -DRDKIT_INCLUDE_EXT_DIR=${RDBASE}/include/rdkit/External\
        -DRDKIT_LIBRARY_DIR=${RDBASE}/lib\
        -DCMAKE_CXX_STANDARD=14 \
        -DINCLUDE_TESTS=${BUILD_CFM_TEST}\
        -DINCLUDE_TRAIN=${BUILD_CFM_TRAIN};
    make install
```

This should produce the executable files in the ```cfm/bin``` directory.  Change to this directory. Run the programs from a command line as detailed on <https://sourceforge.net/p/cfm-id/wiki/Home/>
(Note: replace ~ with the paths where you've installed Boost or RDKit or lpsolve respectively.)

## COMPILING FOR MACOS

### Install Homebrew

Some packages such as libboost and mpi can be installed via homebrew: [#<https://brew.sh/>]

### Install LPSolve (Required if build CFM-Train or CFM-Test)

```bash
brew install lp_solve
```

### Compile OpenMPI (Required if build CFM-Train or CFM-Test)

By Default OpenMPI build does not include C++ binding. That is openmpi installed from homebrew will not work.

```bash
configure --enable-mpi-cxx --disable-mpi-fortran --prefix=/usr/local/
make install all
```

### Compile Boost lib

Download boost-1.71.0 and assume we are going to install it in the ```/opt/boost_1_17_0```

```
./bootstrap.sh --prefix=/opt/boost_1_17_0 --without-libraries=python
sudo ./b2 cxxflags=-std=c++17 cxxflags="-stdlib=libc++" linkflags="-stdlib=libc++" -j 8 install
```

### Compile Rdkit

Download Rdkit-2017_09_3 and assume we are going to install it in the ```/opt/RDKit_2017_09_3```

```
   wget https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz;
   tar -zxvf Release_2017_09_3.tar.gz
   cd ../Release_2017_09_3
   mkdir build
   cd build
   cmake .. -DRDK_PGSQL_STATIC=OFF -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_CPP_TESTS=OFF -DRDK_BUILD_DESCRIPTORS3D=OFF -DRDK_INSTALL_STATIC_LIBS=OFF   -DRDK_INSTALL_INTREE=OFF -DCMAKE_INSTALL_PREFIX=/opt/RDKit_2017_09_3 -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_OPTIMIZE_NATIVE=ON -DCMAKE_CXX_STANDARD=11 -DCMAKE_BUILD_TYPE=Release -DCMAKE_MACOSX_RPATH=TRUE
   make install
```

Note that  ```-DRDK_INSTALL_INTREE=ON``` install RDKit lib within its source file, while ```-DRDK_INSTALL_INTREE=OFF``` will install RDKit in the ```/usr/local/``` or specficed by ```DCMAKE_INSTALL_PREFIX``` However, RDKit will not automaticlly install  InChI Extension in the  ```/usr/local/```. You can move InChI Extension with:

```
sudo mkdir  -p /opt/RDKit_2017_09_3/include/rdkit/External/INCHI-API\
sudo cp  ../External/INCHI-API/*.h  /opt/RDKit_2017_09_3/include/rdkit/External;\
```

### Setup Libraries

Set ```DYLD_LIBRARY_PATH``` to include Boost, RDKit and LPSolve library locations.

```
export DYLD_LIBRARY_PATH=/opt/boost_1_17_0/lib:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/opt/RDKit_2017_09_3/lib:$DYLD_LIBRARY_PATH
```

### Build CFM-ID

Download CFM-ID Release version from <https://bitbucket.org/wishartlab/cfm-id-code/downloads/?tab=tags>, or clone git repo. Run cmake ```CFM_ROOT``` where ```CFM_ROOT``` is the location of the cfm directory
e.g. if you are in cfm/build, you can use cmake .. , setting the ```LPSOLVE_INCLUDE_DIR``` and ```LPSOLVE_LIBRARY_DIR``` values appropriately. Use  ```INCLUDE_TESTS``` and ```INCLUDE_TRAIN``` to enable or disbale```cfm-tests``` and ```cfm-train```

Following command assume libraries are installed via Homebrew

```
    cd cfm
    mkdir  build;
    cd  build;
    cmake ..  -DINCLUDE_TESTS=ON -DINCLUDE_TRAIN=ON -DLPSOLVE_INCLUDE_DIR=/usr/local/Cellar/lp_solve/5.5.2.11/include -DLPSOLVE_LIBRARY_DIR=/usr/local/Cellar/lp_solve/5.5.2.11/lib -DRDKIT_INCLUDE_DIR=/opt/RDKit_2017_09_3/include/rdkit/ -DRDKIT_INCLUDE_EXT_DIR=/opt/RDKit_2017_09_3/include/rdkit/External -DRDKIT_LIBRARY_DIR=/opt/RDKit_2017_09_3/lib -DCMAKE_CXX_STANDARD=11
    make install
```

This should produce the executable files in the ```cfm/bin``` directory.  Change to this directory. Run the programs from a command line as detailed on <https://sourceforge.net/p/cfm-id/wiki/Home/>
(Note: replace ~ with the paths where you've installed Boost or RDKit or lpsolve respectively.)

------------------

## COMPILING FOR WINDOWS

### Build on windows has NOT been verfied on current version

1. Install CMake.

2. Install Boost (see www.boost.org). At minimum, include the filesystem, system and serialization modules.
   Set an environment variable BOOST_ROOT to the Boost install location.

3. Install RDKit (see <http://rdkit.org/>), including the InChI Extensions (python extensions are not required).
   Set the environment variable RDBASE to the RDKit install location (the directory with Code, lib etc..)..

4. Download and unzip a development version of LPSolve (e.g. lp_solve_5.5.2.0_dev_win32.zip
    * see <https://sourceforge.net/projects/lpsolve>).

5. (optional - if compiling the cfm-train and cfm-test executables)
   Install a version of MPI e.g. Microsoft MPI. <https://www.microsoft.com/en-us/download/details.aspx?id=100593>

6. Start the CMake GUI and set the source code location to the cfm directory (the directory with cfm-code, cfm-id...etc).
   Click Configure. A pop-up should appear asking you to select the generator.
   This code has been tested with VisualStudio 10 (using the free VisualStudio Express 2010 edition) so this is recommended.

7. Update the LPSOLVE_INCLUDE_DIR to the root directory of LPSolve (i.e. where lp_lib.h is) and LPSOLVE_LIBRARY_DIR
   to the same directory (i.e. where liblpsolve55.dll is).

8. If you want to compile the cfm-train and cfm-test modules, click the INCLUDE_TRAIN and INCLUDE_TESTS
   checkboxes respectively. Otherwise make sure these are unchecked.

9. Once configration is complete, click Generate. This should generate the relevant project or makefiles.
    For Visual Studio, cfm.sln will be generated. Open this file in Visual Studio and build the INSTALL project. Any other generator, you're on your own!

10. This should produce the executable files in the cfm/bin directory.  Either add this directory to your path
    or start a command prompt in this directory. Run them from a command line as detailed
    on <https://sourceforge.net/p/cfm-id/wiki/Home/>.

------------------

## Build On ComputeCanada Cedar Cluster

### Use System Packages

Use system provide compiler, boost, cmake, and mpi

```

module load StdEnv/2018.3
module load nixpkgs/16.09
module load intel/2018.3
module load boost/1.68.0
module load openmpi/3.1.2
module load cmake/3.12.3
```

## go to home

```cd $HOME
mkdir libs
cd libs
```

### Boost and Cmake

Use -DBoost_NO_BOOST_CMAKE=ON to avoid problem
As for cmake 3.16.3 and Boost 1.70 <https://gitlab.kitware.com/cmake/cmake/blob/master/Modules/FindBoost.cmake>
"If Boost was built using the boost-cmake project or from Boost 1.70.0 on
it provides a package configuration file for use with find_package's config mode.
This module looks for the package configuration file called
``BoostConfig.cmake`` or ``boost-config.cmake`` and stores the result in
``CACHE`` entry "Boost_DIR".  If found, the package configuration file is loaded
and this module returns with no further action.  See documentation of
the Boost CMake package configuration for details on what it provides.

Set ``Boost_NO_BOOST_CMAKE`` to ``ON``, to disable the search for boost-cmake."

### Get lp_solve

Install lp solve in the project directory or home directory. Precompiled version should work out of box at <https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_dev_ux64.tar.gz/download>

```
mkdir lpsolve
cd lpsolve
wget https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_dev_ux64.tar.gz
tar -xzf lp_solve_5.5.2.11_dev_ux64.tar.gz
rm lp_solve_5.5.2.11_dev_ux64.tar.gz
cd ..
```

### Get RDkit

Install RDkit in the project directory or home directory. Make sure set ``Boost_NO_BOOST_CMAKE`` to ``ON``  and
set ``DRDK_INSTALL_INTREE`` to ``ON``

```
wget -O rdkit.tar.gz https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz;\
tar xvzf rdkit.tar.gz;
rm rdkit.tar.gz;
cd rdkit-Release_2017_09_3;
mkdir build;
cd build;
cmake .. -DRDK_PGSQL_STATIC=OFF \
     -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
     -DRDK_BUILD_CPP_TESTS=OFF \
     -DRDK_BUILD_DESCRIPTORS3D=OFF \
     -DRDK_INSTALL_STATIC_LIBS=OFF \
     -DRDK_BUILD_INCHI_SUPPORT=ON \
     -DRDK_OPTIMIZE_NATIVE=ON \
     -DCMAKE_CXX_STANDARD=14 \
     -DCMAKE_BUILD_TYPE=Release \
     -DRDK_INSTALL_INTREE=ON \
     -DBoost_NO_BOOST_CMAKE=ON\
     -DBOOST_ROOT=$BOOST_ROOT;
make install
```

### Build CFM-ID

```
cmake .. -DINCLUDE_TESTS=OFF\
  -DINCLUDE_TRAIN=ON \
  -DLPSOLVE_INCLUDE_DIR=/$HOME/libs/lpsolve\
  -DLPSOLVE_LIBRARY_DIR=/$HOME/libs/lpsolve \
  -DRDKIT_INCLUDE_DIR=/$HOME/libs/rdkit_2017_09_3/Code \
  -DRDKIT_INCLUDE_EXT_DIR=/$HOME/libs/rdkit_2017_09_3/External \
  -DRDKIT_LIBRARY_DIR=/$HOME/libs/rdkit_2017_09_3/lib \
  -DCMAKE_CXX_STANDARD=14 \
  -DBOOST_ROOT=$BOOST_ROOT\
  -DBoost_NO_BOOST_CMAKE=ON\
  -DCMAKE_BUILD_TYPE=Release 
make install
```
