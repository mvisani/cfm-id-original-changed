
## COMPILING FOR LINUX

###  Install Cmake and GCC
Install CMake VERSION 3.7.0 or higher, and GCC 7 or higher. Eariler GCC version may still work.

### Get Boost library
Install Boost VERSION 1.62 or higher. At minimum, include the regex,serialization,filesystem,system,thread,program_options,test modules.
You can install this from prebuild repo via ppa or your favorite package managment tool, here is an example on ubuntu:
```
sudo apt-get install libboost-all-dev
```
If you wish to install this from source code
Download boost_1_62_0.tar.gz from http://www.boost.org/users/history/version_1_62_0.html
```
cd boost_1_62_0
./bootstrap.sh --prefix=. \
                    --with-libraries=regex,\
                    serialization,\
                    filesystem,\
                    system,\
                    thread,
                    program_options,\
                    test 
./b2 address-model=64 install
```
### Get RDKit library
Install RDKit (see http://rdkit.org/), including the InChI Extensions (python extensions are not required). Set the environment variable RDBASE to the RDKit install location (the directory with Code, lib etc..).
Download RDKit_2017_09_3.tgz from https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz;

*NOTE:Newer RDKit may work but we have not test it yet*

```
   tar -zxvf RDKit_2017_09_3.tgz
   cd ../..
   mkdir build
   cd build
   cmake .. \ -DRDK_PGSQL_STATIC=OFF\  -DRDK_BUILD_PYTHON_WRAPPERS=OFF\   -DRDK_BUILD_CPP_TESTS=OFF -DRDK_BUILD_DESCRIPTORS3D=OFF\ -DRDK_INSTALL_STATIC_LIBS=OFF   -DRDK_INSTALL_INTREE=ON \ -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_OPTIMIZE_NATIVE=ON  \ -DCMAKE_CXX_STANDARD=11 \  -DCMAKE_BUILD_TYPE=Release
   make install
```
Note that ```-DRDK_INSTALL_INTREE=ON``` will install RDKit lib within its source file, while ```-DRDK_INSTALL_INTREE=OFF``` will install RDKit in the ```/usr/local/```. However, RDKit will not automaticlly install  InChI Extension in the  ```/usr/local/```. You can move InChI Extension with:
```
mkdir  -p  /usr/local/include/rdkit/External/INCHI-API/;\
cp  ../External/INCHI-API/*.h  /usr/local/include/rdkit/External/INCHI-API/;\
```
### Get LPSolve library
You may be able to use one of the pre-compiled dev versions: https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.5/lp_solve_5.5.2.5_dev_ux64.tar.gz/download

If you wish to build your own, download compile the source code for LPSolve. Download lp_solve_5.5.2.5_source.tar.gz from https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.5
```
    tar -zxvf lp_solve_5.5.2.5_source.tar.gz
    cd lp_solve_5.5/lpsolve55
    ./ccc
```
This should create libs in e.g. lp_solve_5.5/lpsolve55/bin/ux64.
If you encounter a build error with ```./ccc```, please use our  patched version at: https://bitbucket.org/wishartlab/cfm-id-code/downloads/lpsolve55_patched_ccc

### Get MPI library (optional)
If Compiling the cfm-train and cfm-test executables, Install a version of MPI.

###  Setup Libraries
if your libaray installation is not under the stanard ```/usr/``` directory, you will need to set ```LD_LIBRARY_PATH``` to include Boost, RDKit and LPSolve library locations. This can be done in one of following method,
1. Add following command in ```~/.bashrc```
    ```
    export LD_LIBRARY_PATH = $LD_LIBRARY_PATH:~/boost_1_62_0/lib:~/RDKit_2017_09_3/lib:~/lp_solve_5.5/lpsolve55/bin/ux64
    ```
   Then reload by  ```source ~/.bashrc ```
2. go to ```/etc/ld.so.conf.d``` add ```*.conf``` for each library Boost, RDKit and LPSolve library locations
    > boost.conf
      ~/boost_1_64_0/lib
      rdkit.conf
      ~/RDKit_2013_09_1/lib
      lp_solve.conf
      ~/lp_solve_5.5/lpsolve55/bin/ux64
        
   Then reload ld by  ```sudo ldconfig ```

### Build CFM-ID
Download CFM-ID Release version from https://bitbucket.org/wishartlab/cfm-id-code/downloads/?tab=tags, or clone git repo. Run cmake ```CFM_ROOT ``` where ```CFM_ROOT``` is the location of the cfm directory
e.g. if you are in cfm/build, you can use cmake .. , setting the ```LPSOLVE_INCLUDE_DIR``` and ```LPSOLVE_LIBRARY_DIR``` values appropriately. Use  ```INCLUDE_TESTS``` and ```INCLUDE_TRAIN``` to enable or disbale``` cfm-tests``` and ``` cfm-train``` 

```
    cd cfm
    mkdir  build;\
    cd  build;\
    cmake  ..  
        -DINCLUDE_TESTS=${BUILD_CFM_TEST}\
        -DINCLUDE_TRAIN=${BUILD_CFM_TRAIN}\
        -DLPSOLVE_INCLUDE_DIR=/usr/local/include/lp_solve\
        -DLPSOLVE_LIBRARY_DIR=/usr/local/lib\
        -DRDKIT_INCLUDE_DIR=/usr/local/include/rdkit\
        -DRDKIT_INCLUDE_EXT_DIR=/usr/local/include/rdkit/External\
        -DRDKIT_LIBRARY_DIR=/usr/local/lib\
        -DCMAKE_CXX_STANDARD=11;\
    make install
```

This should produce the executable files in the ```cfm/bin``` directory.  Change to this directory. Run the programs from a command line as detailed on https://sourceforge.net/p/cfm-id/wiki/Home/
(Note: replace ~ with the paths where you've installed Boost or RDKit or lpsolve respectively.)


COMPILING FOR MAC
------------------
Please refers to the same intrcution as for liunx with follow change
###  Install Libraries
Some packages such as libboost and mpi can be installed via homebrew: https://brew.sh/
###  Setup Libraries
Set ```DYLD_LIBRARY_PATH ``` to include Boost, RDKit and LPSolve library locations. 
```
export LD_LIBRARY_PATH = $LD_LIBRARY_PATH:~/boost_1_62_0/lib:~/RDKit_2017_09_3/lib:~/lp_solve_5.5/lpsolve55/bin/ux64
```
## COMPILING FOR WINDOWS 

### Build on windows has been verfied on current version

1. Install CMake.

2. Install Boost (see www.boost.org). At minimum, include the filesystem, system and serialization modules.
   Set an environment variable BOOST_ROOT to the Boost install location.

3. Install RDKit (see http://rdkit.org/), including the InChI Extensions (python extensions are not required).
   Set the environment variable RDBASE to the RDKit install location (the directory with Code, lib etc..)..

4. Download and unzip a development version of LPSolve (e.g. lp_solve_5.5.2.0_dev_win32.zip
    - see https://sourceforge.net/projects/lpsolve).

5. (optional - if compiling the cfm-train and cfm-test executables)
   Install a version of MPI e.g. Microsoft MPI.

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
    on https://sourceforge.net/p/cfm-id/wiki/Home/.
    
## Build On ComputeCanada Cedar Cluster

###  Use System Packages
Use system provide compiler, boost, cmake, and mpi
```
module load intel/2018.3
module load boost/1.68.0
module load openmpi/3.1.2
module load cmake/3.12.3
```

### Boost and Cmake 
Use -DBoost_NO_BOOST_CMAKE=ON to avoid problem
As for cmake 3.16.3 and Boost 1.70 https://gitlab.kitware.com/cmake/cmake/blob/master/Modules/FindBoost.cmake
"If Boost was built using the boost-cmake project or from Boost 1.70.0 on
it provides a package configuration file for use with find_package's config mode.
This module looks for the package configuration file called
``BoostConfig.cmake`` or ``boost-config.cmake`` and stores the result in
``CACHE`` entry "Boost_DIR".  If found, the package configuration file is loaded
and this module returns with no further action.  See documentation of
the Boost CMake package configuration for details on what it provides.

Set ``Boost_NO_BOOST_CMAKE`` to ``ON``, to disable the search for boost-cmake."

### Get lp_solve
Install lp solve in the project directory or home directory. Precompiled version should work out of box at https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_dev_ux64.tar.gz/download

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
     -DCMAKE_CXX_STANDARD=11 \
     -DCMAKE_BUILD_TYPE=Release \
     -DRDK_INSTALL_INTREE=ON \
     -DBOOST_ROOT=$BOOST_ROOT \
     -DBoost_NO_BOOST_CMAKE=ON;
```

### Build CFM-ID
```
cmake .. -DINCLUDE_TESTS=ON\
  -DINCLUDE_TRAIN=ON -DLPSOLVE_INCLUDE_DIR=/home/feiw/libs/lpsolve\
  -DLPSOLVE_LIBRARY_DIR=/home/feiw/libs/lpsolve \
  -DRDKIT_INCLUDE_DIR=/home/feiw/libs/rdkit_2017_09_3/Code \
  -DRDKIT_INCLUDE_EXT_DIR=/home/feiw/libs/rdkit_2017_09_3/External \
  -DRDKIT_LIBRARY_DIR=/home/feiw/libs/rdkit_2017_09_3/lib \
  -DCMAKE_CXX_STANDARD=11 \
  -DBOOST_ROOT=$BOOST_ROOT\
  -DBoost_NO_BOOST_CMAKE=ON;
```
