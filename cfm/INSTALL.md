# COMPILING FOR WINDOWS #
---------------------

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

8. (optional - if compiling the cfm-train and cfm-test executables)
   Update the LBFGS_INCLUDE_DIR and LIBFGS_LIBRARY_DIR variables to the locations of lbfgs.h and lbfgs.lib respectively.

9. If you want to compile the cfm-train and cfm-test modules, click the INCLUDE_TRAIN and INCLUDE_TESTS
   checkboxes respectively. Otherwise make sure these are unchecked.

10. Once configration is complete, click Generate. This should generate the relevant project or makefiles.
    For Visual Studio, cfm.sln will be generated. Open this file in Visual Studio and build the INSTALL project. Any other generator, you're on your own!

11. This should produce the executable files in the cfm/bin directory.  Either add this directory to your path
    or start a command prompt in this directory. Run them from a command line as detailed
    on https://sourceforge.net/p/cfm-id/wiki/Home/.

# COMPILING FOR LINUX #

1. Install CMake VERSION 3.7.0 or higher

2. Install Boost VERSION 1.62 or higher. At minimum, include the regex,serialization,filesystem,system,thread,program_options,test modules.

> e.g:
> Download boost_1_64_0.tar.gz from http://www.boost.org/users/history/version_1_64_0.html
> tar -zxvf boost_1_64_0.tar.gz
> cd boost_1_64_0
> ./bootstrap.sh --prefix=. --with-libraries=regex,serialization,filesystem,system,thread,program_options,test
>
> ./b2 address-model=64 cflags=-fPIC cxxflags=-fPIC install
> export BOOST_ROOT=~/boost_1_64_0

3. Install RDKit (see http://rdkit.org/), including the InChI Extensions (python extensions are not required). Set the environment variable RDBASE to the RDKit install location (the directory with Code, lib etc..).

> e.g.
> Download RDKit_2017_09_3.tgz from https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz;\
> tar -zxvf RDKit_2017_09_3.tgz
> cd ../..
> mkdir build
> cd build
> cmake .. -DRDK_PGSQL_STATIC=OFF -DRDK_BUILD_PYTHON_WRAPPERS=OFF  -DRDK_BUILD_CPP_TESTS=OFF  -DRDK_BUILD_DESCRIPTORS3D=OFF  -DRDK_INSTALL_STATIC_LIBS=OFF  > > > -DRDK_INSTALL_INTREE=ON  -DRDK_BUILD_INCHI_SUPPORT=ON  -DRDK_OPTIMIZE_NATIVE=ON  -DCMAKE_CXX_STANDARD=11  -DCMAKE_BUILD_TYPE=Release
>make install
> export RDBASE=~/RDKit_2017_09_3

4. Download and compile the source code for LPSolve. Note: you may be able to use one of the pre-compiled dev versions (e.g.lp_solve_5.5.2.0_dev_ux64.tar.gz) but compiling from source is probably more reliable in terms of getting a correct match.

e.g.
Download lp_solve_5.5.2.0_source.tar.gz from https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.0
tar -zxvf lp_solve_5.5.2.0_source.tar.gz
cd lp_solve_5.5/lpsolve55
./ccc
(should create libs in e.g. lp_solve_5.5/lpsolve55/bin/ux64)

5. (optional - if compiling the cfm-train and cfm-test executables)
  Install a version of MPI.

6. Download or check out the cfm code and create a new directory
where you want the build files to appear and move to that directory.

e.g.
svn checkout svn://svn.code.sf.net/p/cfm-id/code/cfm cfm
mkdir build
cd build

8. Run cmake CFM_ROOT  where CFM_ROOT is the location of the cfm directory
e.g. if you are in cfm/build, you can use cmake .. , setting the LPSOLVE_INCLUDE_DIR and LPSOLVE_LIBRARY_DIR values appropriately.

cmake .. -DLPSOLVE_INCLUDE_DIR=~/lp_solve_5.5 -DLPSOLVE_LIBRARY_DIR=~/lp_solve_5.5/lpsolve55/bin/ux64

9. (optional - if compiling the cfm-train and cfm-test executables), Use e.g.
> cmake .. -DINCLUDE_TESTS=ON -DINCLUDE_TRAIN=ON -DLPSOLVE_INCLUDE_DIR=~/lp_solve_5.5 -DLPSOLVE_LIBRARY_DIR=~/lp_solve_5.5/lpsolve55/bin/ux64

10. make install

11. This should produce the executable files in the cfm/bin directory.  Change to this directory.

12. Set LD_LIBRARY_PATH to include Boost, RDKit and LPSolve library locations
e.g.
> export LD_LIBRARY_PATH = $LD_LIBRARY_PATH:~/boost_1_64_0/lib:~/RDKit_2013_09_1/lib:~/lp_solve_5.5/lpsolve55/bin/ux64

Or go to /etc/ld.so.conf.d add *.conf for each library Boost, RDKit and LPSolve library locations
e.g.
> boost.conf
  ~/boost_1_64_0/lib
  rdkit.conf
  ~/RDKit_2013_09_1/lib
  lp_solve.conf
  ~/lp_solve_5.5/lpsolve55/bin/ux64
reload ld
> sudo ldconfig

13. (optional - if compiling the cfm-train and cfm-test executables) Also add the libLBFGS location.
> export LD_LIBRARY_PATH = $LD_LIBRARY_PATH:~/liblbfgs-1.10/bin/lib

14. Run the programs from a command line as detailed on https://sourceforge.net/p/cfm-id/wiki/Home/
(Note: replace ~ with the paths where you've installed Boost or RDKit or lpsolve respectively.)


COMPILING FOR MAC
------------------

1. Install CMake (or check it's already there by running cmake at the command line).

2. Install Boost (see www.boost.org). At minimum, include the filesystem, system and serialization modules. Set an environment variable BOOST_ROOT to the Boost install location.

e.g:
Download boost_1_64_0.tar.gz from http://www.boost.org/users/history/version_1_64_0.html
tar -zxvf boost_1_64_0.tar.gz
cd boost_1_64_0
./bootstrap.sh --prefix=. --with-libraries=regex,serialization,filesystem,system
./b2 address-model=64 cflags=-fPIC cxxflags=-fPIC install
export BOOST_ROOT=~/boost_1_64_0

3. Install RDKit (see http://rdkit.org/), including the InChI Extensions (python extensions are not required). Set the environment variable RDBASE to the RDKit install location (the directory with Code, lib etc..).

e.g.
Download RDKit_2013_09_1.tgz from https://sourceforge.net/projects/rdkit/files/rdkit/Q3_2013/
tar -zxvf RDKit_2013_09_1.tgz
cd RDKit_2013_09_1/External/INCHI-API
bash download-inchi.sh
cd ../..
mkdir build
cd build
cmake .. -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DBOOST_ROOT=~/boost_1_64_0 
make install
export RDBASE=~/RDKit_2013_09_1

4. Download and compile the source code for LPSolve.

e.g.
Download lp_solve_5.5.2.0_source.tar.gz from https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.0
tar -zxvf lp_solve_5.5.2.0_source.tar.gz
cd lp_solve_5.5/lpsolve55
./ccc.osx
(should create libs in e.g. lp_solve_5.5/lpsolve55/bin/osx64)

5. (optional - if compiling the cfm-train and cfm-test executables) Install a version of MPI.

6. (optional - if compiling the cfm-train and cfm-test executables) Download liblbfgs-1.10.tar.gz
from https://github.com/downloads/chokkan/liblbfgs/liblbfgs-1.10.tar.gz
> tar -zxvf liblbfgs-1.10.tar.gz
> cd liblbfgs-1.10
> ./configure
> make
> make install

7. Download or check out the cfm code and create a new directory where you want the build files
to appear and move to that directory.

e.g.
svn checkout svn://svn.code.sf.net/p/cfm-id/code/cfm cfm
mkdir build
cd build

8. Run cmake CFM_ROOT  where CFM_ROOT is the location of the cfm directory e.g. if you are in cfm/build, you can use cmake .. , setting the LPSOLVE_INCLUDE_DIR and LPSOLVE_LIBRARY_DIR values appropriately.

cmake .. -DLPSOLVE_INCLUDE_DIR=~/lp_solve_5.5 -DLPSOLVE_LIBRARY_DIR=~/lp_solve_5.5/lpsolve55/bin/osx64


9. (optional - if compiling the cfm-train and cfm-test executables), Use
e.g.
> cmake .. -D INCLUDE_TESTS=ON -D INCLUDE_TRAIN=ON -DLPSOLVE_INCLUDE_DIR=~/lp_solve_5.5 -DLPSOLVE_LIBRARY_DIR=~/lp_solve_5.5/lpsolve55/bin/osx64 -DLBFGS_INCLUDE_DIR=~/liblbfgs-1.10/bin/include -DLBFGS_LIBRARY_DIR=~/liblbfgs-1.10/bin/lib

10. make install

11. This should produce the executable files in the cfm/bin directory.  Change to this directory.

12. Set DYLD_LIBRARY_PATH to include Boost, RDKit and LPSolve library locations
e.g.
> export DYLD_LIBRARY_PATH = $DYLD_LIBRARY_PATH:~/boost_1_64_0/lib:~/RDKit_2013_09_1/lib:~/lp_solve_5.5/lpsolve55/bin/osx64

13. (optional - if compiling the cfm-train and cfm-test executables) Also add the libLBFGS location.
> export DYLD_LIBRARY_PATH = $DYLD_LIBRARY_PATH:~/liblbfgs-1.10/bin/lib

14. Run the programs from a command line as detailed on https://sourceforge.net/p/cfm-id/wiki/Home/
(Note: replace ~ with the paths where you've installed Boost or RDKit or lpsolve respectively.)

###### ComputeCanada Cedar Install Note ####################
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

#Get Boost
wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz
tar -zxvf boost_1_66_0.tar.gz
cd boost_1_66_0
./bootstrap.sh --prefix=. --with-libraries=regex,serialization,filesystem,system,thread,program_options,test
./b2 address-model=64 cflags=-fPIC cxxflags=-fPIC install
export BOOST_ROOT=$pwd

#Get Rdkit
wget -O rdkit.tar.gz https://github.com/rdkit/rdkit/archive/Release_2017_09_3.tar.gz;\
tar xvzf rdkit.tar.gz;
rm -rf * && cmake .. -DRDK_PGSQL_STATIC=OFF -DRDK_BUILD_PYTHON_WRAPPERS=OFF  -DRDK_BUILD_CPP_TESTS=OFF  -DRDK_BUILD_DESCRIPTORS3D=OFF  -DRDK_INSTALL_STATIC_LIBS=OFF   -DRDK_BUILD_INCHI_SUPPORT=ON  -DRDK_OPTIMIZE_NATIVE=ON  -DCMAKE_CXX_STANDARD=11 -DCMAKE_BUILD_TYPE=Release -DRDK_INSTALL_INTREE=ON -DBOOST_ROOT=$BOOST_ROOT -DBoost_NO_BOOST_CMAKE=ON;

# Build CFM-ID
cmake .. -DINCLUDE_TESTS=ON  -DINCLUDE_TRAIN=ON -DLPSOLVE_INCLUDE_DIR=/home/feiw/libs/lpsolve -DLPSOLVE_LIBRARY_DIR=/home/feiw/libs/lpsolve -DRDKIT_INCLUDE_DIR=/home/feiw/libs/rdkit_2017_09_3/Code -DRDKIT_INCLUDE_EXT_DIR=/home/feiw/libs/rdkit_2017_09_3/External -DRDKIT_LIBRARY_DIR=/home/feiw/libs/rdkit_2017_09_3/lib -DCMAKE_CXX_STANDARD=11 -DBOOST_ROOT=$BOOST_ROOT -DBoost_NO_BOOST_CMAKE=ON;
