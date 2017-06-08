In case we cannot access root account...

# Environment Setting: add the following to the end of ~/.bashrc, run source ~/.bashrc
```
export PATH=/home/hli2/libs/gcc/bin:/home/hli2/libs/gawk/bin:/home/hli2/libs/binutils/bin:/home/hli2/libs/glibc/bin:/home/hli2/libs/mpi/bin:$PATH
export LD_LIBRARY_PATH=/home/hli2/libs/binutils/lib:/home/hli2/libs/gcc/lib:/home/hli2/libs/gcc/lib64:/home/hli2/libs/gawk/lib:/home/hli2/libs/gmp/lib:/home/hli2/libs/mpfr/lib:/home/hli2/libs/mpc/lib:/home/hli2/libs/glibc/lib
```

#install m4
```
wget ftp://ftp.gnu.org/gnu/m4/m4-1.4.17.tar.gz
tar -xvf m4-1.4.17.tar.gz
cd m4-1.4.17/
./configure --prefix=$HOME/libs/m4
make
make install
cd ..
rm m4-1.4.17 -rf
rm m4-1.4.17.tar.bz2
```

#install gmp
```
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/gmp-6.1.0.tar.bz2
tar -xvf gmp-6.1.0.tar.bz2
cd gmp-6.1.0
./configure --prefix=$HOME/libs/gmp M4=$HOME/libs/m4/bin/m4
make
make install
cd ..
rm -rf gmp-6.1.0
rm gmp-6.1.0.tar.bz2
```

#install mpfr
```
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpfr-3.1.4.tar.bz2
tar -xvf mpfr-3.1.4.tar.bz2 
cd mpfr-3.1.4/
../configure --prefix=$HOME/libs/mpfr with_gmp=$HOME/libs/gmp 
make
make install
cd ..
rm mpfr-3.1.4 -rf
rm mpfr-3.1.4.tar.bz2
```

#install mpc
```
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/mpc-1.0.3.tar.gz
tar -xvf mpc-1.0.3.tar.gz 
cd mpc-1.0.3/
./configure --prefix=$HOME/libs/mpc --with-mpfr=$HOME/libs/mpfr with_gmp=$HOME/libs/gmp 
make
make install
cd ..
rm -rf mpc-1.0.3
rm mpc-1.0.3.tar.gz
```
#install GCC with gcc and gfortran
```
export LD_LIBRARY_PATH=/home/hli2/libs/gmp/lib:/home/hli2/libs/mpfr/lib:/home/hli2/libs/mpc/lib
wget ftp://gd.tuwien.ac.at/gnu/gcc/releases/gcc-6.1.0/gcc-6.1.0.tar.bz2
tar -xvf gcc-6.1.0.tar.bz2
cd gcc-6.1.0/
./configure --prefix=$HOME/libs/gcc --enable-languages=c,c++,fortran --disable-multilib --enable-checking=release --disable-bootstrap --with-gmp=$HOME/libs/gmp  --with-mpfr=$HOME/libs/mpfr --with-mpc=$HOME/libs/mpc
make
make install
cd ..
rm gcc-6.1.0 -rf
rm gcc-6.1.0.tar.bz2 
```

#install OpenBLAS
```
wget http://github.com/xianyi/OpenBLAS/archive/v0.2.18.tar.gz
tar -xvf v0.2.18.tar.gz 
cd OpenBLAS-0.2.18/
```
Edit Makefile.rule and change gcc and fc to above gcc folder.
```
make
make PREFIX=$HOME/libs/OpenBLAS install
cd ..
rm OpenBLAS-0.2.18/ -rf
rm v0.2.18.tar.gz
```

#install LAPACK
```
wget http://www.netlib.org/lapack/lapack-3.6.0.tgz
tar -xvf lapack-3.6.0.tgz 
cd lapack-3.6.0
```
Edit make.inc.Example and rename it to make.inc 
Change FORTRAN, LOADER, cc
Change MakeFile:
```
lib: lapacklib tmglib
#lib: blaslib variants lapacklib tmglib
```
to:
```
#lib: lapacklib tmglib
lib: blaslib variants lapacklib tmglib
```
Then type:
```
make
```

#install Armadillo
```
wget http://sourceforge.net/projects/arma/files/armadillo-7.200.2.tar.xz
tar -xvf armadillo-7.200.2.tar.xz
cd armadillo-7.200.2
```
Edit CMakeLists.txt and add:
```
set(OpenBLAS_LIBRARY "/home/hli2/libs/OpenBLAS/lib/libopenblas.a")
set(OpenBLAS_INCLUDE_DIR "/home/hli2/libs/OpenBLAS/include")
set(LAPACK_LIBRARY "/home/hli2/libs/lapack-3.6.0/liblapack.a")
set(LAPACK_LIBRARY "/home/hli2/libs/lapack-3.6.0/librefblas.a")
set(LAPACK_LIBRARY "/home/hli2/libs/lapack-3.6.0/libtmglib.a")
```
at the beginning and add:
```
set(BUILD_SHARED_LIBS OFF)
```
in the middle before the message. Then：
```
mkdir build
cd build
cmake ../
make
```
Copy include in the main  folder and build/libarmadillo.a to the project folder and add them to CMakeList.txt.
Or:
```
make
make install DESTDIR=$HOME/libs/armadillo
cd ..
rm armadillo-7.200.2 -rf
rm armadillo-7.200.2.tar.xz 
```

#install Boost
```
sudo apt-get update
sudo apt-get install libboost-all-dev
```

```
tar -xvf boost_1_61_0.tar.bz2
cd boost_1_61_0
echo "using gcc : 6.1 : /home/hli2/libs/gcc/bin/g++ ; " >> tools/build/src/user-config.jam
./bootstrap.sh --prefix="/home/hli2/libs/boost"
./b2 install
cd ..
rm boost_1_61_0 -rf
rm boost_1_61_0.tar.bz2 
```

#install mlpack
```
wget http://mlpack.org/files/mlpack-2.0.2.tar.gz
tar -xvf mlpack-2.0.2.tar.gz 
cd mlpack-2.0.2/
```
Edit CMakeLists.txt and add followings to the top:
```
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lgfortran")
set(ARMADILLO_LIBRARY "/home/hli2/libs/armadillo/usr/lib/libarmadillo.so")
set(ARMADILLO_INCLUDE_DIR "/home/hli2/libs/armadillo/usr/include")

set(BOOST_ROOT "/home/hli2/libs/boost")
set(BOOST_INCLUDEDIR "/home/hli2/libs/boost/include")
set(BOOST_LIBRARYDIR "/home/hli2/libs/boost/lib")
set(Boost_NO_SYSTEM_PATHS TRUE)
```
Then:
```
mkdir build
cd build
export LD_LIBRARY_PATH=/home/hli2/libs/gmp/lib:/home/hli2/libs/mpfr/lib:/home/hli2/libs/mpc/lib
cmake -DCMAKE_C_COMPILER='/home/hli2/libs/gcc/bin/gcc' -DCMAKE_CXX_COMPILER='/home/hli2/libs/gcc/bin/g++' -DCMAKE_EXE_LINKER_FLAGS='-Wl,-rpath,/home/hli2/libs/gcc/lib64' ../
make install DESTDIR=$HOME/libs/mlpack
```
Remove const from the error file as in the warming.

#install libxml2
```
wget ftp://xmlsoft.org/libxml2/libxml2-2.9.4.tar.gz
tar xvf libxml2-2.9.4.tar.gz
cd libxml2-2.9.4
./configure --prefix=$HOME/libs/libxml2
make
make install
cd ..
rm -rf libxml2-2.9.4
rm libxml2-2.9.4.tar.gz
```

#install GSL
```
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.1.tar.gz
tar xvf gsl-2.1.tar.gz
cd gsl-2.1
./configure --prefix=$HOME/libs/gsl
make
make install
rm -rf gsl-2.1
rm gsl-2.1.tar.gz
```

#install GKlib
```
mkdir tmp
tar xvf GKlib.tgz -C /tmp
cd /tmp/GKlib
make config prefix=$HOME/libs/GKlib
make
make install
cd ../../
rm -rf tmp
rm GKlib.tgz
```

#install binutils
```
wget http://ftp.gnu.org/gnu/binutils/binutils-2.26.tar.gz
cd binutils-2.26
CC=/home/petrie/Workspace/hli/libs/gcc/bin/gcc ./configure --prefix=/home/petrie/Workspace/hli/libs/binutil
make
make install
```

#Install tbb

This cannot support using with MPI:
```
sudo apt-get update
sudo apt-get install libtbb-dev
```

Or
Download the version for linux: 
https://www.threadingbuildingblocks.org/sites/default/files/software_releases/linux/tbb2017_20161128oss_lin_0.tgz
Extract files to /home/hli/Libs/tbb, edit /tbb/bin/tbbvars.sh, set:
```
TBBROOT=/home/hli/Libs/tbb
```
then run:
```
. /home/hli/Libs/tbb/bin/tbbvars.sh intel64
```
Add the above to .bashrc.

#Install autoconf
```
    cd /home/hli2/libs
    mkdir autoconf
    wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz
    tar -xvf autoconf-2.69.tar.gz
    cd autoconf-2.69
    ./configure CC=/home/hli2/libs/gcc/bin/gcc CCX=/home/hli2/libs/gcc/bin/g++ M4=/home/hli2/libs/m4/bin/m4 --prefix=/home/hli2/libs/autoconf
    make 
    make install
    cd ..
    rm -rf autoconf-2.69
    rm autoconf-2.69.tar.gz
```

#Install automake
```
    cd /home/hli2/libs
    mkdir automake
    wget http://ftp.gnu.org/gnu/automake/automake-1.15.tar.gz
    tar -xvf automake-1.15.tar.gz
    cd automake-1.15
    ./configure CC=/home/hli2/libs/gcc/bin/gcc CCX=/home/hli2/libs/gcc/bin/g++ AUTOCONF=/home/hli2/libs/autoconf/bin/autoconf prefix=/home/hli2/libs/automake
    make 
    make install 
    cd ..
    rm -rf automake-1.15
    rm automake-1.15.tar.gz
```

#Install MPI
MPICH3 cannot be installed using apt.


root:
```
wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
tar -xvf mpich-3.2.tar.gz
cd mpich-3.2/
./configure 2>&1 | tee c.txt
make 2>&1 | tee m.txt
sudo make install |& tee mi.txt
cd ../
rm mpich-3.2/ -rf
```

galaxy:
```
cd /home/hli2/libs
mkdir /mpi
wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
tar -xvf mpich-3.2.tar.gz
cd mpich-3.2/
./configure --disable-fortran CC=/home/hli2/libs/gcc/bin/gcc CCX=/home/hli2/libs/gcc/bin/g++ AUTOMAKE=/home/hli2/libs/automake/bin/automake --disable-fortran --prefix=/home/hli2/libs/mpi 2>&1 | tee c.txt
make 2>&1 | tee m.txt
make install 2>&1 | tee mi.txt
cd ../
rm -rf mpich-3.2
rm mpich-3.2.tar.gz
```

#Install libaio
```
cd /home/hli2/libs
wget https://git.fedorahosted.org/cgit/libaio.git/snapshot/libaio-0.3.109.tar.gz
tar -xvf libaio-0.3.109.tar.gz
cd libaio-0.3.109
make
make install prefix=/home/hli2/libs/libaio
cd ../
rm -rf libaio-0.3.109
rm libaio-0.3.109.tar.gz
```

#Install gawk
```
cd /home/hli2/libs
wget https://ftp.gnu.org/gnu/gawk/gawk-4.1.4.tar.gz
tar -xvf gawk-4.1.4.tar.gz
cd gawk-4.1.4
./configure prefix=/home/hli2/libs/gawk
make 
make install
cd ../
rm gawk-4.1.4 -rf
rm gawk-4.1.4.tar.gz
```

#Install binutils
```
cd /home/hli2/libs
wget https://ftp.gnu.org/gnu/binutils/binutils-2.27.tar.gz
tar -xvf binutils-2.27.tar.gz
cd binutils-2.27
./configure prefix=/home/hli2/libs/binutils
make
make install
cd ../

rm binutils-2.27 -rf
rm binutils-2.27.tar.gz
```

#Install glibc
```
cd /home/hli2/libs
wget https://ftp.gnu.org/gnu/libc/glibc-2.24.tar.gz
tar -xvf glibc-2.24.tar.gz
cd glibc-2.24
mkdir build
cd build
../configure prefix=/home/hli2/libs/glibc
make
make install
cd ../../
rm glibc-2.24.tar.gz
rm glibc-2.24 -rf
```