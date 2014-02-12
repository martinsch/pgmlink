ILASTIK_BUILD=$PWD

sudo apt-get install build-essential g++ gcc git cmake cmake-curses-gui vim gfortran libxext-dev &&

mkdir -p ilastik-build/build &&
cd ilastik-build &&

cmake -D BUILDEM_DIR=$ILASTIK_BUILD/build -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc $ILASTIK_BUILD &&

make &&

BUILDEM_DIR=$ILASTIK_BUILD/build 

cmake -D BUILDEM_DIR=$ILASTIK_BUILD/build $ILASTIK_BUILD &&

make -j8
