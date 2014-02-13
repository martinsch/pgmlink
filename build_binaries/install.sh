ILASTIK_BUILD=$PWD
CPLEX_ROOT_DIR=$HOME/local

sudo apt-get install build-essential g++ gcc git cmake cmake-curses-gui gfortran libxext-dev &&

mkdir -p ilastik-build/build &&
cd ilastik-build &&

cmake -D BUILDEM_DIR=$ILASTIK_BUILD/build -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc -D CMAKE_CPLEX_ROOT_DIR=$CPLEX_ROOT_DIR $ILASTIK_BUILD &&

make &&

BUILDEM_DIR=$ILASTIK_BUILD/build 

cmake -D BUILDEM_DIR=$ILASTIK_BUILD/build $ILASTIK_BUILD &&

make -j8
