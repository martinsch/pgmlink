ILASTIK_BUILD=$(dirname $0)/ilastik-build
CPLEX_ROOT_DIR=$HOME/local

sudo apt-get install build-essential g++ gcc git cmake cmake-curses-gui gfortran libxext-dev libfontconfig1-dev&&

mkdir -p $ILASTIK_BUILD/build &&
cp run_ilastik.sh $ILASTIK_BUILD/build &&

cd $ILASTIK_BUILD &&

BUILDEM_DIR=$ILASTIK_BUILD/build 

cmake -D BUILDEM_DIR=$BUILDEM_DIR -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc -D CMAKE_CPLEX_ROOT_DIR=$CPLEX_ROOT_DIR $ILASTIK_BUILD/../.. &&

make &&

cmake -D BUILDEM_DIR=$ILASTIK_BUILD/build $ILASTIK_BUILD/../.. &&

make -j
