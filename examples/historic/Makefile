VERSION = 1
SUBVERSION = 3
ANN_VERSION = 1.1.2
LP_SOLVE_VERSION = 5
LP_SOLVE_SUBVERSION = 5
ARCHITECTURE = ux64

TRACK = track2_${VERSION}.${SUBVERSION}
TRACK_CPLEX = track2_cplex_${VERSION}.${SUBVERSION}
TESTLIBS = test_libs_${VERSION}.${SUBVERSION}
TRACKOBJS = ${TRACK}.o
TRACK_CPLEXOBJS = ${TRACK_CPLEX}.o
TESTLIBSOBJS = test_libs_${VERSION}.${SUBVERSION}.o
SRC = src/track2_${VERSION}.${SUBVERSION}.cpp
SRC_TESTLIBS = src/test_libs_${VERSION}.${SUBVERSION}.cpp
FLAGS_STD = 
FLAGS_CPLEX = -DUSE_CPLEX

INCLUDES = -I${PWD}/../3rd_party/include

RPATH = ${PWD}/../3rd_party/lib

LIBPATH = -L${RPATH}

LIBS = -lhdf5 -lhdf5_hl -lANN 
LIBS_STD = -llpsolve${LP_SOLVE_VERSION}${LP_SOLVE_SUBVERSION}
LIBS_CPLEX = -lcplex  -lpthread

all : ${TRACK} ${TRACK_CPLEX} ${TESTLIBS}
	mv ${TRACK} bin
	mv ${TRACK_CPLEX} bin
	mv ${TESTLIBS} bin

${TRACK} : ${TRACKOBJS}
	${CXX} -o $@ ${TRACKOBJS} ${LIBPATH} ${LIBS} ${LIBS_STD} -Wl,-rpath,${RPATH}

${TRACK_CPLEX} : ${TRACK_CPLEXOBJS}
	${CXX} -o $@ ${TRACK_CPLEXOBJS} ${LIBPATH} ${LIBS} ${LIBS_CPLEX} -Wl,-rpath,${RPATH}

${TESTLIBS} : ${TESTLIBSOBJS}
	${CXX} -o $@ ${TESTLIBSOBJS} ${LIBPATH} ${LIBS} ${LIBS_STD} ${LIBS_CPLEX} -Wl,-rpath,${RPATH}

${TRACKOBJS} : ${SRC}
	${CXX} -c ${INCLUDES} ${FLAGS_STD} -o $@ $<

${TESTLIBSOBJS} : ${SRC_TESTLIBS}
	${CXX} -c ${INCLUDES} -o $@ $<

${TRACK_CPLEXOBJS} : ${SRC}
	${CXX} -c ${INCLUDES} ${FLAGS_CPLEX} -o $@ $<

clean :
	! (test -e ${TRACKOBJS}) || rm ${TRACKOBJS} 
	! (test -e ${TRACK_CPLEXOBJS}) || rm ${TRACK_CPLEXOBJS} 
	! (test -e ${TESTLIBSOBJS}) || rm ${TESTLIBSOBJS} 
	! (test -e bin/${TRACK}) || rm bin/${TRACK}
	! (test -e bin/${TRACK_CPLEX}) || rm bin/${TRACK_CPLEX}
	! (test -e bin/${TESTLIBS}) || rm bin/${TESTLIBS}
