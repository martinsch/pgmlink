# pgmLink - A Tracking-by-Assignment Software
pgmLink is about *tracking by assignment with probabilistic graphical models and other approaches*. 
It is written in C++ (core) and Python (high level, optional). 

A graphical user interace (GUI) for our methods is freely available on github [martinsch/ilastik](http://github.com/martinsch/ilastik)
and binaries for Windows/Mac/Linux are provided upon request.


pgmlink provides an implementation of

1. **Conservation tracking**, as it is described in

   M. Schiegg, P. Hanslovsky, B. X. Kausler, L. Hufnagel, F. A. Hamprecht. 
   [**Conservation Tracking**](http://hci.iwr.uni-heidelberg.de/Staff/mschiegg/schiegg_13_conservation.pdf). Proceedings of the IEEE International Conference 
   on Computer Vision (ICCV 2013), 2013.

2. **Chaingraph tracking**, as it is described in

   B. X. Kausler, M. Schiegg, B. Andres, M. Lindner, H. Leitte, L. Hufnagel, 
   U. Koethe, F. A. Hamprecht. [**A Discrete Chain Graph Model for 3d+t Cell 
   Tracking with High Misdetection Robustness**](http://hci.iwr.uni-heidelberg.de//Staff/bkausler/kausler_12_discrete.pdf). Proceedings of the European 
   Conference on Computer Vision (ECCV 2012), 2012.




Please cite the appropriate paper if you use this software.
Feel free to contact us if you have any questions or suggestions.


## Dependencies
Dependencies that have to be built and installed manually:

- LEMON http://lemon.cs.elte.hu/trac/lemon
- VIGRA https://github.com/ukoethe/vigra
- openGM https://github.com/opengm/opengm
- cplex/concert http://www-01.ibm.com/software/integration/optimization/cplex-optimizer/
- dlib http://dlib.net/
- mlpack (>= 1.0.8) http://mlpack.org/

Dependencies that should be available as packages:

- ann
- boost
  - boost-serialization
  - boost-random
  - boost-program-options
  - boost-python (optional)
  - boost-test (optional)
- armadillo http://arma.sourceforge.net/
- libxml2-dev
- doxygen (optional)

### CPLEX
pgmLink depends on [IBM ILOG CPLEX](http://www-01.ibm.com/software/integration/optimization/cplex-optimization-studio/), in particular `libcplex`, `libilocplex` and `libconcert`. If you are an academic you can obtain a free license from the [IBM Academic Initiative](http://www-03.ibm.com/ibm/university/academic/pub/page/academic_initiative).

#### Installation fails with an *internal LaunchAnywhere application error*
Some people encounter the following error during the CPLEX installation:

```
Launching installer...

An internal LaunchAnywhere application error has occured and this application cannot proceed. (LAX)

Stack Trace:
java.lang.IllegalArgumentException: Invalid Unicode sequence: illegal character
    at java.util.Properties.loadImpl(Properties.java:356)
    at java.util.Properties.load(Properties.java:288)
    at com.zerog.common.java.util.PropertiesUtil.loadProperties(DashoA10*..)
	at com.zerog.lax.LAX.<init>(DashoA10*..)
	at com.zerog.lax.LAX.main(DashoA10*..)

```
This installer bug happens when the installer is called from a shell with a modified appearance (in particular, a modified prompt via PS1). **Install CPLEX from a unaltered shell to circumvent the bug**.

#### CPLEX shared libraries


The CPLEX package does not provide shared versions of all required libraries, but only static variants (luckily with PIC). You can link your own shared libraries using the following commands

on Linux:
```
g++ -fpic -shared -Wl,-whole-archive libcplex.a -Wl,-no-whole-archive -o libcplex.so
g++ -fpic -shared -Wl,-whole-archive libilocplex.a -Wl,-no-whole-archive -o libilocplex.so
g++ -fpic -shared -Wl,-whole-archive libconcert.a -Wl,-no-whole-archive -o libconcert.so
```
on Linux and

```
g++ -fpic -shared -Wl,-all_load libcplex.a -Wl,-noall_load -o libcplex.dylib
g++ -fpic -shared -Wl,-all_load libilocplex.a -Wl,-noall_load -o libilocplex.dylib
g++ -fpic -shared -Wl,-all_load libconcert.a -Wl,-noall_load -o libconcert.dylib
```
on Mac respectively.

on Mac:
```
g++ -fpic -shared -Wl,-all_load libcplex.a -Wl,-noall_load -o libcplex.dylib
g++ -fpic -shared -Wl,-all_load libconcert.a -Wl,-noall_load -o libconcert.dylib
g++ -fpic -shared -Wl,-all_load libilocplex.a -Wl,-noall_load -L/PATH/TO/CONCERT/LIB -L/PATH/TO/CPLEX/LIB -lconcert -lcplex -o libilocplex.dylib
```


## How to Build
pgmLink is built with [cmake](www.cmake.org):
```
cmake <path-to-pgmlink>/.
make
```

There are several build options which you can set - among others - with `ccmake`:

- **CMAKE_BUILD_TYPE**: Set to *RELEASE* for maximal performance. *DEBUG* disables optimization, adds debugging symbols, and tolerates warnings (usually interpreted as errors).

- **LOGGING_LEVEL**: The logging level is set at compile time. You have no performance impact due to logging at all when setting it to *NO_LOGGING*. 

- **WITH_CHECKED_STL**: Use a safety-enhanced STL (GCC only). May impact performance.

- **WITH_PYTHON**: Add Python wrappers (requires boost-python).

- **WITH_TESTS**: Compile unit tests. You can execute them via `make test`


## Build instructions for armadillo
Before compilation, modify `include/armadillo_bits/config.hpp` to include
```
#define ARMA_USE_LAPACK   (needs liblapack)
#define ARMA_USE_BLAS     (needs libblas)
#define ARMA_USE_WRAPPER  (use wrapper for lapack and blas)
```
(uncomment if necessary).

## Build instructions for mlpack:
After installation of `libxml2-dev`, create a softlink to `libxml` in your `include` directory. If installed globally, the appropriate command is
```
sudo ln -s /usr/include/libxml2/libxml /usr/include/
```


# Build GUI

Pgmlink has been integrated into a [fork of ilastik](https://github.com/martinsch/ilastik). In particular, for both Chaingraph tracking and Conservation tracking, ilastik workflows have been implemented to ease use of these algorithms in a user-friendly software. For user documentation, please see the [official ilastik homepage](http://www.ilastik.org).

Binaries for an ilastik copy where Conservation tracking is integrated, are available on request. However, with the tools provided here, it is easy to automatically build ilastik-for-conservation-tracking yourself:

For this, it is enough to clone the [(modified) ilastik build repository](https://github.com/martinsch/ilastik-build-Linux.git) and follow the instructions in the readme there. It is essential to install CPLEX before building the binaries and to specify the CPLEX location as a build macro, as described in the said readme.
