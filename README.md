# pgmLink - A Tracking-by-Assignment Software
pgmLink is about *tracking by assignment with probabilistic graphical models and other approaches*. 
It is written in C++ (core) and Python (high level, optional). 

## Dependencies
Dependencies that have to be built and installed manually:

- LEMON http://lemon.cs.elte.hu/trac/lemon
- VIGRA https://github.com/ukoethe/vigra
- openGM https://github.com/opengm/opengm
- cplex/concert http://www-01.ibm.com/software/integration/optimization/cplex-optimizer/

Dependencies that should be available as packages:

- ann
- boost
  - boost-serialization
  - boost-python (optional)
  - boost-test (optional)
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


The CPLEX package does not provide shared versions of all required libraries, but only static variants (luckily with PIC). You can link your own shared libraries using the following commands:

```
g++ -fpic -shared -Wl,-whole-archive libcplex.a -Wl,-no-whole-archive -o libcplex.so
g++ -fpic -shared -Wl,-whole-archive libilocplex.a -Wl,-no-whole-archive -o libilocplex.so
g++ -fpic -shared -Wl,-whole-archive libconcert.a -Wl,-no-whole-archive -o libconcert.so
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
