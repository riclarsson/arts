[![Build](https://github.com/atmtools/arts/workflows/Build/badge.svg?branch=master)](https://github.com/atmtools/arts/commits/master)


Welcome to ARTS
===============

ARTS is free software. Please see the file COPYING for details.

If you use data generated by ARTS in a scientific publication, then please
mention this and cite the most appropriate of the ARTS publications that are
summarized on http://www.radiativetransfer.org/docs/

[CONTRIBUTING.md](CONTRIBUTING.md) provides information on contributing
to ARTS on GitHub.
 
For documentation, please see the files in the doc subdirectory.

For building and installation instructions please read below.


Building ARTS
-------------

Build Prerequisites:

- gcc/g++ >=8 (or llvm/clang >=8) older versions might work, but are untested
- cmake (>=3.1.0)
- zlib
- openblas
- netcdf (optional)
- Python3 (>=3.6)
  - required modules:
    docutils
    lark-parser
    matplotlib
    netCDF4
    numpy
    pytest
    scipy
    setuptools

To build the documentation you also need:

- pdflatex (optional)
- doxygen (optional)
- graphviz (optional)


Using cmake
-----------

Here are the steps to use cmake to build ARTS:

```
cd arts
mkdir build
cd build
cmake ..
make
```

If you only want to build the arts executable you can just run 'make arts'
instead of 'make'.

If you have a multi-core processor or multiprocessor machine, don't forget
to use the -j option to speed up the compilation:

```
make -jX
```

Where X is the number of parallel build processes.
X=(Number of Cores)+1 gives you usually the fastest compilation time.

Developer install of the PyARTS Python package:

```
cd python
pip install --user -e .
```

You only have to do the package install once. If the ARTS source has changed,
update the PyARTS package by running:

```
make -jX pyarts
```

Build configurations
--------------------
To build a release version without assertions or debugging symbols use:

```
cmake -DCMAKE_BUILD_TYPE=Release ..
make clean
make
```

To switch back to the debug version use:

```
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make clean
make
```

This is also the default configuration if you run cmake without
options in an empty build directory.

Native build (EXPERIMENTAL)
---------------------------

Finally, an experimental build type for Native infrastructures exists on GCC:

```
cmake -DCMAKE_BUILD_TYPE=Native ..
make clean
make
```

This option should make the executable slightly faster, more so on better
systems, but not portable. Note that since this build-mode is meant for
fast-but-accurate computations, some IEEE rules will be ignored. For now only
complex computations are IEEE incompatible running this mode of build.

Tests
-----

'make check' will run several test cases to ensure that ARTS is working
properly. Use 'make check-all' to run all available controlfiles, including
computation time-intensive ones.

Some tests depend on the arts-xml-data package. cmake automatically looks if it
is available in the same location as ARTS itself. If necessary, a custom path
can be specified:

```
cmake -DARTS_XML_DATA_PATH=/home/myname/arts-xml-data ..
```

If arts-xml-data cannot be found, those tests are ignored.

By default, the tests are executed serially.
If you want to run them concurrently, you can use:

```
cmake -DTEST_JOBS=X ..
```

X is the number of tests that should be started in parallel.

You can also use the ctest command directly to run the tests:

```
ctest -j4
```

To run specific tests, use the -R option and specify part of the test case name
you want to run. The following command will run all tests that have 'ppath' in
their name, e.g. arts.ctlfile.fast.ppath1d ...:

```
ctest -R ppath
```

To see the output of ARTS, use the -V option:

```
ctest -V -R fast.doit
```

By default, ctest will not print any output from ARTS to the screen. The option
--output-on-failure can be passed to ctest to see output in the case an error
occurs. If you want to always enable this, you can set the environment variable
CTEST_OUTPUT_ON_FAILURE:

```
export CTEST_OUTPUT_ON_FAILURE=1
```


HITRAN catalog support
----------------------

By default, ARTS only supports the latest HITRAN 2012 catalog version. Because
isotopologues have been renamed between different catalog versions, ARTS needs
to be compiled for one specific HITRAN version. If you want to use HITRAN 2008,
you have to recompile ARTS with:

```
cmake -DWITH_HITRAN2008=1 ..
make arts
```

To switch back to HITRAN 2012, run:

```
cmake -DWITH_HITRAN2008=0 ..
make arts
```


Optional features
-----------------

To include features that rely on Fortran code located in the 3rdparty
subdirectory use:

```
cmake -DENABLE_FORTRAN=1 -DCMAKE_Fortran_COMPILER=gfortran ..
```

This enables Disort, Fastem, Refice and Tmatrix.

If necessary, certain Fortran modules can be selectively disabled:

```
cmake -DENABLE_FORTRAN=1 -DNO_DISORT=1 ..
cmake -DENABLE_FORTRAN=1 -DNO_REFICE=1 ..
```

IMPORTANT: Only gfortran and Intel Fortran are currently supported.
Also, a 64-bit system is required (size of long type must be 8 bytes).


Enable NetCDF: The basic matpack types can be read from NetCDF files, if NetCDF
support is enabled:

```
cmake -DENABLE_NETCDF=1 ..
```


Disabling features
------------------

By default, a library to use the ARTS workspace interface in typhon is built.
You can disable building the C API:

```
cmake -DNO_C_API=1 ..
```

Disable assertions:
```
cmake -DNO_ASSERT=1 ..
```

Disable OpenMP:
```
cmake -DNO_OPENMP=1 ..
```

Disable the built-in documentation server:
```
cmake -DNO_DOCSERVER=1 ..
```

Treat warnings as errors:
```
cmake -DWERROR=1 ..
```

Disable FFTW autodetection:

ARTS automatically detects the availability of the FFTW3 library needed to
speed up the calculation of HITRAN cross section species . If you need to
disable FFTW support for any reason, you can do so with the following cmake
option:

```
cmake -DNO_FFTW=1 ..
```


TMatrix precision
-----------------

By default, ARTS uses double-precision for the T-matrix calculations.
When using the Intel compiler, quad-precision can be enable with cmake:

```
cmake -DENABLE_FORTRAN=1 -DENABLE_TMATRIX_QUAD=1 ..
```

Note that quad-precision is software emulated. T-matrix calculations will
around 10x slower.


ccache support
--------------

To utilize ccache when available use:

```
cmake -DENABLE_CCACHE=1 ..
```

For details see https://ccache.samba.org/


Intel compiler
--------------

If you want to compile with the Intel compiler[1], start with an empty build
directory and run:

```
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc ..
```

[1] http://software.intel.com/c-compilers


LLVM/Clang compiler
-------------------

If you want to compile with the LLVM/Clang compiler[1], start with an empty
build directory and run:

```
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..
```

You might also have to explicitly pick the right Fortran compiler since clang
doesn't have one:

```
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
-DENABLE_FORTRAN=1 -DCMAKE_Fortran_COMPILER=gfortran ..
```

Note that at this point, on OS X the default Apple Clang compiler does not
support OpenMP. Other versions of Clang support it via libomp.

[1] http://clang.llvm.org
[2] http://libcxx.llvm.org


macOS / Xcode
-----------

If you're on a Mac and have the Apple Xcode development environment installed,
you can generate a project file and use Xcode to build ARTS:

```
cmake -G Xcode ..
open ARTS.xcodeproj
```


Experimental features (ONLY USE IF YOU KNOW WHAT YOU'RE DOING)
--------------------------------------------------------------

Enable C++17 (only for compatibility testing, do not use C++17 features in your
code):

```
cmake -DENABLE_CXX17=1 ..
```


Valgrind profiling
------------------

The callgrind plugin included in valgrind is the recommended profiling method
for ARTS.

Due to limitations of valgrind, you need to disable the tmatrix code
(-DNO_TMATRIX=1) when compiling ARTS with Fortran support.

Certain things should be taken into account when calling ARTS with valgrind.
Since recursion (cycles) will lead to wrong profiling results it is
important to use the following settings to obtain profile data for ARTS:

```
valgrind --tool=callgrind --separate-callers=10 --separate-recs=3 arts -n1 ...
```

For detail on these options consult the valgrind manual:

http://valgrind.org/docs/manual/cl-manual.html#cl-manual.cycles

-n1 should be passed to ARTS because parallelisation can further scew the
results. Since executing a program in valgrind can lead to 50x slower
execution, it is recommended to create a dedicated, minimal controlfile for
profiling.

After execution with valgrind, the resulting callgrind.out.* file can be
opened in kcachegrind[1] for visualization. It is available as a package for
most Linux distributions.

Note that you don't have to do a full ARTS run. You can cancel the program
after some time when you think you have gathered enough statistics.

[1] https://kcachegrind.github.io/


Linux perf profiling
--------------------

The [Performance Counters for Linux](https://perf.wiki.kernel.org/) offer a
convenient way to profile any program with basically no runtime overhead.
Profiling works for all configurations (Debug, RelWithDebInfo and Release). To
ensure that the calltree can be analyzed correctly, compile ARTS without frame
pointers. This has minimal impact on performance:

```
cmake -DCMAKE_C_FLAGS="-fno-omit-frame-pointer" \
      -DCMAKE_CXX_FLAGS="-fno-omit-frame-pointer" ..
```

Prepend the perf command to your arts call to record callgraph information:

```
perf record -g src/arts MYCONTROLFILE.arts
```

This can also be applied to any test case:

```
perf record -g ctest -R TestDOIT$
```

After recording, use the report command to display an interactive view of the
profiling information:

```
perf report -g graph,0.5,callees
```

This will show a reverse call tree with the percentage of time spent in each
function. The function tree can be expanded to expose the calling functions.

