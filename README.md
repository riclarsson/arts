[![Build](https://github.com/atmtools/arts/workflows/Build/badge.svg?branch=arts3-dev)](https://github.com/atmtools/arts/commits/arts3-dev)


> [!CAUTION]
> This branch is under heavy development. Many features are missing or not fully implemented.
> If you are not an ARTS developer, you might want to stay away.

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


Dependencies
------------

Build Prerequisites (provided by Miniforge3):

- gcc/g++ >=12 (or llvm/clang >=16)
- cmake (>=3.18)
- zlib
- openblas
- libc++ (only for clang)
- libmicrohttpd (>=0.9, optional, for documentation server)
- netcdf (optional)
- Python3 (3.9, 3.10, 3.11, 3.12)
  - required modules:
    docutils
    lark-parser
    matplotlib
    netCDF4
    numpy
    pytest
    scipy
    setuptools
    xarray

To build the documentation you also need:

- pdflatex (optional)
- doxygen (optional)
- graphviz (optional)


Building ARTS
-------------

The following instructions assume that you are using Miniforge3 as a build environment.  The installer is available at
[the project's Github page](https://github.com/conda-forge/miniforge#miniforge).

Use the provided `environment-dev-{linux,mac}.yml` files to install
all required dependencies into your current conda environment.

Optionally, a separate environment for development can be created, if you want to keep your current environment clean:

```
mamba create -n pyarts-dev python=3.10
mamba activate pyarts-dev
```

Install dependencies on Linux:
```
mamba env update -f environment-dev-linux.yml
```

Install dependencies on macOS:
```
mamba env update -f environment-dev-mac.yml
```

Next, follow these steps to use `cmake` to build ARTS:
```
cd arts
cmake --preset=default-gcc-conda  # On macOS use default-clang-conda
cmake --build build -jX
```

X is the number of parallel build processes.
**X=Number of Cores** gives you usually the fastest compilation time.

WARNING: The compilation is very memory intensive. If you have 16GB of RAM,
don't use more than 6-8 cores. With 8GB, don't use more than 2-3 cores.

Development install of the PyARTS Python package:

```
python3 -m pip install --user -e build/python
```

You only have to do the python package install once.
If the ARTS source has changed, update the PyARTS package by running:

```
cmake --build build -jX --target pyarts
```


Build configurations
--------------------

By default, ARTS is built in release mode with optimizations enabled and
assertions and debugging symbols turned off.

Whenever you change the configuration, remove your build directory first:

```
rm -rf build
```

To build with assertions and debugging symbols use:

```
cmake --preset=default-gcc-conda -DCMAKE_BUILD_TYPE=RelWithDebInfo
```

This configuration offers a good balance between performance and debugging
capabilities. Since this still optimizes out many variables, it can be
necessary for some debugging cases to turn off all optimizations. For those
cases, the full debug configuration can be enabled. Note that ARTS runs a lot
slower in this configuration:

```
cmake --preset=default-gcc-conda -DCMAKE_BUILD_TYPE=Debug
```


Installing PyARTS
-----------------

To install the PyARTS Python package, you need to build it and install it with
pip. Create your build directory and configure ARTS with cmake as described in
the previous sections. Then, run the following commands inside your build
directory:

```
cmake --build build --target pyarts
python3 -m pip install --user -e build/python
```

This will not mess with your system's Python installation.
A link to the pyarts package is created in your home directory, usually
`$HOME/.local/lib/python3.X/site-packages/pyarts.egg-link`.

You don't need to reinstall the package with pip after updating ARTS.
You only need to run `cmake --build build --target pyarts` again.


Tests
-----

'cmake --build build --target check' will run several test cases to ensure that
ARTS is working properly. Use 'check-all' to run all available controlfiles,
including computation time-intensive ones.

Some tests depend on the arts-xml-data package. cmake automatically looks if it
is available in the same location as ARTS itself. If necessary, a custom path
can be specified.

```
cmake --preset=default-gcc-conda -DARTS_XML_DATA_PATH=/home/myname/arts-xml-data
```

If arts-xml-data cannot be found, those tests are ignored.

By default, 4 tests are executed in parallel.
If you change the number of concurrently run test, you can add this option to your `cmake --preset=....` call:

```
-DTEST_JOBS=X
```

X is the number of tests that should be started in parallel.

You can also use the ctest command directly to run the tests:

First, change to the `build` directory:
```
cd build
```

This runs all test with 4 jobs concurrently:
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


Native build
------------

To squeeze out every last drop of performance, you can also build a version
specifically optimized for your machine's processor:

```
-DCMAKE_BUILD_TYPE=Native
```

This option should make the executable slightly faster, more so on better
systems, but not portable. Note that since this build-mode is meant for
fast-but-accurate computations, some IEEE rules will be ignored. For now only
complex computations are IEEE incompatible running this mode of build.


Optional features
-----------------

Features that rely on Fortran code located in the 3rdparty
subdirectory are enabled by default, but can be disabled by passing the
following option to the `cmake --preset=...` command:

```
-DENABLE_FORTRAN=0
```

This disables Disort, Fastem and Tmatrix.

If necessary, certain Fortran modules can be selectively disabled:

```
-DNO_DISORT=1
```
or
```
-DENABLE_FORTRAN=1 -DNO_TMATRIX=1
```

IMPORTANT: Only gfortran is currently supported.
Also, a 64-bit system is required (size of long type must be 8 bytes).


Enable NetCDF: The basic matpack types can be read from NetCDF files, if NetCDF
support is enabled:

```
cmake --preset=default-gcc-conda -DENABLE_NETCDF=1
```

Precompiled headers: PCH can speed up builds significantly. However, it hampers
the ability for ccache to properly skip unnecessary compilations, potentially
increasing rebuild times. Tests have shown that it only speeds up the build
considerably for Clang, but not for GCC.

```
cmake --preset=default-clang-conda -DENABLE_PCH=1 ..
```

If you enable PCH and also use ccache, you need to set the `CCACHE_SLOPPINESS`
environment variable properly:

```
export CCACHE_SLOPPINESS=pch_defines,time_macros
```


Disabling features
------------------

Disable assertions: `-DNO_ASSERT=1`

Disable OpenMP: `-DNO_OPENMP=1`

Disable the built-in documentation server: `-DNO_DOCSERVER=1`


ccache support
--------------

The build utilizes ccache automatically when available, it can be
turned off with the option `-DENABLE_CCACHE=0`

For details see https://ccache.samba.org/


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
ensure that the calltree can be analyzed correctly, compile ARTS with frame
pointers. This has minimal impact on performance. Use the following preset to
enable this setting:

```
cmake --preset=perf-gcc-conda
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

