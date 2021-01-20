# Modified Routing Framework

This repository contains C++14 implementation for congestion-aware traffic assignment and autonmous mobility-on-demand routing. The code is based on the repository [Routing Framework](https://github.com/vbuchhold/routing-framework). The main additions that are provided in the current repository include the following:

* Support for computing rebalancing flows for autonomous mobility-on-demand routing, as described in Solovey, Salazar, and Pavone: “Scalable and Congestion-aware Routing for Autonomous Mobility-on-Demand via Frank-Wolfe Optimization”, in *Robotics: Science and Systems*, 2019.

* Batched OD-pairs that specify the number of passengers travelling between every origin-destination, rather than having exactly one passenger for each pair.

* Network graphs can be specified as CSV paths.

* Exogeneous flows can be specified to simulate exising traffic flow that is not controlled by the algorithm.

## Prerequisites

You need to have some tools and libraries installed. On Debian and its derivatives (such as Ubuntu)
the `apt-get` tool can be used:

```
$ sudo apt-get install build-essential
$ sudo apt-get install scons
$ sudo apt-get install python3
$ sudo apt-get install libboost-all-dev
$ sudo apt-get install libcairo2-dev
$ sudo apt-get install libnuma-dev
$ sudo apt-get install libproj-dev
$ sudo apt-get install zlib1g-dev
```

Next, you need to clone, build and install the libraries in the `External` subdirectory. To do so,
type the following commands at the top-level directory of the framework:

```
$ git submodule update --init
$ cd External
$ cd fast-cpp-csv-parser && sudo cp *.h /usr/local/include && cd ..
$ cd nanoflann && sudo cp -r include /usr/local && cd ..
$ cd RoutingKit && make && sudo cp -r include lib /usr/local && cd ..
$ cd vectorclass && sudo mkdir /usr/local/include/vectorclass && sudo cp *.h special/* $_ && cd ..
```

## Building

Once you installed the packages, simply type `scons` at the top-level directory of the framework to run in Devel mode:

```
$ scons
```
For Release or Debug modes use either `scons -Q variant-Release` or `scons -Q variant-Debug`.

### SCons Integration for Eclipse

The plugin SConsolidator provides tool integration for SCons in Eclipse.
Install it via its update site `http://www.sconsolidator.com/update`.

At the time of writing, there is a bug in the SConsolidator sources that you have to fix by hand.
If you are lucky, the bug will have been fixed in the official sources at the time you read this.
If not, try the following:

1. Find out where Eclipse installed the SConsolidator sources.
2. Open file `ch.hsr.ifs.sconsolidator.core/scons_files/BuildInfoCollector.py` with a text editor.
3. Change the definition of the function `get_compiler_flags` from

```
def get_compiler_flags(lang, environ):
    return (environ['CXXFLAGS'] if lang == 'c++' else environ['CCFLAGS'])
```

to

```
def get_compiler_flags(lang, environ):
    return [environ['CXXFLAGS'] if lang == 'c++' else environ['CCFLAGS']]
```

To import the framework as a SCons project into Eclipse, follow these steps:

1. Choose the menu `File`, `Import...`.
2. Select `C/C++`, `New SCons project from existing source`.
3. Choose the top-level directory of the framework.
4. Choose the menu `Project`, `Properties`, `C/C++ Build`, `Settings`, `Binary Parsers`.
5. Enable `Elf Parser`.
