# Frank-Wolfe for Traffic Assignment

This repository contains C++14 implementation for congestion-aware traffic assignment (TA) and autonmous mobility-on-demand routing (AMoD) using the Frank-Wolfe algorithm. The code is based on the repository [Routing Framework](https://github.com/vbuchhold/routing-framework). The main additions that are provided in the current repository include the following:

* Computation of rebalancing flows for autonomous mobility-on-demand routing (AMoD), as described in Solovey, Salazar, and Pavone: “Scalable and Congestion-aware Routing for Autonomous Mobility-on-Demand via Frank-Wolfe Optimization”, in *Robotics: Science and Systems*, 2019.

* Computation of Interpolated Traffic Assignment (I-TAP), which balance fairness and efficiencity, as described in Jalota, Solovey, Zoepf, and Pavone: “Who Gets What and Why? A Fast Method to Balance Efficiency and Fairness in Traffic Routing”, *ArXiv*, 2021. Additional code can be found in [git:fair-routing](https://github.com/StanfordASL/fair-routing).

* Computation of Constrained System Optimum (CSO), which balance fairness and efficiencity, as described in Jahn, Möhring, Schulz, and Stier-Moses: “System-Optimal Routing of Traffic Flows with User Constraints in Networks with Congestion“, *Operations Research*, 2005.

* Path-based solution representations.

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
For Release or Debug modes use either `scons -Q variant=Release` or `scons -Q variant=Debug`.

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

## Useage
The main program for computing TA and AMoD is provided in `AssignTraffic`. Basic useage is specified in the `help` option:

```
Build/Devel/Launchers/./AssignTraffic -help
```

In addition to the running parameters, the main inputs consist of a graph that is given in a binary representation and a CSV file describing the OD pairs. The program `ConvertGraph` can be used to convert graphs represented in other formats (including CSV) into a binary format.

### Input representation
* A CSV graph is defined via a vertex and edge CSV files, which can be then converted into a binary format. The vertex file specifies vertex IDs and their corresponding xy coordinates. E.g.,

|vert_id|xcoord|ycoord|
|-------|------|------|
|0|-74.0168143|40.7051367|
|1|-74.0164164|40.7048817|
|2|-74.0164046|40.7047995|
|...|...|...|

* The edge file specifies the edges (with respect to the vertex IDs) and their attributes, where edge length is given in meters, capacity specifies the number of passing vehicles in an hour, and speed denotes the free-flow speed (km/h).  E.g.,

|edge_tail|edge_head|length|capacity|speed|
|---------|---------|------|--------|-----|
|3|0|172|14783|56|
|2|1|9|14783|56|
|8|1|121|5279|40|
|0|2|51|10559|40|
|...|...|...|...|...|

* The od-pairs file specifies origin and destination vertices for a given demand, and the volume (number of vehicles or people). E.g.,

|origin|destination|volume|
|------|-----------|------|
|20|31|96|
|20|41|120|
|20|57|72|
|20|65|88|
|...|...|...|

### Output representation
Below we describe the various outputs that the algorithm can produce according to the flags given to `AssignTraffic'.

* `-o`: Summary output file that specifies running time and solution quality according with respect to the algorithm iteration.

* `-fp`: Flow pattern obtained after each algorithm iteration, which specifies the total number of vehicles going through any given edge in the graph.

* `-paths`: Specifies the path computed for every OD pair for a given iteration of the algorithm.

* `-w`: Specifies the weights of the aforementioned paths in the final solution. 
