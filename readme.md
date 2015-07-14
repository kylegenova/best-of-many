An Implementation of the Best-of-Many Christofides' Algorithm for the Traveling Salesman Problem
=======================================

Description
------------
This project contains the implementations of the algorithms described in [this paper][arXiv], which compile into a single executable called Best-of-Many. The procedures available are Christofides', Column Generation, Max Entropy Sampling, Splitting Off and Tree Packing, Column Generation + SwapRound, and Splitting Off and Tree Packing + SwapRound.  

The executable supports .tsp and .tsv formatted files for input, with a few exceptions (i.e. files supplying custom distance functions). 'Program file' inputs specifying a subset of the procedures on a set of files are also supported, with .csv spreadsheet output. Finally, the executable can be used in conjunction with the provided Python scripts to generate plots of the performance of the algorithms over time. Right now only a Windows version is available.


[arXiv]:http://arxiv.org/abs/1506.07776


Dependencies and Installation:
------------------------------

After checking out the code on this repository, install the following dependencies on your machine:
+ [Cygwin][cyg] 32bit (the minimum set of packages is fine; only the core is a depended upon). After installing, you must add `<cygwin directory>\bin` (which is `C:\cygwin\bin` by default) to either your user or system Path environment variable (with a semicolon separator from the rest of the path). This is required so that Concorde has access to Cygwin. After adding the reference to your environment variables and before using best-of-many, you must unfortunately restart your computer so that the change is propogated to system services.
+ [Visual Studio][vs] (Express 2013 or higher). This is used to compile Best-of-Many as well as the dependencies provided as source files.
+ [Gurobi Optimizer 5.6.3, 32-Bit][gur]. Gurobi is a fast LP solver used only in the Column Generation procedures. Install the program to the default `C:\gurobi563\` directory. Make sure to install a valid license before proceeding further (instructions are available on Gurobi's site). Free licenses for research purposes are available. If you can't get a license for some reason, you can install Gurobi anyway- the Column Generation functions won't work (and you'll have to follow the Usage section guidelines on how not to run it), but the rest of the program should compile and run without issue.

[cyg]: https://www.cygwin.com/
[vs]: https://www.visualstudio.com/en-us/products/free-developer-offers-vs.aspx
[gur]: http://www.gurobi.com/downloads/download-center

Next, download the following zips/archives:
+ The [Concorde TSP-Solver][cc] command-line executable. Select the 'windows-cygwin' version.
+ Kolmogorov's [blossomV][bV] code for computing minimum perfect matchings.
+ [Triangle][tri], a library for generating exact delaunay triangulations. It is used both by blossomV and by Best-of-Many.
+ [Eigen][eig], a Linear Algebra Library.

[cc]: http://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm
[bV]: http://pub.ist.ac.at/~vnk/software.html
[tri]: http://www.cs.cmu.edu/~quake/triangle.html
[eig]: http://eigen.tuxfamily.org/index.php?title=Main_Page

Fully unzip all four of the archives that were just downloaded. Name the 4 folders as follows: `eigen\`, `blossomV\`, `delaunay\` and `concorde.exe\` (Each extracted folder should contain item(s) other than just another folder in its immediate subdirectory). Place the four folders inside `installation\`.Run the [Python 2.x][py] script `installer.py`, which is also located in `installation\`. If you get any errors, just delete the directory `dependencies\`, which the script creates, and try again after making the necessary changes.

[py]: https://www.python.org/downloads/

Now, you may open the file `<repository root>\windows\best-of-many.sln` in Visual Studio and build it in Release or Debug mode, generating best-of-many.exe.  An important final requirement is that the directory in which Best-of-Many is located also contains `concorde.exe`; this executable is currently located in the directory `dependencies\` (which is created by `installer.py`).


Usage
-------
##### Basic Usage
The most basic way to run Best-of-Many is by passing it the paths to one or more TSPLIB files, each followed by their optimal cost (for the purpose of error calculation). It will run Christofides', Column Generation, Splitting Off, Maximum Entropy, and SwapRound on the trees from both Column Generation and Splitting Off, outputting the results as a CSV file. It will also save all the result and log files. Example:

`best-of-many afile.tsp 7 ../anotherfile.tsp 31 ../path/to/athirdfile.tsp 127`

##### Programs
For more complex operations, Best-of-Many runs "program" files, which are instruction sets for instances and algorithms to run on them. A program which reruns the experiments (excepting improvements to Best-of-Many since that date) from our paper is provided in `best-of-many\programs`. If you want this to program run, you must get the [TSPLIB][lib], [VSLI][circ], and [KONECT][unw] instances they reference, and extract them all to `best-of-many\instances`. Below is an example program file:
[lib]: http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/tsp/
[circ]:http://www.math.uwaterloo.ca/tsp/vlsi/
[unw]:http://konect.uni-koblenz.de/

```
# Welcome to the example program file for Best-of-Many.
# This is a comment because of the leading #.

# Whitespace is okay too.

# The first five lines determine if the corresponding algorithms will be run.
# SWAPROUND applies to both Column Generation and Splitting Off.
CHRISTOFIDES:TRUE
MAX_ENTROPY:TRUE
COLUMN_GENERATION:TRUE
SPLITTING_OFF:TRUE
SWAPROUND:TRUE
# Any value other than TRUE will be considered FALSE.
# When one of the above lines is ommitted, the variable defaults to TRUE.
# (so the five lines above are functionally irrelevant)

# This is the number of samples to do for each of the sample-based algorithms.
# The default is 100 if no count is provided. A malformed count might crash
# the program.
SAMPLE_COUNT:100

# These two constants determine how Column Generation cutoff is done.
# Column Generation halts when the total improvement over the last
# CG_HISTORY_SIZE iterations is less than CG_EPSILON.
# The values here are the defaults. A malformed input might crash the program.
CG_EPSILON:0.1
CG_HISTORY_SIZE:100

# These three arguments determine what is done with the output of the program.
# If SAVE_PLOTS is TRUE, some visualization files (in custom file-formats) will
# be generated. The files can be read be view_cg_plot.py and view_histogram.py.
# If SAVE_SPREADSHEET is TRUE, a csv spreadsheet with the results of the
# program will be generated. This is the main way to get data out of the program.
# If SAVE_LOG_FILES is TRUE, some log files generated by best-of-many and
# concorde will be saved.
SAVE_PLOTS:FALSE
SAVE_SPREADSHEET:TRUE
SAVE_LOG_FILES:FALSE

# This is the instances section. Every instance must be accompanied by the
# optimal tour cost.
# If the optimal isn't known the dependency concorde.exe can calculate it.
# Alternatively, just supply a value of 1; the error% statistics in the
# spreadsheet will be off, but the algorithms will run normally.

# The instances line is mandatory and must be exactly as it is below.
INSTANCES:

# Use Unix-style relative paths for the instance locations:
../instances/berlin52.tsp 7542
../instances/pr439.tsp 107217

# The paths are relative to the location of the program file.
# When running Best-of-Many, use a Unix-style relative path for the
# location of the program file as well.
```
