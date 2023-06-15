# FindDominationNumber 
This repository contains a program created for the article "J. Renders, S. Tokunaga, C.T. Zamfirescu: Dominating maximal outerplane graphs, manuscript".

The program uses McKay's graph6 format to read and write graphs. See <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>. As well as Brinkmann and Mckay's planar code format. The definition can be found in the `plantri` manual. See <https://users.cecs.anu.edu.au/~bdm/plantri/plantri-guide.txt>.

### Short manual
This program can be used to determine the domination number of a given input graph. One can also count the number of minimal dominating sets and output graphs with a specific domination number. Moreover, the program contains methods for filtering graphs which satisfy certain conditions involving disjoint paths. In particular, the conditions of Lemma 2, Lemma 3 and Lemma 4 of the manuscript can be checked.   

This program supports graphs with less than 64 vertices.

### Installation

This requires a working shell and `make`. Navigate to the folder containing findDominationNumber.c and compile using `make` to create a binary. 

Use `make clean` to remove the binary created in this way.

### Usage of findDominationNumber

All options can be found by executing `./findDominationNumber -h`.

Usage: `./findDominationNumber [-2|-3|-4] [-c] [-d] [-o#] [-v] [-h]`

By default this program computes the domination number of the input graphs. With `-o#` graphs with domination number `#` are output. When using either of `-2`, `-3` or `-4` the graphs which satisfy the conditions of Lemma 1, 2 and 3, respectively are output.

Graphs are read from stdin in graph6 format by default or in planar code format if `-d` is present. Graphs are sent to stdout in graph6 format or planar code depending on the presence of `-d`. If the input graph had a header, so will the output graph.

The order in which the arguments appear does not matter.
```
  -2, --K2-lemma                Check if the input graph satisfies the 
                                 conditions of Lemma 1; with -d the input 
                                 graphs are assumed to be triangulated 
                                 discs; Otherwise behavior is unexpected;
                                 without -d there may be false positives
  -3, --P3-lemma                Check if the input graph satisfies the 
                                 conditions of Lemma 2; with -d the input 
                                 graphs are assumed to be 3-connected 
                                 triangulated discs; Otherwise behavior 
                                 is unexpected; without -d there may be
                                 false positives
  -4, --P4-lemma                Check if the input graph satisfies the 
                                 conditions of Lemma 3; with -d the input
                                 graphs are assumed to be triangulated 
                                 discs; Otherwise behavior is unexpected;
                                 without -d there may be false positives
  -c, --count                   Count the number of minimal dominating
                                 sets; The min, max and avg of all input
                                 graphs will be output; cannot be used with 
                                 -2, -3 or -4
  -d, --triangulated-disc       Read the input graphs in planar code; the 
                                 boundary cycle will be computed; will lead
                                 to unexpected behavior if it is not a 
                                 triangulated disc and -2, -3 or -4 is used
  -o#, --output=#               Output graphs with domination number #; 
                                 cannot be used with -2, -3 or -4   
  -v, --verbose                 Give more detailed output
  -h, --help                    Print this help text
```

### Examples

`./findDominationNumber`
Compute the domination number of the input graphs. One can see how many graphs have a specific domination number.

`./findDominationNumber -c`
Compute the domination number of the input graphs and count the minimal dominating sets. The maximum, minimum and average number of minimal dominating sets of all input graphs is given.

`./findDominationNumber -2`
Compute the domination number of the input graphs and send to stdout if they satisfy the conditions of Lemma 1. Inputting the triangulated discs in planar code and using `-d` is faster as the boundary cycle is known in this case. Note that without using `-d` the algorithm looks at all copies of K2 in the graph, which may lead to false positive. These need to be manually inspected to figure out if the copy of K2 lies on the boundary cycle.

`./findDominationNumber -2 -v`
Same as the previous example, but also prints out which copies of K2 satisfy the conditions.

`./findDominationNumber -o#`
Compute the domination number of the input graphs and send those with domination number `#` to stdout. 

Graphs in graph6 or planar code format can be sent to stdin via a file:
`./findDominationNumber < location/of/file.g6`
`./findDominationNumber -d < location/of/file.pl`

Or directly via their graph6 code:
`./findDominationNumber <<< 'IsP@OkWHG'`

Or sent from the stdout of another program e.g. `plantri`, which outputs graphs in graph6 or planar code format:
`./plantri | ./findDominationNumber -d`