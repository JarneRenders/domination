/**
 * findDominationNumber.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 *
 */

#define USAGE "Usage: ./findDominationNumber [-2|-3|-4] [-c] [-d] [-o#] [-v] [-h]"
#define HELPTEXT "Helptext:\n\
By default this program computes the domination number of the input graphs.\n\
With `-o#` graphs with domination number `#` are output. When using either of\n\
`-2`, `-3` or `-4` the graphs which satisfy the conditions of Lemma 1, 2 and 3,\n\
respectively, are output.\n\
\n\
Graphs are read from stdin in graph6 format by default or in planar code format\n\
if `-d` is present. Graphs are sent to stdout in graph6 format or planar code\n\
depending on the presence of `-d`. If the input graph had a header, so will the\n\
output graph.\n\
\n\
The order in which the arguments appear does not matter.\n\
\n\
  -2, --K2-lemma                Check if the input graph satisfies the\n\
                                 conditions of Lemma 1; with -d the input\n\
                                 graphs are assumed to be triangulated\n\
                                 discs; Otherwise behavior is unexpected;\n\
                                 without -d there may be false positives\n\
  -3, --P3-lemma                Check if the input graph satisfies the\n\
                                 conditions of Lemma 2; with -d the input\n\
                                 graphs are assumed to be 3-connected\n\
                                 triangulated discs; Otherwise behavior\n\
                                 is unexpected; without -d there may be\n\
                                 false positives\n\
  -4, --P4-lemma                Check if the input graph satisfies the\n\
                                 conditions of Lemma 3; with -d the input\n\
                                 graphs are assumed to be triangulated\n\
                                 discs; Otherwise behavior is unexpected;\n\
                                 without -d there may be false positives\n\
  -c, --count                   Count the number of minimal dominating\n\
                                 sets; The min, max and avg of all input\n\
                                 graphs will be output; cannot be used with\n\
                                 -2, -3 or -4\n\
  -d, --triangulated-disc       Read the input graphs in planar code; the\n\
                                 boundary cycle will be computed; will lead\n\
                                 to unexpected behavior if it is not a\n\
                                 triangulated disc and -2, -3 or -4 is used\n\
  -o#, --output=#               Output graphs with domination number #;\n\
                                 cannot be used with -2, -3 or -4\n\
  -v, --verbose                 Give more detailed output\n\
  -h, --help                    Print this help text\n"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "readGraph/readGraph6.h"
#include "bitset.h"

//******************************************************************************
//
//                          Structs
//
//******************************************************************************

struct graph {
    bitset *adjacencyList;
    char *graphString;
    int numberOfVertices;
    int *orderedNbrs;
    int *outerCycle;
    int lenOuterCycle;
};

struct options {
    bool verboseFlag;
    bool countFlag;
    bool deleteK1Flag;
    bool deleteK2Flag;
    bool deleteP3Flag;
    bool deleteP4Flag;
    bool triangulatedDiscFlag;
    bool writeHeader;
};

// For recursively building paths between path.start and path.end
struct path {
    int start;
    int end;
    bitset vertices;
    int lastAdded;
    bool isFinished;
};

//******************************************************************************
//
//                  Macros for dealing with graphs
//
//******************************************************************************

#define addEdge(adjacencyList, v, w) {\
    add(adjacencyList[v], w);\
    add(adjacencyList[w], v);\
}

#define removeEdge(adjacencyList, v, w){\
    removeElement(adjacencyList[v], w);\
    removeElement(adjacencyList[w], v);\
}

//  Used when the graph struct includes a planar embedding (ordering of the
//  edges incident to a vertex).
#define getNextEdge(g, i, j) (g)->orderedNbrs[(g)->numberOfVertices*(i) + (j)];


//******************************************************************************
//
//                      Writing graphs
//
//******************************************************************************

//  Graphical representation of graph to stderr
void printGraph(bitset adjacencyList[], int numberOfVertices) {
    for(int i = 0; i < numberOfVertices; i++) {
        fprintf(stderr, "%d: ", i);
        forEach(nbr, adjacencyList[i]) {
            fprintf(stderr, "%d ", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

//  Graphical representation of embedding to stderr
void printOrderedGraph(struct graph *g) {
    for(int i = 0; i < g->numberOfVertices; i++) {
        fprintf(stderr, "%d: ", i);
        int start = next(g->adjacencyList[i], -1);
        int curr = start;
        do {
            fprintf(stderr, "%d ", curr);
            curr = getNextEdge(g, i, curr);
        } while(curr != start);
        fprintf(stderr, "\n");
    }
}

//  Graph6 code to stdout
void writeToG6(bitset adjacencyList[], int numberOfVertices) {
    char graphString[8 + numberOfVertices*(numberOfVertices - 1)/2];
    int pointer = 0;

    //  Save number of vertices in the first one, four or 8 bytes.
    if(numberOfVertices <= 62) {
        graphString[pointer++] = (char) numberOfVertices + 63;
    }
    else if(numberOfVertices <= 258047) {
        graphString[pointer++] = 63 + 63;
        for(int i = 2; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else if(numberOfVertices <= 68719476735) {
        graphString[pointer++] = 63 + 63;
        graphString[pointer++] = 63 + 63;
        for(int i = 5; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else {
        fprintf(stderr, "Error: number of vertices too large.\n");
        exit(1);
    }

    // Group upper triangle of adjacency matrix in groups of 6. See B. McKay's 
    // graph6 format.
    int counter = 0;
    char charToPrint = 0;
    for(int i = 1; i < numberOfVertices; i++) {
        for(int j = 0; j < i; j++) {
            charToPrint = charToPrint << 1;
            if(contains(adjacencyList[i], j)) {
                charToPrint |= 1;
            }
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
                charToPrint = 0;
                counter = 0;
            }
        }
    }

    //  Pad final character with 0's.
    if(counter != 0) {
        while(counter < 6) {
            charToPrint = charToPrint << 1;
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
            }
        }
    }

    //  End with newline and end of string character.
    graphString[pointer++] = '\n';
    graphString[pointer++] = '\0';
    printf("%s", graphString);
}

//******************************************************************************
//
//             Simple methods for determining the connectivity
// 
//******************************************************************************

void connectivityDFS(bitset adjacencyList[], int numberOfVertices,
bitset *remainingVertices, bitset *checkedVertices, int vertexToCheck) {
    removeElement(*remainingVertices, vertexToCheck);
    add(*checkedVertices, vertexToCheck);
    forEach(nbr,
     intersection(adjacencyList[vertexToCheck], *remainingVertices)) {
        connectivityDFS(adjacencyList, numberOfVertices, remainingVertices,
         checkedVertices, nbr);
    }
}

bool is1Connected(bitset adjacencyList[], int numberOfVertices, bitset
deletedVertices) {
    bitset remainingVertices = complement(EMPTY, numberOfVertices);
    remainingVertices = difference(remainingVertices, deletedVertices);
    bitset checkedVertices = EMPTY;

    connectivityDFS(adjacencyList, numberOfVertices, &remainingVertices,
     &checkedVertices, next(remainingVertices, -1));
    if(isEmpty(remainingVertices)) {
        return true;
    }
    return false;
}

bool isKConnected(bitset adjacencyList[], int numberOfVertices, bitset
deletedVertices, int K) {
    if(K == 1) {
        return is1Connected(adjacencyList, numberOfVertices, deletedVertices);
    }

    if(!is1Connected(adjacencyList, numberOfVertices, deletedVertices)) {
        return false;
    }

    for(int i = 0; i < numberOfVertices; i++) {
        if(contains(deletedVertices, i)) {
            continue;
        }
        if(!isKConnected(adjacencyList, numberOfVertices,
         union(deletedVertices, singleton(i)), K-1)) {
            return false;
        }
    }
    return true;
}

//******************************************************************************
//
//                 Methods for finding domination number
//
//******************************************************************************

// Could maybe be faster if we dynamically keep track of dominated vertices
bool isDominatingSet(struct graph* g, bitset includedVertices, bitset set) {

    //  Loop over includedVertices not in set (those in set are of course
    //  dominated).
    forEach(vertex, difference(includedVertices, set)) {
        if(isEmpty(intersection(g->adjacencyList[vertex], set))) { 
            return false;
        }
    }
    return true;
}

void generateAllDominatingSets(struct graph *g, struct options *options, bitset
includedVertices, bitset currentSet, int *sizeOfSmallestDomSet, int lastAdded,
long long unsigned int *nOfCoveringsForDominationNumber) {

    if(isDominatingSet(g, includedVertices, currentSet)) {
        if(size(currentSet) < *sizeOfSmallestDomSet) { // Keep track of smallest
            *sizeOfSmallestDomSet = size(currentSet);
            (*nOfCoveringsForDominationNumber) = 0;
        }
        if(size(currentSet) == *sizeOfSmallestDomSet) {
            (*nOfCoveringsForDominationNumber)++;
        }

        // Do not return here if we want all dominating sets, can return here if
        // we only want all dominating sets of minimal size.
        return; 
    }

    //  If we only want to determine the domination number and not count, we can
    //  prune if the set is non-dominating and has size one less then the
    //  smallest found so far.
    if(!options->countFlag && size(currentSet) >= *sizeOfSmallestDomSet - 1) {
        return;
    }

    //  Otherwise return when size is same as smallest.
    if(size(currentSet) >= *sizeOfSmallestDomSet) {
        return;
    }

    //  Pruning criterium: if v is not yet in set and < lastAdded it cannot be
    //  in set
    for(int v = 0; v < lastAdded; v++) {
        if(contains(currentSet, v)) {
            continue;
        }

        //  v has no neighbours in set
        if(!isEmpty(intersection(g->adjacencyList[v], currentSet))) {
            continue;
        }

        //  no included neighbour of v will be in set in this branch of the
        //  recursion tree
        if(next(intersection(g->adjacencyList[v], includedVertices), lastAdded)
         == -1) {
            return;
        }
    }

    //  Add vertices which are included and not yet in set. lastAdded is the
    //  vertex of highest index in the set. Each iteration starts the branch in
    //  which vertex>lastAdded is added and no element between vertex and
    //  lastAdded can be present.
    forEachAfterIndex(vertex, difference(includedVertices, currentSet),
     lastAdded) {
        add(currentSet, vertex);
        generateAllDominatingSets(g, options, includedVertices, currentSet,
         sizeOfSmallestDomSet, vertex, nOfCoveringsForDominationNumber);
        removeElement(currentSet, vertex);
    }
}

int findDominationNumber(struct graph *g, struct options *options, bitset
excludedVertices, long long unsigned int *nOfCoveringsForDominationNumber) {
    int dominationNumber = INT_MAX;
    bitset includedVertices = complement(excludedVertices, g->numberOfVertices);
    generateAllDominatingSets(g, options, includedVertices, EMPTY,
     &dominationNumber, -1, nOfCoveringsForDominationNumber);
    return dominationNumber;
}

//******************************************************************************
//
//         Methods for checking Lemma conditions (see manuscript)
//
//******************************************************************************

//  Check if Lemma 1 holds, i.e. given a triangulated disk, is there a K2 on the
//  boundary such that its removal has the same domination number. Warning: If
//  g was not read in planar code we return true if it holds for any K2
//  (not necessarily on the boundary.)
bool K2LemmaHolds(struct graph *g, struct options *options,
 int dominationNumber) {
    bool holds = false;
    for(int i = 0; i < g->numberOfVertices; i++) {
        forEachAfterIndex(nbr, g->adjacencyList[i], i) {
            bitset excludedVertices = union(singleton(i), singleton(nbr));
            long long unsigned int dummy;
            if(findDominationNumber(g, options, excludedVertices, &dummy) ==
             dominationNumber) {
                if(options->verboseFlag) {
                    fprintf(stderr, "After removing %d %d\n", i, nbr);
                }
                holds = true;
            }
        }
    }
    return holds;
}

//  Given a set of pairs of start and end points (the array of path structs)
//  start with the first pair and recursively build a path from start to end,
//  with the remaining vertices do the same for the second pair, then the
//  third, etc. If it is not possible, backtrack until it is or until we have
//  tried all possibilities.
bool existPaths(struct graph *g, bitset remainingVertices, struct
path *pathArray[], int pathArrayLength, int currentPathIdx) {

    //  If last added is neighbour of path, the path is complete.
    if(contains(g->adjacencyList[pathArray[currentPathIdx]->lastAdded],
     pathArray[currentPathIdx]->end)) {

        //  If not the final path, start building the next one.
        if(currentPathIdx < pathArrayLength - 1) {
            return existPaths(g, remainingVertices, pathArray, pathArrayLength,
             currentPathIdx + 1);
        }
        return true;
    }

    //  Add remaining nbr to the last added and see if this can complete the
    //  path. 
    struct path *currentPath = pathArray[currentPathIdx];
    int currentVertex = currentPath->lastAdded;
    forEach(nbr,
     intersection(g->adjacencyList[currentVertex],remainingVertices)) {

        // Add nbr to path
        currentPath->lastAdded = nbr;
        bitset newRemainingVertices = 
         difference(remainingVertices, singleton(nbr));
        add(pathArray[currentPathIdx]->vertices, nbr);
        if(existPaths(g, newRemainingVertices, pathArray, pathArrayLength,
         currentPathIdx)) {
            return true;
        }

        //  Could not find a path, so try a new neighbour 
        removeElement(pathArray[currentPathIdx]->vertices, nbr);
        currentPath->lastAdded = currentVertex;
    }
    return false;
}

//  For every vertex i in ends, create a path struct with end i and start v for
//  each vertex v in starts. Then execute the recursive function to see if
//  these paths exist such that they do not share any vertex except the end. 
bool containsDisjointPaths(struct graph *g, bitset starts, bitset ends) {
    forEach(i, ends) {

        //  Create a path for every start which is not i.
        struct path *pathArray[size(difference(starts, singleton(i)))];
        bitset remainingVertices = difference(
         complement(singleton(i), g->numberOfVertices), starts);
        int j = 0;
        forEach(v, starts) {
            struct path *path = malloc(sizeof(struct path));
            path->start = v;
            path->end = i;
            path->vertices = union(singleton(v), singleton(i));
            path->lastAdded = v;
            pathArray[j++] = path;
        }
        bool disjointPathsExist = existPaths(g, remainingVertices,
         pathArray, size(difference(starts, singleton(i))), 0);
        for(int k = 0; k < size(difference(starts, singleton(i))); k++) {
            free(pathArray[k]);
        }
        if(!disjointPathsExist) {
            return false;
        }
    }
    return true;
}

//  Check if the domination number of G-X is the same as dominationNumber for
//  each subset X of excludedVertices. Uses a bitset trick to iterate over the
//  subsets of excluded vertices, but this breaks the usual 128 bit
//  implementations. (Can be fixed if we define some macro for -1).
bool removalOfSubsetsGivesSameDominationNumber(struct graph *g,
 struct options *options, bitset excludedVertices, int dominationNumber) {
    bitset set = excludedVertices;
    while(!isEmpty(set)) {
        long long unsigned int temp = 0;
        int newDom = findDominationNumber(g, options, set, &temp);
        if(newDom != dominationNumber) {
            return false;
        }
        set = intersection((set - 1), excludedVertices);
    }
    return true;
}

//  Check if Lemma 2 holds, under the assumption graphs are 3-connected
//  triangulated disks. Warning: If g was not read in planar code we return
//  true if it holds for any P3 (not necessarily on the boundary.) 
bool P3LemmaHolds(struct graph *g, struct options *options,
 int dominationNumber) {
    bool found = false;

    //  If -d is present and graph is not a triangulation.
    if(g->lenOuterCycle > 3) { 

        //  Loop over all 3-paths on boundary cycle.
        for (int i = 0; i < g->lenOuterCycle; i++) {

            //  (i) Check if removal of every subset of the vertices of the
            //  3-path leaves the same domination number.
            bitset excludedVertices = singleton(g->outerCycle[i]);
            add(excludedVertices, g->outerCycle[(i+1)%g->lenOuterCycle]);
            add(excludedVertices, g->outerCycle[(i+2)%g->lenOuterCycle]);
            if(!removalOfSubsetsGivesSameDominationNumber(g, options,
             excludedVertices, dominationNumber)) {
                continue;
            }


            //  (iii) Check if for every vertex x not on the path whether there
            //  are paths from each vertex on P3 to x such that they are
            //  internally disjoint.  
            bitset ends = difference(complement(EMPTY, g->numberOfVertices),
             excludedVertices);
            if(!containsDisjointPaths(g, excludedVertices, ends)) {
                continue;
            }
            if(options->verboseFlag) {
                fprintf(stderr, "P3: %d %d %d \n", 
                 g->outerCycle[i], g->outerCycle[(i+1)%g->lenOuterCycle], 
                 g->outerCycle[(i+2)%g->lenOuterCycle]);
            }
            found = true;
        }
        return found;
    }

    //  If we are here then either -d was not given so we do not know which is
    //  the outer cycle or g is a triangulation and then any face can be the
    //  outer cycle. Either way, we check all 3-paths. 

    if(options->verboseFlag) {
        fprintf(stderr,
         "No large enough outer cycle found, checking all 3-paths:\n");
    }
    for(int v = 0; v < g->numberOfVertices; v++) {
        forEach(u, g->adjacencyList[v]) {
            forEachAfterIndex(w, g->adjacencyList[v], u) {
                bitset excludedVertices = singleton(v);
                add(excludedVertices, u);
                add(excludedVertices, w);

                if(!removalOfSubsetsGivesSameDominationNumber(g, options,
                 excludedVertices, dominationNumber)) {
                    continue;
                }

                bitset ends = difference(complement(EMPTY, g->numberOfVertices),
                 excludedVertices);
                if(!containsDisjointPaths(g, excludedVertices, ends)) {
                    continue;
                }
                if(options->verboseFlag) {
                    fprintf(stderr, "P3: %d %d %d\n",u,v,w);
                }
                found = true;
            }
        }
    }
    return found;
}

//  Check if Lemma 3 holds. Warning: If g was not read in planar code we return
//  true if it holds for any P4 (not necessarily on the boundary.) 
bool P4LemmaHolds(struct graph *g, struct options *options,
 int dominationNumber) {
    bool found = false;

    //  Check 4-connectivity. (Plantri cannot output these graphs without some
    //  filter.)
    if(!isKConnected(g->adjacencyList, g->numberOfVertices, EMPTY, 4)) {
        if(options->verboseFlag) {
            fprintf(stderr, "Graph is not 4-connected.\n");
        }
        return false;
    }

    //  If -d and not a triangulation.
    if(g->lenOuterCycle > 3) {
        for (int i = 0; i < g->lenOuterCycle; i++) {
            bitset excludedVertices = singleton(g->outerCycle[i]);
            add(excludedVertices, g->outerCycle[(i+1)%g->lenOuterCycle]);
            add(excludedVertices, g->outerCycle[(i+2)%g->lenOuterCycle]);
            add(excludedVertices, g->outerCycle[(i+3)%g->lenOuterCycle]);

            //  (i) Check if removal of every subset of the vertices of the
            //  4-path leaves the same domination number.
            if(!removalOfSubsetsGivesSameDominationNumber(g, options,
             excludedVertices, dominationNumber)) {
                continue;
            }


            //  (iii) Check if for every vertex x not on the path whether there
            //  are paths from each vertex on P4 to x such that they are
            //  internally disjoint.  
            bitset ends = difference(complement(EMPTY, g->numberOfVertices),
             excludedVertices);
            if(!containsDisjointPaths(g, excludedVertices, ends)) {
                continue;
            }


            // (iv) Check if for every vertex y on the 4-path, there exist paths
            // to each of the other vertices which are internally disjoint.
            if(!containsDisjointPaths(g, excludedVertices, excludedVertices)) {
                continue;
            }

            if(options->verboseFlag) {
                fprintf(stderr, "P4: %d %d %d %d \n",
                 g->outerCycle[i], 
                 g->outerCycle[(i+1)%g->lenOuterCycle], 
                 g->outerCycle[(i+2)%g->lenOuterCycle], 
                 g->outerCycle[(i+3)%g->lenOuterCycle]);
            }
            found = true;
        }
        return found;
    }

    if(g->lenOuterCycle == 3) {
        if(options->verboseFlag) {
            fprintf(stderr, "Triangulation has no P4 on outer cycle.\n");
        }
        return false;
    }
    if(options->verboseFlag) {
        fprintf(stderr,
         "No outer cycle found, checking all 4-paths:\n");
    }
    for(int v = 0; v < g->numberOfVertices; v++) {
        forEach(u, g->adjacencyList[v]) {
            forEach(w, g->adjacencyList[v]) {
                forEach(x, difference(g->adjacencyList[w], singleton(v))) {
                    bitset excludedVertices = singleton(u);
                    add(excludedVertices, v);
                    add(excludedVertices, w);
                    add(excludedVertices, x);
                    if(!removalOfSubsetsGivesSameDominationNumber(g, options,
                     excludedVertices, dominationNumber)) {
                        continue;
                    }
                    bitset ends = 
                     difference(complement(EMPTY, g->numberOfVertices), 
                     excludedVertices);
                    if(!containsDisjointPaths(g, excludedVertices, ends)) {
                        continue;
                    }

                    if(!containsDisjointPaths(g, excludedVertices,
                     excludedVertices)) {
                        continue;
                    }
                    
                    fprintf(stderr, "%d %d %d %d \n", u,v,w,x);
                    found = true;
                }
            }
        }
    }
    return found;
}

//******************************************************************************
//
//                          Methods for I/O
//
//******************************************************************************

//  Find the outer cycle of a triangulated disc if its length is at least 4.
//  Assumes g is a triangulated disc.
int findOuterCycle(struct graph *g, struct options *options) {
    if(!options->triangulatedDiscFlag) {
        fprintf(stderr, "Error: not the right format for triangulated disc.\n");
        exit(1);
    }
    int outerCycle[g->numberOfVertices];
    int lenCycle = 0;
    for(int i = 0; i < g->numberOfVertices; i++) {
        int start = next(g->adjacencyList[i], -1);
        int prev = start;
        do { 
            int curr = getNextEdge(g, i, prev);
            if(!contains(g->adjacencyList[prev], curr) ||
             size(g->adjacencyList[i]) == 2) {
                outerCycle[0] = prev;
                outerCycle[1] = i;
                outerCycle[2] = curr;
                lenCycle = 2;
                break;
            }
            prev = curr;
        } while(prev != start);
        if(lenCycle > 0) {
            break;
        } 
    }
    if(lenCycle == 0) {
        if(options->verboseFlag) {
            fprintf(stderr,
             "This is a triangulation so every face can be outer cycle.\n");
        }
        return 3;
    }

    // Get clockwise order.
    int next = getNextEdge(g, outerCycle[2], outerCycle[1]);
    if(contains(g->adjacencyList[outerCycle[1]], next) &&
     size(g->adjacencyList[outerCycle[2]]) != 2) {
        int temp = outerCycle[0];
        outerCycle[0] = outerCycle[2];
        outerCycle[2] = temp;
    }
    while(!contains(g->adjacencyList[outerCycle[0]], outerCycle[lenCycle])) {
        next = getNextEdge(g, outerCycle[lenCycle], outerCycle[lenCycle - 1]);
        outerCycle[++lenCycle] = next;
    }
    lenCycle++;

    if(options->verboseFlag) {
        fprintf(stderr, "Outer Cycle:\n");
        for(int i = 0; i < lenCycle; i++) {
            fprintf(stderr, "%d ", outerCycle[i]);
        }
        fprintf(stderr, "\n");
    }
    g->outerCycle = malloc((lenCycle) * sizeof(int));
    memcpy(g->outerCycle, outerCycle, (lenCycle) * sizeof(int));
    return lenCycle;
}

int readPlanarCode(struct graph *g, struct options *options) {
    unsigned char character;
    if(fread(&character, sizeof(unsigned char), 1, stdin) == 1) {
        // Not necessarily a head, might be a 62 vertex graph without header
        if(character == '>') {  
            unsigned short header[6];
            header[0] = character;
            header[1] = (unsigned short) getc(stdin);
            header[2] = (unsigned short) getc(stdin);
            if((header[1] == '>') && (header[2] == 'p')) { 

                // Should be a planarcode header now
                while((character = getc(stdin)) != '<');
                character = getc(stdin); // Read second time
                if(character != '<') {
                    fprintf(stdout, "Problems with header -- single '<'\n");
                    exit(1);
                }
                if(!fread(&character, sizeof(unsigned char), 1, stdin)) {
                    exit(0);
                } 
            }
        }

        if(character != 0) {
            unsigned short code[6];
            code[0] = character;
            g->numberOfVertices = code[0];
            g->adjacencyList = calloc(g->numberOfVertices,sizeof(bitset));
            g->orderedNbrs =
             malloc(g->numberOfVertices*g->numberOfVertices*sizeof(int*));
            for(int v = 0; v < code[0]; v++) {
                int nbr = (int) getc(stdin);
                if(!nbr) continue;
                add(g->adjacencyList[v], nbr - 1);
                add(g->adjacencyList[nbr - 1], v);
                int prevNbr = nbr;
                int firstNbr = nbr;
                while((nbr = (int) getc(stdin))) {
                    add(g->adjacencyList[v], nbr - 1);
                    add(g->adjacencyList[nbr - 1], v);
                    g->orderedNbrs[v*g->numberOfVertices + prevNbr - 1] =
                     nbr - 1;
                    prevNbr = nbr; 
                }
                g->orderedNbrs[v*g->numberOfVertices + prevNbr - 1] =
                 firstNbr - 1;
            }
        }
        int lenCycle = findOuterCycle(g, options);
        g->lenOuterCycle = lenCycle;

        return 0;
    }
    return 1;
}

int readGraph6(struct graph *g, struct options *options, char *graphString) {
    g->numberOfVertices = getNumberOfVertices(graphString);
    g->lenOuterCycle = 0;
    if(g->numberOfVertices == -1 || g->numberOfVertices > MAXVERTICES) {
        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }
        return 1;
    }
    g->adjacencyList = calloc(g->numberOfVertices, sizeof(bitset));
    if(loadGraph(graphString, g->numberOfVertices, g->adjacencyList) == -1) {
        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }
        return 1;
    }
    return 0;
}

void writePlanarCode(struct graph *g, bool header) {
    // Write header
    if(header) {
        printf(">>planar_code<<");
    }
    // Write number of vertices
    printf("%c",(unsigned char)g->numberOfVertices);
    // Write adjacency list
    for (int v = 0; v < g->numberOfVertices; v++) {
        int start = next(g->adjacencyList[v], -1);
        int curr = start;
        do {
            printf("%c",(unsigned char) (curr + 1));
            curr = getNextEdge(g, v, curr);
        } while(curr != start);
        printf("%c",0);
    }
}

void writeGraph(struct graph *g, struct options *options) {
    if(options->triangulatedDiscFlag) {
        // printOrderedGraph(g);
        writePlanarCode(g, options->writeHeader);
        options->writeHeader = false;
    }
    else {
        printf("%s", g->graphString);
    }
}

void freeGraph(struct graph *g, struct options *options) {
    free(g->adjacencyList);
    if(options->triangulatedDiscFlag) {
        free(g->orderedNbrs);
        if(g->lenOuterCycle > 3) {
            free(g->outerCycle);
        }
    }
    else {
        free(g->graphString);
    }
}

//******************************************************************************
//
//                              MAIN
//
//******************************************************************************

int main(int argc, char ** argv) {
    struct options options = {0};
    options.writeHeader = true;
    int opt;
    int output = -1;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"K2-lemma", no_argument, NULL, '2'},
            {"P3-lemma", no_argument, NULL, '3'},
            {"P4-lemma", no_argument, NULL, '4'},
            {"count", no_argument, NULL, 'c'},
            {"triangulated-disc", no_argument, NULL, 'd'},
            {"help", no_argument, NULL, 'h'},
            {"output", required_argument, NULL, 'o'},
            {"verbose", no_argument, NULL, 'v'}
        };

        opt = getopt_long(argc, argv, "234cdho:v", long_options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case '2':
                options.deleteK2Flag = true;
                fprintf(stderr,
                 "Outputting graphs in which the removal of a K2 has the same domination number as the original.\n");
                break;
            case '3':
                options.deleteP3Flag = true;
                fprintf(stderr,
                 "Outputting graphs in which Lemma 2 properties are satisfied.\n");
                break;
            case '4':
                options.deleteP4Flag = true;
                fprintf(stderr,
                 "Outputting graphs in which Lemma 3 properties are satisfied.\n");
                break;
            case 'c':
                options.countFlag = true;
                break;
            case 'd':
                options.triangulatedDiscFlag = true;
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'o':
                output = (int) strtol(optarg, (char **)NULL, 10);
                break;
            case 'v':
                options.verboseFlag = true;
                break;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,
                 "Use ./findDominationNumber --help for more detailed instructions.\n");
                return 1;
        }
    }

    unsigned long long int counter = 0;
    unsigned long long int skippedGraphs = 0;
    unsigned long long int passedGraphs = 0;
    unsigned long long int frequencies[MAXVERTICES] =
     { [ 0 ... MAXVERTICES-1 ] = 0 };
    unsigned long long int minimumNOfCoverings = INT_MAX;
    unsigned long long int maximumNOfCoverings = 0;
    unsigned long long int sumOfAllNsOfCoverings = 0;
    clock_t start = clock();

    // Check if more than one of these variables is set.
    if(((int)(options.deleteK2Flag ? 1:0) + (int)(options.deleteP3Flag ? 1:0) +
     (int)(options.deleteP4Flag ? 1:0) + 
     (int)((output != -1 || options.countFlag) ? 1:0)) >= 2) {
        fprintf(stderr,
         "Error: only select one of -1, -2, -3, -4. Do not combine with -c or -o#.\n");
        exit(1);
    }

    if((options.deleteK2Flag || options.deleteP3Flag || options.deleteP4Flag) &&
     !options.triangulatedDiscFlag) {
        fprintf(stderr,
         "Warning: the search is not restricted to the boundary cycle.");
        fprintf(stderr, " Can give false positives.\n");
    }

    //  Start looping over lines of stdin.
    while(true) {
        struct graph g;
        if(options.triangulatedDiscFlag) {
            if(readPlanarCode(&g, &options) == 1) {
                break;
            }
            if(options.verboseFlag) {
                fprintf(stderr, "Looking at: \n");
                printOrderedGraph(&g);
            }
        }
        else {
            g.graphString = NULL;
            size_t size;
            if(getline(&(g.graphString), &size, stdin) == -1) {
                free(g.graphString);
                break;
            }
            if(readGraph6(&g, &options, g.graphString) == 1) {
                skippedGraphs++;
                continue;
            }
            if(options.verboseFlag) {
                fprintf(stderr, "Looking at %s\n", g.graphString);
            }
        }

        //  Loaded g at this point, determine its domination number.
        counter++;
        long long unsigned int nOfCoveringsForDominationNumber = 0;

        int dominationNumber = findDominationNumber(&g, &options, EMPTY,
         &nOfCoveringsForDominationNumber);
        frequencies[dominationNumber]++;

        if(dominationNumber == output) {
            passedGraphs++; 
            writeGraph(&g, &options);
        }

        if(options.deleteK2Flag) {
            if(K2LemmaHolds(&g, &options, dominationNumber)) {
                writeGraph(&g, &options);
                passedGraphs++;
            }
        }
        else if(options.deleteP3Flag) {
            if(P3LemmaHolds(&g, &options, dominationNumber)) {
                passedGraphs++;
                writeGraph(&g, &options);
            }
        }
        else if(options.deleteP4Flag) {
            if(P4LemmaHolds(&g, &options, dominationNumber)) {
                passedGraphs++;
                writeGraph(&g, &options);
            }
        }

        if(options.countFlag) {
            sumOfAllNsOfCoverings += nOfCoveringsForDominationNumber;
            if(nOfCoveringsForDominationNumber < minimumNOfCoverings) {
                minimumNOfCoverings = nOfCoveringsForDominationNumber;
            }
            if(nOfCoveringsForDominationNumber > maximumNOfCoverings) {
                maximumNOfCoverings = nOfCoveringsForDominationNumber;
            }
        }        
        freeGraph(&g, &options);
    }

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    for (int i = 0; i < MAXVERTICES; ++i) {
        if(frequencies[i] != 0) {
            fprintf(stderr, "\n \t%16lld graphs with domination number %d",
             frequencies[i], i);
        }
    }
    fprintf(stderr, "\n");

    if(options.countFlag) {
        fprintf(stderr,
         "Minimum number of minimal dominating sets: %llu\n", 
         minimumNOfCoverings);
        fprintf(stderr,
         "Maximum number of minimal dominating sets: %llu\n", 
         maximumNOfCoverings);
        fprintf(stderr,
         "Average number of minimal dominating sets: %llu\n", 
         sumOfAllNsOfCoverings/(counter-skippedGraphs));
    }

    fprintf(stderr,"\rChecked %lld graphs in %f seconds.\n", counter, time_spent);
    fprintf(stderr, "Written %llu\n", passedGraphs);

    return 0;
}