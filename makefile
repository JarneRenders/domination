compiler=gcc
flags=-std=gnu11 -march=native -Wall -Wno-missing-braces -O3
profileflags=-std=gnu11 -march=native -Wall -g -pg -fsanitize=address

# The 64-bit version of this program is faster but only supports graphs up to 64 vertices.
64bit: findDominationNumber.c readGraph/readGraph6.c bitset.h 
	$(compiler) -DUSE_64_BIT -o findDominationNumber findDominationNumber.c readGraph/readGraph6.c $(flags)

profile: findDominationNumber.c readGraph/readGraph6.c bitset.h 
	$(compiler) -DUSE_64_BIT -o findDominationNumber-pr findDominationNumber.c readGraph/readGraph6.c $(profileflags)

all: 64bit 

.PHONY: clean
clean:
	rm -f findDominationNumber findDominationNumber-pr

