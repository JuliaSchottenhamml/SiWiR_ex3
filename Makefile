CC=mpic++
CFLAGS= -march=native -mtune=native -std=c++0x -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG
LDFLAGS=
SOURCES=CGFD.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MatrixMult
COMMON=Timer.h

all: 
	$(CC) $(CFLAGS) $(SOURCES) -o cg
test: 
	$(CC) $(CFLAGS) $(SOURCES) -o cg

clean:
	rm -f *.o cg

.PHONY : all clean
