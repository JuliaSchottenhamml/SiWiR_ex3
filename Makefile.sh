CC=g++
CFLAGS= -march=native -mtune=native -std=c++0x -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG
LDFLAGS=
SOURCES=MatrixMult.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MatrixMult
COMMON=Timer.h

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
        $(CC) $(CFLAGS) $(OBJECTS) -o $@

.cpp.o: $(COMMON)
	$(CC) $(CFLAGS) $< -o $@
