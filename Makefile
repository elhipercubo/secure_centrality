CXX ?= g++
TARGET = secure_inf
SRC = secure_inf.cpp

OPT ?= -O3 -pipe
WARN = -Wextra -Wall -Wno-shadow -Wno-vla -Wno-sign-compare -Wno-unused-result -Wstrict-aliasing=1 -pedantic-errors
CXXFLAGS ?= -std=c++11 $(WARN) $(OPT)
LDLIBS ?= -pthread -lm -lntl

all: $(TARGET)

$(TARGET): $(SRC) defs.h Makefile
	$(CXX) $(SRC) $(CXXFLAGS) $(LDLIBS) -o $@

clean:
	rm -f $(TARGET) *.o
