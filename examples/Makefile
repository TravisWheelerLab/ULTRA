# examples/Makefile — build the ULTRA example program

# C++ compiler and flags
CXX        := g++
CXXFLAGS   := -std=c++11 -Wall -I../src

# Path to the ULTRA static library
ULTRA_LIB  := ../build/libultra_core.a

# Libraries to link against
LDLIBS     := $(ULTRA_LIB) -pthread

# Your example source & target
SRC        := library_example.cpp
OBJ        := $(SRC:.cpp=.o)
TARGET     := library_example

# Default target: build the example executable
all: $(TARGET)

# Link step
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

# Compile step
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
