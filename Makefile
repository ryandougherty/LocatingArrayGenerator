# Compiler
CXX = g++

# Flags
# Use -std=c++20 to fix the operator== error in phase2.h
CXXFLAGS = -std=c++20 -Wall -O3

# Source files
# ADD phase2_ce/phase2_ce.cpp to the list
SRCS = LocAG.cpp \
       phase1/phase1.cpp \
       phase2/phase2.cpp \
       phase2_ce/phase2_ce.cpp \
       utils/utils.cpp

# Object files (derived from SRCS)
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = LocAG

# Default rule
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET) LocAG.exe