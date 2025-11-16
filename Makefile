# Compiler to use
CXX = g++

# C++ flags: C++20 standard, all warnings, high optimization
CXXFLAGS = -std=c++20 -Wall -O3

# --- OS Specific Configuration ---
# Check for Unix-like systems (Linux, macOS, WSL, Git Bash, etc.)
# The 'shell' command is used to run 'uname'
UNAME_S := $(shell uname -s)

# Default to Linux/Unix settings
TARGET = LocAG
LDFLAGS = -ltbb
RM = rm -f
# We will use / for all paths, as g++ and make handle this well
OBJS_LIST = LocAG.o \
            phase1/phase1.o \
            phase2/phase2.o \
            phase2_ce/phase2_ce.o \
            phase2_greedy/phase2_greedy.o \
            utils/utils.o

# --- OS Overrides ---
# Check for MINGW (e.g., Git Bash on Windows)
ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
    TARGET = LocAG.exe
    RM = rm -f
# Check for Cygwin
else ifeq ($(findstring CYGWIN,$(UNAME_S)),CYGWIN)
    TARGET = LocAG.exe
    RM = rm -f
# Check for native Windows (cmd.exe)
else ifeq ($(OS),Windows_NT)
    TARGET = LocAG.exe
#     LDFLAGS = -ltbb
    LDFLAGS = 
    RM = del /F /Q
    # Use backslashes for native Windows 'del' command
    OBJS_LIST = LocAG.o \
                phase1\phase1.o \
                phase2\phase2.o \
                phase2_greedy\phase2_greedy.o \
                phase2_ce\phase2_ce.o \
                utils\utils.o
endif

# All the source files (using Unix-style paths)
SRCS = LocAG.cpp \
       phase1/phase1.cpp \
       phase2/phase2.cpp \
       phase2_greedy/phase2_greedy.cpp \
       phase2_ce/phase2_ce.cpp \
       utils/utils.cpp

# Object files are derived from source files, replacing .cpp with .o
OBJS = $(SRCS:.cpp=.o)

# By default, the first target is built. 'all' is a common name.
all: $(TARGET)

# The linking rule: combines all .o files into the final executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# --- Generic Compilation Rules ---
# VPATH tells 'make' to look for .cpp files in these directories
VPATH = phase1:phase2:phase2_ce:phase2_greedy:utils

# This single generic rule builds .o files from .cpp files
# It automatically finds .cpp files in the VPATH directories
# $< is the .cpp file (prerequisite)
# $@ is the .o file (target)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# This rule handles the main LocAG.cpp in the root directory
LocAG.o: LocAG.cpp
	$(CXX) $(CXXFLAGS) -c LocAG.cpp -o LocAG.o

# Rule to clean up build files
clean:
	-$(RM) $(TARGET)
	-$(RM) $(OBJS_LIST)

# Tell 'make' that 'all' and 'clean' are not actual files
.PHONY: all clean
