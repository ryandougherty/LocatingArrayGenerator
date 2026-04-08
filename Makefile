# Compiler to use
CXX = g++

# C++ flags: C++20 standard, all warnings, high optimization
CXXFLAGS = -std=c++20 -Wall -O3

# --- TBB (optional parallel support) ---
# To enable: make USE_TBB=1
# Requires libtbb-dev (Linux) or mingw-w64-ucrt-x86_64-tbb (MSYS2)
ifdef USE_TBB
    CXXFLAGS += -DHAS_TBB
    LDFLAGS = -ltbb
else
    LDFLAGS =
endif

# --- OS Specific Configuration ---
UNAME_S := $(shell uname -s 2>/dev/null)

# Default to Linux/Unix settings
TARGET = LocAG
RM = rm -f
OBJS_LIST = LocAG.o \
            phase1/phase1.o \
            phase2/phase2.o \
            phase2_ce/phase2_ce.o \
            phase2_greedy/phase2_greedy.o \
            utils/utils.o

# --- OS Overrides ---
ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
    TARGET = LocAG.exe
    RM = rm -f
else ifeq ($(findstring CYGWIN,$(UNAME_S)),CYGWIN)
    TARGET = LocAG.exe
    RM = rm -f
else ifeq ($(OS),Windows_NT)
    TARGET = LocAG.exe
    RM = del /F /Q
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
