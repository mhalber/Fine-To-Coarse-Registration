#Base Settings
UNAME = $(shell uname -s)

ifeq ($(UNAME),Darwin)
CC = clang++
else
CC = g++
endif

CPPFLAGS = -O3 -g0 -std=c++11 -Wall -Wno-\#warnings -Wno-sign-compare

# Executable name
EXE_NAME = fetbenchmark

# Useful Directories
BIN_DIR = ../bin/
LIB_DIR = /usr/local/lib/
OUT_DIR = ../out/
BSC_DIR = basics/

# List of source files
SRCS = fetbenchmark.cpp
OBJS = $(SRCS:%.cpp=$(OUT_DIR)%.o)

# Compile and link options
CPPFLAGS += -I. -I/usr/local/include  -I../$(BSC_DIR) 

# Make targets
all: clean fetbenchmark

fetbenchmark: $(OBJS)
		$(CC) $(FRAMEWORKS) $(CPPFLAGS) $(OBJS) -o ${BIN_DIR}${EXE_NAME}

clean:
		rm -f $(OUT_DIR)*.a $(OUT_DIR)*.o ${BIN_DIR}${EXE_NAME} ${BIN_DIR}${EXE_NAME}.exe


# Compile command
$(OUT_DIR)%.o: %.cpp
		$(CC) $(CPPFLAGS) -c $< -o $@


# GNU Make: targets that don't build files
.PHONY: all clean
