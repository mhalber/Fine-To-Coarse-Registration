#Base Settings
UNAME = $(shell uname -s)

ifeq ($(UNAME),Darwin)
CC = clang++
else
CC = g++
endif

CPPFLAGS = -O3 -g0 -std=c++11 -Wall 

#System specific warnings
ifeq ($(UNAME),Linux)
CPPFLAGS += -Wno-sign-compare -Wno-unused-result #TODO: gcc complains about these, add fixes in code
endif

ifeq ($(UNAME),Darwin)
CPPFLAGS += -Wno-\#warnings
endif


# Executable name
EXE_NAME = fetregister

# Useful Directories
BIN_DIR = ../bin/
LIB_DIR = /usr/local/lib/
OUT_DIR = ../out/
BSC_DIR = basics/
GAPS_DIR = gaps/
THIRD_PARTY_DIR = $(BSC_DIR)third_party/
IMGUI_DIR = $(THIRD_PARTY_DIR)imgui/

# List of source files
SRCS = fetregister.cpp
OBJS = $(SRCS:%.cpp=$(OUT_DIR)%.o)
OBJS += ../$(IMGUI_DIR)imgui.o ../$(IMGUI_DIR)imgui_draw.o ../$(IMGUI_DIR)imgui_demo.o

# Compile and link options
CPPFLAGS += -I. -I/usr/include -I/usr/local/include -I../$(IMGUI_DIR) -I../$(BSC_DIR) -I../$(GAPS_DIR)pkgs -DRN_USE_CSPARSE -DRN_USE_SPLM -DNDEBUG

 # Libraries
ifeq ($(UNAME),Linux)
OPENGL_LIBS =  -L/usr/lib64 -lGLEW -lGL -lGLU 

ifeq ($(USE_WINDOW), 1)
LIBS = -L$(LIB_DIR) -lglfw3 -lrt -lm -ldl -lXrandr -lXinerama -lXext -lXcursor -lXrender -lXfixes -lX11 -lpthread -lxcb -lXau -lXdmcp
CPPFLAGS += -DFETREGISTER_USE_WINDOW
LIBS += -L../$(GAPS_DIR)lib/x86_64/ -lFET -lR3Surfels -lR3Graphics -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -ljpeg -lpng -lsplm -lCSparse
else
LIBS = -L../$(GAPS_DIR)lib/x86_64/ -lFET -lR3Surfels -lR3Graphics -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -ljpeg -lpng -lsplm -lCSparse
endif

LIBS += $(OPENGL_LIBS)

endif

ifeq ($(UNAME),Darwin)
LIBS = -L../$(GAPS_DIR)/lib/x86_64 -lFET -lR3Surfels -lR3Graphics -lR3Shapes -lR2Shapes -lRNMath -lRNBasics -ljpeg -lpng -lsplm -lCSparse 
FRAMEWORKS = -framework opengl 

ifeq ($(USE_WINDOW), 1)
LIBS += -L$(LIB_DIR) -lglfw3 
CPPFLAGS += -DFETREGISTER_USE_WINDOW
endif

endif

# Make targets
all: clean fetregister

fetregister: $(OBJS)
		$(CC) $(FRAMEWORKS) $(CPPFLAGS) $(OBJS) -o ${BIN_DIR}${EXE_NAME} $(LIBS)

clean:
		rm -f $(OUT_DIR)*.a $(OUT_DIR)*.o ${BIN_DIR}${EXE_NAME} ${BIN_DIR}${EXE_NAME}.exe


# Compile command
$(OUT_DIR)%.o: %.cpp
		$(CC) $(CPPFLAGS) -c $< -o $@


# GNU Make: targets that don't build files
.PHONY: all clean
