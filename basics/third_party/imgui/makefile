ifeq ($(UNAME),Darwin)
CC = clang++
else
CC = g++
endif

CPPFLAGS = -O3 -g0 -std=c++11 -Wall -Wno-\#warnings 
SRCS = imgui.cpp imgui_draw.cpp imgui_demo.cpp
OBJS = $(SRCS:.cpp=.o)

dear_imgui: $(OBJS)

clean:
		rm *.o 

# Compile command
%.o: %.cpp
		$(CC) $(CPPFLAGS) -c $< -o $@
