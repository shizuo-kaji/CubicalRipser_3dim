.DEFAULT_GOAL := all
#CC = g++-9
CC = c++
#CC = cl     # for windows
CFLAGS = -O3 -std=c++11 -march=native -W
TARGET = cubicalripser
SRCS = cubicalripser.cpp dense_cubical_grids.cpp cube.cpp coboundary_enumerator.cpp joint_pairs.cpp compute_pairs.cpp
DEPS=$(SRCS:.cpp=.d)
OBJS=$(SRCS:.cpp=.o)

-include $(DEPS)

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS)

.cpp.o:
	$(CC) -MMD -MP $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f $(TARGET) *.o
