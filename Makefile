#CC = g++-9
CC = c++
#CC = cl     # for windows
CFLAGS = -O3 -std=c++11 -march=native
TARGET = cubicalripser
SRCS = cubicalripser.cpp dense_cubical_grids.cpp cube.cpp coboundary_enumerator.cpp joint_pairs.cpp compute_pairs.cpp
OBJS = cubicalripser.o dense_cubical_grids.o cube.o coboundary_enumerator.o joint_pairs.o compute_pairs.o

.PHONY: all
all: $(TARGET)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(TARGET) *.o

$(TARGET): $(OBJS) $(SRCS)
	$(CC) -o $@ $(OBJS)

