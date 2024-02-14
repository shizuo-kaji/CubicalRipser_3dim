.DEFAULT_GOAL := all
#CC = g++-9
CC = c++
#CC = cl     # for windows
CFLAGS = -O3 -std=c++11 -march=native -W
TARGET1 = cubicalripser
TARGET2 = tcubicalripser
SRCS1 = cubicalripser.cpp dense_cubical_grids.cpp cube.cpp coboundary_enumerator.cpp joint_pairs.cpp compute_pairs.cpp
SRCS2 = cubicalripser.cpp dense_cubical_grids_T.cpp cube.cpp coboundary_enumerator_T.cpp joint_pairs.cpp compute_pairs.cpp
DEPS=$(SRCS:.cpp=.d)
OBJS1=$(SRCS1:.cpp=.o)
OBJS2=$(SRCS2:.cpp=.o)

-include $(DEPS)

.PHONY: all
all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJS1)
	$(CC) -o $(TARGET1) $(OBJS1)

$(TARGET2): $(OBJS2)
	$(CC) -o $(TARGET2) $(OBJS2)

.cpp.o:
	$(CC) -MMD -MP $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f $(TARGET1) $(TARGET2) *.o
