CC = c++
#CC = cl     # for windows
CFLAGS = -O3 -std=c++11 -march=native
TARGET = cubicalripser
SRCS = cubicalripser.cpp dense_cubical_grids.cpp array_index.cpp vertices.cpp birthday_index.cpp columns_to_reduce.cpp simplex_coboundary_enumerator.cpp write_pairs.cpp union_find.cpp joint_pairs.cpp compute_pairs.cpp
OBJS = cubicalripser.o dense_cubical_grids.o array_index.o vertices.o birthday_index.o columns_to_reduce.o simplex_coboundary_enumerator.o write_pairs.o union_find.o joint_pairs.o compute_pairs.o

.PHONY: all
all: $(TARGET)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(TARGET) *.o

$(TARGET): $(OBJS) $(SRCS)
	$(CC) -o $@ $(OBJS)

cubicalripser_3dim.o: cubicalripser.cpp
dense_cubical_grids.o: dense_cubical_grids.cpp dense_cubical_grids.h
array_index.o: array_index.cpp array_index.h
vertices.o: vertices.cpp vertices.h
birthday_index.o: birthday_index.cpp birthday_index.h
columns_to_reduce.o: columns_to_reduce.cpp columns_to_reduce.h
simplex_coboundary_enumerator.o: simplex_coboundary_enumerator.cpp simplex_coboundary_enumerator.h
write_pairs.o: write_pairs.cpp write_pairs.h
union_find.o: union_find.cpp union_find.h
joint_pairs.o: joint_pairs.cpp joint_pairs.h
compute_pairs.o: compute_pairs.cpp compute_pairs.h
