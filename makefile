NEW  := MMult0

%.o: %.c
	gcc -O2 -Wall -c $< -o $@

all:
	make clean;
	make compare_matrix_multi.x;

compare_matrix_multi.x: compare_matrix_multi.o $(NEW).o utils.o
	gcc compare_matrix_multi.o $(NEW).o utils.o -o compare_matrix_multi.x

run:
	make all
	# echo "version = '$(NEW)';" > output_$(NEW).m
	./compare_matrix_multi.x >>output_$(NEW).csv

clean:
	rm -f *.o *~ core *.x