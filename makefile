# NEW  := MMult0
# NEW  := MMult_1x4_reg
# NEW  := MMult_4x4_vecreg
# NEW  := MMult_4x4_vecreg_subblock
NEW  := MMult_4x4_vecreg_subblock_cache
# NEW  := Strassen
# NEW := Strassen_multithread

%.o: %.c
	gcc -O2 -Wall -msse3 -c $< -o $@

all:
	make clean;
	make compare_matrix_multi.x;

compare_matrix_multi.x: compare_matrix_multi.o $(NEW).o utils.o
	gcc -pthread compare_matrix_multi.o $(NEW).o utils.o -o compare_matrix_multi.x

run:
	make all
	echo "Size,Gflops,Diff" > output_$(NEW).csv
	./compare_matrix_multi.x >> output_$(NEW).csv

clean:
	rm -f *.o *~ core *.x
