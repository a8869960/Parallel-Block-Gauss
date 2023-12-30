FLAGS=-O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -c

all: main.o matrix_input.o process_gauss.o matrix_operations.o matrix_output.o results.o other_functions.o inverse_matrix.o for_gauss.o gauss_func.o
	g++ -pthread main.o matrix_input.o process_gauss.o matrix_operations.o matrix_output.o results.o other_functions.o inverse_matrix.o for_gauss.o gauss_func.o

main.o: main.cpp functions.h
	g++ $(FLAGS) main.cpp

matrix_input.o: matrix_input.cpp functions.h
	g++ $(FLAGS) matrix_input.cpp

process_gauss.o: process_gauss.cpp functions.h
	g++ $(FLAGS) process_gauss.cpp

matrix_operations.o: matrix_operations.cpp functions.h
	g++ $(FLAGS) matrix_operations.cpp

matrix_output.o: matrix_output.cpp functions.h
	g++ $(FLAGS) matrix_output.cpp

results.o: results.cpp functions.h
	g++ $(FLAGS) results.cpp

other_functions.o: other_functions.cpp functions.h
	g++ $(FLAGS) other_functions.cpp

inverse_matrix.o: inverse_matrix.cpp functions.h
	g++ $(FLAGS) inverse_matrix.cpp

for_gauss.o: for_gauss.cpp functions.h
	g++ $(FLAGS) for_gauss.cpp

gauss_func.o: gauss_func.cpp
	g++ $(FLAGS) gauss_func.cpp

clean:
	rm -f *.out *.o *.gch
