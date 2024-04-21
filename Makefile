# CC = g++
# CFLAGS = -std=c++11 -Wall
# LAPACK_FLAGS = -std=c++17 -I/opt/homebrew/Cellar/lapack/3.12.0/include -L/opt/homebrew/lib -llapacke -llapack -lblas

# all: sparseMatrix lapacke1

# sparseMatrix: sparseMatrix.o denseMatrix.o
# 	$(CC) -o $@ $^

# lapacke1: lapacke1.cpp
# 	$(CC) $(CFLAGS) -o $@ $^ $(LAPACK_FLAGS)

# %.o: %.cpp
# 	$(CC) $(CFLAGS) -c $< -o $@

# clean:
# 	rm -f sparseMatrix lapacke1 *.o
CC = g++
CFLAGS = -std=c++11 -Wall
LAPACK_FLAGS = -std=c++17 -I/opt/homebrew/Cellar/lapack/3.12.0/include -L/opt/homebrew/lib -llapacke -llapack -lblas

all: lapacke1

lapacke1: lapacke1.o sparseMatrix.o denseMatrix.o
	$(CC) $(CFLAGS) -o $@ $^ $(LAPACK_FLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f lapacke1 *.o
