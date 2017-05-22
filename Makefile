all:
	g++ -std=c++11 -fopenmp -g main.cpp landau-wang.cpp mersenne.cpp -o landau-wang -lboost_system -lboost_filesystem -xhost -O3
clean:
	rm landau-wang
