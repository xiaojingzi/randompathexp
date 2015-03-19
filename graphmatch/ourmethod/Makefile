main_debug: main.cc
	g++ main.cc -o main -std=c++11 -Wall -g -ggdb

main_release: main.cc
	g++ main.cc -o $@ -std=c++11 -Wall -O3

test: main
	echo "./prog dataset T D Epsilon"
	./main_release kdd 5 50 0.01
clean:
	rm -rf a.out main random_project
