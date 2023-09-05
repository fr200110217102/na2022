main1: main1.cpp
	g++ main1.cpp -o main1 --std=c++20

main2: main2.cpp
	g++ main2.cpp -o main2 --std=c++20

main3: main3.cpp
	g++ main3.cpp -o main3 --std=c++20

story:
	xelatex report.tex
	xelatex report.tex
	rm *.log
	rm *.aux
	rm *.gz

run: main
	./main1 32 48 1-32.csv
	./main1 64 96 1-64.csv
	./main1 128 192 1-128.csv
	./main1 256 384 1-256.csv
	./main2 16 32 2-16.csv
	./main2 32 64 2-32.csv
	./main2 64 128 2-64.csv
	./main2 128 256 2-128.csv
	./main3 16 16 0.1
	./main3 16 16 0.01
	./main3 16 16 0.001
	./main3 16 16 0.0001
	./main3 32 32 0.1
	./main3 32 32 0.01
	./main3 32 32 0.001
	./main3 32 32 0.0001
	./main3 64 64 0.1 3-1.csv
	./main3 64 64 0.01 3-2.csv
	./main3 64 64 0.001 3-3.csv
	./main3 64 64 0.0001 3-4.csv
	./main3 64 128 0.1 3-5.csv
	./main3 64 128 0.01 3-6.csv
	./main3 64 128 0.001 3-7.csv
	./main3 64 128 0.0001 3-8.csv
	./main3 128 128 0.1 3-9.csv
	./main3 128 128 0.01 3-10.csv
	./main3 128 128 0.001 3-11.csv
	./main3 128 128 0.0001 3-12.csv
	./main3 128 256 0.1 3-13.csv
	./main3 128 256 0.01 3-14.csv
	./main3 128 256 0.001 3-15.csv
	./main3 128 256 0.0001 3-16.csv