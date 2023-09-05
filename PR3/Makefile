main: main.cpp
	g++ main.cpp -o main --std=c++20

main2: main2.cpp
	g++ main2.cpp -o main2 --std=c++20

story:
	xelatex report.tex
	xelatex report.tex
	rm *.log
	rm *.aux
	rm *.gz

run: main
	./main Adam-Bashforth 1 500000 0
	./main Adam-Bashforth 2 500000 0
	./main Adam-Bashforth 3 500000 0
	./main Adam-Bashforth 4 500000 0
	./main Adam-Bashforth 1 1000000 0
	./main Adam-Bashforth 2 1000000 0
	./main Adam-Bashforth 3 1000000 0
	./main Adam-Bashforth 4 1000000 0
	./main Adam-Bashforth 1 2000000 0
	./main Adam-Bashforth 2 2000000 0
	./main Adam-Bashforth 3 2000000 0
	./main Adam-Bashforth 4 2000000 0
	./main Adam-Moulton 1 500000 0
	./main Adam-Moulton 2 500000 0
	./main Adam-Moulton 3 500000 0
	./main Adam-Moulton 4 500000 0
	./main Adam-Moulton 1 1000000 0
	./main Adam-Moulton 2 1000000 0
	./main Adam-Moulton 3 1000000 0
	./main Adam-Moulton 4 1000000 0
	./main Adam-Moulton 1 2000000 0
	./main Adam-Moulton 2 2000000 0
	./main Adam-Moulton 3 2000000 0
	./main Adam-Moulton 4 2000000 0
	./main BDF 1 500000 0
	./main BDF 2 500000 0
	./main BDF 3 500000 0
	./main BDF 4 500000 0
	./main BDF 1 1000000 0
	./main BDF 2 1000000 0
	./main BDF 3 1000000 0
	./main BDF 4 1000000 0
	./main BDF 1 2000000 0
	./main BDF 2 2000000 0
	./main BDF 3 2000000 0
	./main BDF 4 2000000 0
	./main Classical-RK 0 50000 0
	./main Classical-RK 0 100000 0
	./main Classical-RK 0 200000 0
	./main Classical-RK 0 400000 0
	./main ESDIRK 1 50000 0
	./main ESDIRK 2 100000 0
	./main ESDIRK 3 200000 0
	./main ESDIRK 4 400000 0
	./main Gauss-Legendre 1 50000 0
	./main Gauss-Legendre 2 50000 0
	./main Gauss-Legendre 3 50000 0
	./main Gauss-Legendre 1 100000 0
	./main Gauss-Legendre 2 100000 0
	./main Gauss-Legendre 3 100000 0
	./main Gauss-Legendre 1 200000 0
	./main Gauss-Legendre 2 200000 0
	./main Gauss-Legendre 3 200000 0
	./main Gauss-Legendre 1 400000 0
	./main Gauss-Legendre 2 400000 0
	./main Gauss-Legendre 3 400000 0
	./main Fehlberg 0 100000 0
	./main Dormand-Prince 0 100000 0