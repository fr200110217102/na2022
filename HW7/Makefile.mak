test: test.cpp
	g++ test.cpp -o test -O2

run: test
	./test u Normal Dirichlet 8 > test1.txt
	./test u Normal Dirichlet 16 >> test1.txt
	./test u Normal Dirichlet 32 >> test1.txt
	./test u Normal Dirichlet 64 >> test1.txt
	./test u Normal Neumann 8 >> test1.txt
	./test u Normal Neumann 16 >> test1.txt
	./test u Normal Neumann 32 >> test1.txt
	./test u Normal Neumann 64 >> test1.txt
	./test u Normal DDNN 8 >> test1.txt
	./test u Normal DDNN 16 >> test1.txt
	./test u Normal DDNN 32 >> test1.txt
	./test u Normal DDNN 64 >> test1.txt
	./test u Irnormal Dirichlet 8 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal Dirichlet 16 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal Dirichlet 32 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal Dirichlet 64 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal Neumann 8 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal Neumann 16 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal Neumann 32 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal Neumann 64 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal DDNNN 8 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal DDNNN 16 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal DDNNN 32 0.5 0.5 0.2 >> test1.txt
	./test u Irnormal DDNNN 64 0.5 0.5 0.2 >> test1.txt
	./test v Normal Dirichlet 8 > test2.txt
	./test v Normal Dirichlet 16 >> test2.txt
	./test v Normal Dirichlet 32 >> test2.txt
	./test v Normal Dirichlet 64 >> test2.txt
	./test v Normal Neumann 8 >> test2.txt
	./test v Normal Neumann 16 >> test2.txt
	./test v Normal Neumann 32 >> test2.txt
	./test v Normal Neumann 64 >> test2.txt
	./test v Normal DDNN 8 >> test2.txt
	./test v Normal DDNN 16 >> test2.txt
	./test v Normal DDNN 32 >> test2.txt
	./test v Normal DDNN 64 >> test2.txt
	./test v Irnormal Dirichlet 8 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal Dirichlet 16 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal Dirichlet 32 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal Dirichlet 64 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal Neumann 8 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal Neumann 16 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal Neumann 32 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal Neumann 64 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal DDNNN 8 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal DDNNN 16 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal DDNNN 32 0.5 0.5 0.2 >> test2.txt
	./test v Irnormal DDNNN 64 0.5 0.5 0.2 >> test2.txt
	./test w Normal Dirichlet 8 > test3.txt
	./test w Normal Dirichlet 16 >> test3.txt
	./test w Normal Dirichlet 32 >> test3.txt
	./test w Normal Dirichlet 64 >> test3.txt
	./test w Normal Neumann 8 >> test3.txt
	./test w Normal Neumann 16 >> test3.txt
	./test w Normal Neumann 32 >> test3.txt
	./test w Normal Neumann 64 >> test3.txt
	./test w Normal DDNN 8 >> test3.txt
	./test w Normal DDNN 16 >> test3.txt
	./test w Normal DDNN 32 >> test3.txt
	./test w Normal DDNN 64 >> test3.txt
	./test w Irnormal Dirichlet 8 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal Dirichlet 16 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal Dirichlet 32 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal Dirichlet 64 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal Neumann 8 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal Neumann 16 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal Neumann 32 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal Neumann 64 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal DDNNN 8 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal DDNNN 16 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal DDNNN 32 0.5 0.5 0.2 >> test3.txt
	./test w Irnormal DDNNN 64 0.5 0.5 0.2 >> test3.txt