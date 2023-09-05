main1: main1.cpp
	g++ main1.cpp -o main1 --std=c++20

main2: main2.cpp
	g++ main2.cpp -o main2 --std=c++20

extra: extra.cpp
	g++ extra.cpp -o extra --std=c++20

run: main1 main2
	./main1 32 I L V 5 5 > test1.txt
	./main1 32 F L V 5 5 >> test1.txt
	./main1 32 F Q V 5 5 >> test1.txt
	./main1 32 I L F 5 5 >> test1.txt
	./main1 32 F L F 5 5 >> test1.txt
	./main1 32 F Q F 5 5 >> test1.txt
	./main1 64 F Q V 5 5 >> test1.txt
	./main1 64 F Q V 5 5 >> test1.txt
	./main1 128 F Q V 5 5 >> test1.txt
	./main1 128 F Q F 5 5 >> test1.txt
	./main1 256 F Q V 5 5 >> test1.txt
	./main1 256 F Q F 5 5 >> test1.txt
	./main1 512 F Q V 5 5 >> test1.txt
	./main1 512 F Q F 5 5 >> test1.txt
	./main2 32 I L V 5 5 > test2.txt
	./main2 32 F L V 5 5 >> test2.txt
	./main2 32 F Q V 5 5 >> test2.txt
	./main2 32 I L F 5 5 >> test2.txt
	./main2 32 F L F 5 5 >> test2.txt
	./main2 32 F Q F 5 5 >> test2.txt
	./main2 64 F Q V 5 5 >> test2.txt
	./main2 64 F Q F 5 5 >> test2.txt
	./main2 128 F Q V 5 5 >> test2.txt
	./main2 128 F Q F 5 5 >> test2.txt
	./main2 256 F Q V 5 5 >> test2.txt
	./main2 256 F Q F 5 5 >> test2.txt
	./main2 512 F Q V 5 5 >> test2.txt
	./main2 512 F Q F 5 5 >> test2.txt
	./extra 32 F Q V 5 5 >> test3.txt
	./extra 32 F Q F 5 5 >> test3.txt
	./extra 64 F Q V 5 5 >> test3.txt
	./extra 64 F Q F 5 5 >> test3.txt
	./extra 128 F Q V 5 5 >> test3.txt
	./extra 128 F Q F 5 5 >> test3.txt
	./extra 256 F Q V 5 5 >> test3.txt
	./extra 256 F Q F 5 5 >> test3.txt
	./extra 512 F Q V 5 5 >> test3.txt
	./extra 512 F Q F 5 5 >> test3.txt