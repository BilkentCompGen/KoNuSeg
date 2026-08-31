CC=/usr/bin/g++-11
SRC=./src
INCL=./incl
LIB=./lib
BIN=./bin


build:
	$(CC) -O3 -std=c++17 -fopenmp $(SRC)/*.cpp -I $(INCL) -L $(LIB) -l btllib -o $(BIN)/konuseg.out

clean:
	rm -rf $(BIN)/konuseg.out

debug:
	$(CC) -std=c++17 -g -fopenmp $(SRC)/*.cpp -I $(INCL) -L $(LIB) -l btllib -o $(BIN)/konuseg.out

run_test:
	$(BIN)/konuseg.out -fa test/test.fa -fq test/test.fq -o test/test_result.bed -g 100000000 -s 11111111111111111111111111111111 -k 32 -anum 1174 -knum 6000 -fpr 0.01 -w 5 -minseg 64 -c_thresh 1
	
