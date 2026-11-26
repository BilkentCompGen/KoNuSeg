CXX ?= g++
AR ?= ar
CXXFLAGS ?= -std=c++17 -fopenmp -Wall -Wextra
BTLLIB_VERSION ?= 1.4.10
BTLLIB_URL = https://github.com/bcgsc/btllib/archive/refs/tags/v$(BTLLIB_VERSION).tar.gz

SRC = ./src
INCL = ./incl
LIB = ./lib
BIN = ./bin
BTLLIB_SRC_DIR = $(LIB)/btllib-$(BTLLIB_VERSION)
BTLLIB_OBJ_DIR = $(LIB)/.btllib_obj
BTLLIB_LIB = $(LIB)/libbtllib.a

.PHONY: all build btllib clean distclean debug run_test

all: build

btllib: $(BTLLIB_LIB)

$(BTLLIB_SRC_DIR)/src/btllib/bloom_filter.cpp:
	@echo "[i] Downloading btllib v$(BTLLIB_VERSION) sources..."
	@mkdir -p $(BTLLIB_SRC_DIR)
	@if command -v curl >/dev/null 2>&1; then \
		curl -fsSL $(BTLLIB_URL) | tar -xz -C $(BTLLIB_SRC_DIR) --strip-components=1; \
	elif command -v wget >/dev/null 2>&1; then \
		wget -qO- $(BTLLIB_URL) | tar -xz -C $(BTLLIB_SRC_DIR) --strip-components=1; \
	else \
		echo "[!] Error: curl or wget is required to download btllib sources." >&2; exit 1; \
	fi

$(BTLLIB_LIB): $(BTLLIB_SRC_DIR)/src/btllib/bloom_filter.cpp
	@echo "[i] Compiling libbtllib.a with $(CXX)..."
	@rm -rf $(BTLLIB_OBJ_DIR) && mkdir -p $(BTLLIB_OBJ_DIR) $(LIB)
	cd $(BTLLIB_OBJ_DIR) && $(CXX) -c -O3 -std=c++17 -fopenmp $(abspath $(BTLLIB_SRC_DIR))/src/btllib/*.cpp -I $(abspath $(INCL)) -I $(abspath $(INCL))/cpptoml
	$(AR) rcs $(BTLLIB_LIB) $(BTLLIB_OBJ_DIR)/*.o
	@rm -rf $(BTLLIB_OBJ_DIR)

build: $(BTLLIB_LIB)
	@mkdir -p $(BIN)
	$(CXX) -O3 $(CXXFLAGS) $(SRC)/*.cpp -I $(INCL) -L $(LIB) -l btllib -o $(BIN)/konuseg.out

debug: $(BTLLIB_LIB)
	@mkdir -p $(BIN)
	$(CXX) -g $(CXXFLAGS) $(SRC)/*.cpp -I $(INCL) -L $(LIB) -l btllib -o $(BIN)/konuseg.out

clean:
	rm -rf $(BIN)/konuseg.out $(BTLLIB_OBJ_DIR) test/test_result.bed test/test_result_merged.bed test/test_result.bed_merged.bed

distclean: clean
	rm -rf $(BTLLIB_LIB) $(BTLLIB_SRC_DIR)

run_test: build
	$(BIN)/konuseg.out -fa test/test.fa -fq test/test.fq -o test/test_result.bed -g 100000000 -s 11111111111111111111111111111111 -k 32 -anum 1174 -knum 6000 -fpr 0.01 -w 5 -minseg 64 -c_thresh 1
