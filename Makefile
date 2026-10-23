CXX ?= g++
NVCC ?= nvcc
CUDA_ARCH ?= sm_120
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic -fopenmp

.PHONY: cpu cuda test test-cuda clean
cpu: build/rhscan-cpu
cuda: build/rhscan-cuda

build:
	mkdir -p build

build/rhscan-cpu: src/main.cpp src/riemann.hpp | build
	$(CXX) $(CXXFLAGS) src/main.cpp -o $@

build/rhscan-cuda: src/main.cpp src/riemann.hpp | build
	$(NVCC) -x cu -O3 -std=c++17 -arch=$(CUDA_ARCH) -DRH_CUDA src/main.cpp -o $@

build/probe: tests/probe.cpp src/riemann.hpp | build
	$(CXX) $(CXXFLAGS) tests/probe.cpp -o $@

test: cpu
	python3 tests/check.py

test-cuda: cuda cpu
	RH_SCANNER=./build/rhscan-cuda python3 tests/check.py

clean:
	rm -rf build
