PROJECT_NAME = homework_second

INPUT_FILE = input.txt
OUTPUT_FILE = output.txt

NPROCS ?= $(shell nproc)

IMITATIONS ?= 10
ITERATIONS ?= 100

BUILD_DIR = build

all: build test run

build:
	@echo "==> Configuring the project..."
	@cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -B$(BUILD_DIR) -H.
	@echo "==> Building the project..."
	@cmake --build $(BUILD_DIR) -j $(NPROCS)

test: build
	@echo "==> Running tests..."
	@cd $(BUILD_DIR) && ctest --verbose --parallel $(NPROCS)

run: build
	@echo "==> Running ${PROJECT_NAME}"
	@${BUILD_DIR}/${PROJECT_NAME} ${IMITATIONS} ${ITERATIONS} < ${INPUT_FILE} > ${OUTPUT_FILE}

clang-tidy: build
	@echo "==> Running clang-tidy with $(NPROCS) threads..."
	find include tests -type f \( -name '*.cpp' -o -name '*.hpp' \) -print -o -name 'main.cpp' -print \
	| xargs -P$(NPROCS) -n1 clang-tidy --p=build \
	  --extra-arg=-std=c++20

clean:
	@echo "==> Cleaning up..."
	@rm -rf $(BUILD_DIR)
	@rm -rf coverage matrix.profdata

rebuild: clean build

install:
	sudo apt-get update
	sudo apt-get install -y cmake clang libgtest-dev libboost-all-dev libeigen3-dev graphviz gnuplot

coverage: test
	@llvm-profdata merge -sparse $(BUILD_DIR)/tests/default.profraw -o matrix.profdata
	@llvm-cov show $(BUILD_DIR)/tests/matrix_test -instr-profile=matrix.profdata -format=html -show-branches=count -output-dir=coverage

rename:
	find . -not -path "./third_party/*" -not -path ".git/*" -not -path "./build/*" -not -path "." -type f | xargs sed -i 's/service_template/YOUR_SERVICE_NAME/g'

.PHONY: all build test clean rebuild install
