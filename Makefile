CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Iinclude
AR       := ar

SRC_DIR  := src
CLI_DIR  := cli
TEST_DIR := tests

SRCS     := $(shell find $(SRC_DIR) -name '*.cpp' -type f | sort)
OBJS     := $(SRCS:.cpp=.o)
LIB      := libsolar.a

CLI_SRC  := $(CLI_DIR)/main.cpp
CLI_BIN  := solar

TEST_SRCS := $(shell find $(TEST_DIR) -name 'test_*.cpp' -type f | sort)
TEST_BINS := $(TEST_SRCS:.cpp=)

.PHONY: all clean test

all: $(LIB) $(CLI_BIN)

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CLI_BIN): $(CLI_SRC) $(LIB)
	$(CXX) $(CXXFLAGS) $< -L. -lsolar -o $@

test: $(TEST_BINS)
	@status=0; \
	for t in $(TEST_BINS); do \
		echo "--- $$t ---"; \
		./$$t || status=1; \
	done; \
	exit $$status

$(TEST_DIR)/%: $(TEST_DIR)/%.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $< -L. -lsolar -o $@

clean:
	rm -f $(OBJS) $(LIB) $(CLI_BIN) $(TEST_BINS)
