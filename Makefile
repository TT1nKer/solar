CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Iinclude
DEPFLAGS := -MMD -MP
AR       := ar

SRC_DIR  := src
CLI_DIR  := cli
TEST_DIR := tests

SRCS     := $(shell find $(SRC_DIR) -name '*.cpp' -type f | sort)
OBJS     := $(SRCS:.cpp=.o)
LIB      := libsolar.a

CLI_SRCS := $(shell find $(CLI_DIR) -name '*.cpp' -type f | sort)
CLI_OBJS := $(CLI_SRCS:.cpp=.o)
CLI_BIN  := solar

TEST_SRCS := $(shell find $(TEST_DIR) -name 'test_*.cpp' -type f | sort)
TEST_BINS := $(TEST_SRCS:.cpp=)
DEPFILES  := $(OBJS:.o=.d) $(CLI_OBJS:.o=.d) $(TEST_BINS:=.d)

.PHONY: all clean test test-external-consumer

all: $(LIB) $(CLI_BIN)

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(CLI_DIR)/%.o: $(CLI_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(CLI_BIN): $(CLI_OBJS) $(LIB)
	$(CXX) $(CXXFLAGS) $(CLI_OBJS) -L. -lsolar -o $@

test: $(CLI_BIN) $(TEST_BINS)
	@status=0; \
	for t in $(TEST_BINS); do \
		echo "--- $$t ---"; \
		./$$t || status=1; \
	done; \
	exit $$status

test-external-consumer:
	./tests/test_external_consumer.sh

$(TEST_DIR)/%: $(TEST_DIR)/%.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -MF $@.d $< -L. -lsolar -o $@

$(TEST_DIR)/relativity/test_metric_cli: $(CLI_BIN)

clean:
	rm -f $(OBJS) $(CLI_OBJS) $(LIB) $(CLI_BIN) $(TEST_BINS) \
		$(DEPFILES) $(CLI_BIN).d

-include $(DEPFILES)
