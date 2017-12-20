TARGET	= liquid-rw
INIT_DIR= mk_init
INIT	= $(INIT_DIR)/init.prof
OBJS	= main.o

CC	= gcc
CXX	= g++
CXXFLAGS= -std=c++14

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

all:
	make $(TARGET)

run:$(TARGET) $(INIT)
	./$(TARGET) $(INIT)

$(TARGET):$(OBJS)
	$(CXX) -o $@ $^

$(INIT):$(INIT_DIR)/mk_init.cc
	make -C $(INIT_DIR) run
