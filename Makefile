TARGET	= liquid-rw
OBJS	= main.o

CC	= gcc
CXX	= g++
CXXFLAGS= -std=c++14

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

all:
	make $(TARGET)

run:$(TARGET)
	./$<

$(TARGET):$(OBJS)
	$(CXX) -o $@ $^
