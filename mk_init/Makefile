TARGET	= mk_init
OBJS	= mk_init.o

CXX	= g++
CXXFLAGS= -std=c++14
LDFLAGS	= -lm

%.o:%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

all:
	make $(TARGET)

run:$(TARGET)
	./$<

$(TARGET):$(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)
