TARGET	= mk_init
OBJS	= mk_init.o

CXX	= g++
CXXFLAGS= -std=c++14
LDFLAGS	= -lm

%.o:%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET):$(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)
