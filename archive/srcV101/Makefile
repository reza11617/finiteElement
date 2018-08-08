CC = gcc
CXX = g++
RM = rm -f
IDIR = ../include
CPPFLAGS = -I$(IDIR)
LDFLAGS = -g 
LDLIBS = -L../lib -lm -lcxsparse -lpthread
BIN = ../bin

OBJS = testBench.o \
	/Geometry/Geometry.o /Material/Material.o /Timer/Timer.o \
	/Log/Log.o /Element/Element.o

DEPS =$(subst .o,.h,$(OBJS))

testBench: $(OBJS)
	$(CXX) $(LDFLAGS) -o finite $(OBJS) $(LDLIBS) $(CPPFLAGS)

testBench.o: testBench.cpp 
	cp -r $(DEPS) ../include/
	$(CXX) $(LDFLAGS)  $(CPPFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(LDFLAGS)  $(CPPFLAGS) -c -o $@ $<

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)
	$(RM) *~

distclean: clean
	$(RM) testBench

include .depend
