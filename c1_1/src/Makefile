# input files
TEX=report.tex
PLOTSCRIPTS=plotscript.gpl
# C++ sources are automatically handled

# generated files
REPORT=report.pdf
PLOTS=function.pdf
RESULTS=residual.dat
PROGRAM=myprogram
OBJS=test.o sparse.o

# additional variables
CPPFLAGS=-std=c++11

all: $(PLOTS)

$(PLOTS): $(RESULTS) $(PLOTSCRIPTS)
	gnuplot plotscript.gpl

$(RESULTS): $(PROGRAM)
	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(REPORT) $(RESULTS) $(PLOTS)
