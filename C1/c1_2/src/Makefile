# input files
TEX=report.tex
PLOTSCRIPTS=plotscript.gpl
# C++ sources are automatically handled

# generated files
REPORT=report.pdf
PLOTS=function.pdf
RESULTS=data_solution0.dat data_solution1.dat
PROGRAM=myprogram
OBJS=test.o sparse.o finite.o

# additional variables
CPPFLAGS=-std=c++0x

all: $(REPORT)

$(REPORT): $(PLOTS) $(TEX)
	pdflatex -interaction=batchmode report.tex
	pdflatex -interaction=batchmode report.tex  # do latex twice

$(PLOTS): $(RESULTS) $(PLOTSCRIPTS)
	gnuplot plotscript.gpl

$(RESULTS): $(PROGRAM)
	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -O3 $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(REPORT) $(RESULTS) $(PLOTS)
