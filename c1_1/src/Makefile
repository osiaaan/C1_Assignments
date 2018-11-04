# input files
TEX=report.tex
PLOTSCRIPTS=plotscript.gpl plotscript2.gpl
# C++ sources are automatically handled

# generated files
REPORT=report.pdf
PLOTS=function.pdf function2.pdf
RESULTS=residual.dat
PROGRAM=myprogram
OBJS=test.o sparse.o

# additional variables
CPPFLAGS=-std=c++11

all: $(REPORT)

$(REPORT): $(PLOTS) $(TEX)
	pdflatex -interaction=batchmode report.tex
	pdflatex -interaction=batchmode report.tex  # do latex twice

$(PLOTS): $(RESULTS) $(PLOTSCRIPTS)
	gnuplot plotscript.gpl
	gnuplot plotscript2.gpl

$(RESULTS): $(PROGRAM)
	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(REPORT) $(RESULTS) $(PLOTS)
