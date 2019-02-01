# input files
TEX=report.tex
PLOTSCRIPTS=plotscript.gpl plotscript2.gpl
# C++ sources are automatically handled

# generated files
REPORT=report.pdf
PLOTS=function.pdf function2.pdf
RESULTS=residual_delta0.dat residual_delta1.dat residual_delta2.dat residual_delta3.dat residual_delta4.dat residual_lambda0.dat residual_lambda1.dat residual_lambda2.dat residual_lambda3.dat residual_lambda4.dat
PROGRAM=myprogram
OBJS=main.o sparse.o

# additional variables
CPPFLAGS=-std=c++0x

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
	g++ -O3 $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -O3 $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(REPORT) $(RESULTS) $(PLOTS)
