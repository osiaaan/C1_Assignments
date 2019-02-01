# input files
TEX=report.tex
PLOTSCRIPTS=gnuplot.gpl gnuplot2.gpl
# C++ sources are automatically handled

# generated files
REPORT=report.pdf
PLOTS=plot1.pdf plot2.pdf
RESULTS=plot_1_0.000100.dat plot_1_0.000200.dat plot_1_0.000400.dat plot_1_0.001000.dat plot_2_0.000800.dat plot_3_0.000100.dat plot_3_0.000400.dat plot_3_0.001000.dat stabilityData_1_0.000100.dat stabilityData_1_0.000200.dat stabilityData_1_0.000400.dat stabilityData_1_0.001000.dat stabilityData_2_0.000800.dat stabilityData_3_0.000100.dat stabilityData_3_0.000400.dat stabilityData_3_0.001000.dat
PROGRAM=myprogram
OBJS=main.o vector.o sparse.o modles.o

# additional variables
CPPFLAGS=-std=c++0x

all: $(REPORT)

$(REPORT): $(PLOTS) $(TEX)
	pdflatex -interaction=batchmode report.tex
	pdflatex -interaction=batchmode report.tex  # do latex twice

$(PLOTS): $(RESULTS) $(PLOTSCRIPTS)
	gnuplot gnuplot.gpl
	gnuplot gnuplot2.gpl

$(RESULTS): $(PROGRAM)
	./$(PROGRAM) 1 0.001 16 1
	./$(PROGRAM) 1 0.0004 32 1
	./$(PROGRAM) 1 0.0001 64 1
	./$(PROGRAM) 2 0.0008 16 1
	./$(PROGRAM) 3 0.001 16 1
	./$(PROGRAM) 3 0.0004 32 1
	./$(PROGRAM) 3 0.0001 64 1
	./$(PROGRAM) 3 0.0001 64 1
	./$(PROGRAM) 1 0.0002 16 2

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -O3 $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(REPORT) $(RESULTS) $(PLOTS)
