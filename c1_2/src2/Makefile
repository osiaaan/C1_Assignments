# input files
TEX=report.tex
PLOTSCRIPTS=plotscript.gpl
# C++ sources are automatically handled

# generated files
REPORT=report.pdf
PLOTS=function_1000.pdf function_100.pdf
RESULTS=data_solution_0.000000_100_1.dat data_solution_0.000000_100_2.dat  data_solution_0.000000_1000_1.dat  data_solution_0.000000_1000_2.dat data_solution_2.000000_100_1.dat data_solution_2.000000_100_2.dat data_solution_2.000000_1000_1.dat data_solution_2.000000_1000_2.dat data_solution_20.000000_100_1.dat data_solution_20.000000_100_2.dat data_solution_20.000000_1000_1.dat data_solution_20.000000_1000_2.dat data_solution_100.000000_100_1.dat data_solution_100.000000_100_2.dat data_solution_100.000000_1000_1.dat data_solution_100.000000_1000_2.dat data_solution_200.000000_100_1.dat data_solution_200.000000_100_2.dat data_solution_200.000000_1000_1.dat data_solution_200.000000_1000_2.dat error_.dat
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
	gnuplot plotscript2.gpl

$(RESULTS): $(PROGRAM)
	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -O3 $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(REPORT) $(RESULTS) $(PLOTS)
