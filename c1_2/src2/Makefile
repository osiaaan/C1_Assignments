# input files
TEX=report.tex
PLOTSCRIPTS=plotscript.gpl
# C++ sources are automatically handled

# generated files
REPORT=report.pdf
PLOTS=function.pdf
RESULTS=error0.dat error2.dat error20.dat error200.dat data_solution_0.000000_99_1.dat data_solution_0.000000_99_2.dat data_solution_0.000000_499_1.dat data_solution_0.000000_499_2.dat data_solution_0.000000_999_1.dat data_solution_0.000000_999_2.dat data_solution_2.000000_99_1.dat data_solution_2.000000_99_2.dat data_solution_2.000000_499_1.dat data_solution_2.000000_499_2.dat data_solution_2.000000_999_1.dat data_solution_2.000000_999_2.dat data_solution_20.000000_99_1.dat data_solution_20.000000_99_2.dat data_solution_20.000000_499_1.dat data_solution_20.000000_499_2.dat data_solution_20.000000_999_1.dat data_solution_20.000000_999_2.dat data_solution_200.000000_99_1.dat data_solution_200.000000_99_2.dat data_solution_200.000000_499_1.dat data_solution_200.000000_499_2.dat data_solution_200.000000_999_1.dat data_solution_200.000000_999_2.dat
PROGRAM=myprogram
OBJS=main.o sparse.o finite.o

# additional variables
CPPFLAGS=-std=c++0x

all: $(REPORT)

$(REPORT): $(PLOTS) $(TEX)
	pdflatex -interaction=batchmode report.tex
	pdflatex -interaction=batchmode report.tex  # do latex twice

$(PLOTS): $(RESULTS) $(PLOTSCRIPTS)
	gnuplot plotscript.gpl
	gnuplot plotscript2.gpl
	gnuplot plotscript3.gpl

$(RESULTS): $(PROGRAM)
	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -O3 $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(REPORT) $(RESULTS) $(PLOTS)
