# input files
TEX=report.tex
PLOTSCRIPTS=plotscript.gpl
# C++ sources are automatically handled

# generated files
REPORT=report.pdf
PLOTS=plot2_0.01.pdf plot2_0.1.pdf plot1_0.1.pdf plot1_0.01.pdf
RESULTS=table_1_1 table_1_2 table_1_3 table_1_4 table_1_5 table_2_1 table_2_2 table_2_3 table_2_4 table_2_5 plot_1_0.010000.dat plot_1_0.100000.dat plot_2_0.010000.dat plot_2_0.100000.dat
PROGRAM=myprogram
OBJS=main.o

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
	gnuplot plotscript4.gpl

$(RESULTS): $(PROGRAM)
	./$(PROGRAM) 1 1 0.1 13
	./$(PROGRAM) 1 2 0.1 13
	./$(PROGRAM) 1 3 0.1 13
	./$(PROGRAM) 1 4 0.1 13
	./$(PROGRAM) 1 5 0.1 13
	./$(PROGRAM) 2 1 0.1 13
	./$(PROGRAM) 2 2 0.1 13
	./$(PROGRAM) 2 3 0.1 13
	./$(PROGRAM) 2 4 0.1 13
	./$(PROGRAM) 2 5 0.1 13


$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -O3 $(CPPFLAGS) -c $^ -o $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(REPORT) $(RESULTS) $(PLOTS)
