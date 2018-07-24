all: 
	gfortran -c kkp.f Ctq5Pdf.f
PDF: 
	g++ -c -I/home/tianyu/include/LHAPDF/include pdf_and_ff.cpp 
	g++ -Wall -o pdf_dist kkp.o Ctq5Pdf.o pdf_and_ff.o -lgfortran
FF: 
	LD_LIBRARY_PATH=/home/tianyu/lib/gsl-2.4/lib
	export LD_LIBRARY_PATH
	g++ -c -I/home/tianyu/lib/gsl-2.4/include -I/home/tianyu/lib/LHAPDF/include pdf_and_ff.cpp
	g++ -Wall -o ff_dist kkp.o -L/home/tianyu/lib/gsl-2.4/lib Ctq5Pdf.o pdf_and_ff.o -lgfortran -lgsl -lgslcblas -lm
FF_Q: 
	g++ -c -I/home/tianyu/lib/LHAPDF/include pdf_and_ff.cpp
	g++ -Wall -o ff_Q kkp.o Ctq5Pdf.o pdf_and_ff.o -lgfortran
hadron_production : 
	LD_LIBRARY_PATH=/home/tianyu/lib/gsl-2.4/lib
	export LD_LIBRARY_PATH
	#g++ -Wall -o ghloprod kkp.o Ctq5Pdf.o pdf_and_ff.o parton_cross_sections.o main.cpp -lgsl -lblas -lg2c -lm
	g++ -c -I/home/tianyu/lib/gsl-2.4/include -I/home/tianyu/lib/LHAPDF/include parton_cross_sections.cpp
	g++ -Wall -o cross_sections -L/home/tianyu/lib/gsl-2.4/lib kkp.o Ctq5Pdf.o parton_cross_sections.o -lgfortran -lgsl -lgslcblas -lm
