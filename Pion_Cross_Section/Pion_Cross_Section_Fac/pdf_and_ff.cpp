#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "LHAPDF/LHAPDF.h"
#include "gsl/gsl_integration.h"
using namespace std;
//ofstream fout1("PDF9.txt");
//ofstream fout2("PDF356.txt");
//ofstream fout3("ff20_Qpc_fixed_c_zz.txt");
//ofstream fout4("ff50_Qpc_fixed_c_zz.txt");
int h, p;
double s; 
//p.d.f.

/*
   
The function Ctq5Pdf (Iparton, X, Q)
   returns the parton distribution inside the proton for parton [Iparton] 
   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
              CTEQ5F4 has only 4 flavors and gluon.

*/

extern "C" {
	
	double ctq5pdf_(int*, double*, double*);
	void setctq5_(int*);

}


double cteq5_pdf(int parton, double x, double scale) {
	 
	//Declaration of functions//
	
	//Declarations of variables//
	double pdfval;
	int table=3;
	
	//Initialisation of variables//
	
	//Core//
	
	//Checks
	//10^-5 < x < 1 and 1.0 < Q < 10,000 (GeV)
	if ((x < 1E-5)||(x>1.0)) {
		std::cout << "Momentum fraction outside the range: " << x << " ...\n";
		exit(1);
	}
	if ((scale < 1.0)||(scale>1.0E4)) {
		std::cout << "Scale outside the range: " << scale << " ...\n";
		exit(1);
	}
	if ((parton > 5)||(parton < -5)) {
		std::cout << "Not a valid parton: " << parton << " ...\n";
		exit(1);
	}
	if ((table > 9)||(table < 1)) {
		std::cout << "Not a valid table: " << table << " ...\n";
		exit(1);
	}
	
	//Select the table
	setctq5_(&table);
	
	//
	pdfval=ctq5pdf_(&parton,&x,&scale);
	
	return pdfval;
	
}

extern double cteq5_pdf(int parton, double x, double scale); 

//f.f.

extern "C" {
	void kkp_(int*,int*,double*,double*,double[]);
}


//subroutine kkp(ih,iset,x,qs,dh)
double kkp_ff(int hadron, int parton, double x, double scale) {
	
	//Declaration of functions//
	
	//Declarations of variables//
	double ffval, fflist[11];
	int kkp_parton,kkp_hadron;
	int order=0;
	
	//Initialisation of variables//
	
	//Core//
	
	//Checks
	if ((order != 0)&&(order != 1)) {
		std::cout << "Not a valid order: " << order << " ...\n";
		exit(1);
	}
	if ((x < 0.1)||(x>0.8)) {
		//std::cout << "Momentum fraction outside the range: " << x << " ...\n";
		//exit(1);
	}
	if ((parton < -5)||(parton > 5)) {
		std::cout << "Not a valid parton: " << parton << " ...\n";
		exit(1);
	}
	if ((scale > 100)||((scale < 1.414213562)&&(abs(parton) <= 3))||((scale < 2.9788)&&(4 == abs(parton)))||((scale < 9.46037)&&(5 == abs(parton)))) {
		//std::cout << "Scale outside the range: " << scale << " ...\n";
		//exit(1);
	}

	//
	kkp_hadron=hadron;

	//kkp_parton
	switch(parton) {
			
		case -5:
			kkp_parton=10;
			break;
		case -4:
			kkp_parton=8;
			break;
		case -3:
			kkp_parton=6;
			break;
		case -2:
			kkp_parton=4;
			break;
		case -1:
			kkp_parton=2;
			break;	
		case 0:
			kkp_parton=0;
			break;
		case 1:
			kkp_parton=1;
			break;
		case 2:
			kkp_parton=3;
			break;
		case 3:
			kkp_parton=5;
			break;
		case 4:
			kkp_parton=7;
			break;
		case 5:
			kkp_parton=9;
			break;
			
		default:
			std::cout << "Unknown parton in lhapdf(). Aborting...\n";
			exit(1);
			break;

	}
	
	kkp_(&kkp_hadron,&order,&x,&scale,fflist);
	
	ffval=fflist[kkp_parton];
	
	return ffval;

}

extern double kkp_ff(int hadron, int parton, double x, double scale); 

double integrand1(double x, void * params)
{
	return x*kkp_ff(h, p, x, s); 
}

double integrand2(double x, void * params)
{
	return x*x*kkp_ff(h, p, x, s); 
}

/*
int main(void)
{
	Plot PDF
	int i, parton[11]={0, 1, 2, 3, 4, 5, -1, -2, -3, -4, -5}; 
	float x; 
	for(x=1.1E-05; x<1; x+=1E-05)
	{
		fout1<<x<<" "; 
		for (i=0; i<11; i++)
			fout1<<x*cteq5_pdf(parton[i], x, 9)<<" ";
		fout1<<endl;
		fout2<<x<<" "; 
		for (i=0; i<11; i++)
			fout2<<x*cteq5_pdf(parton[i], x, 356)<<" ";
		fout2<<endl;
	}
	fout1.close();
	fout2.close();*/
	/*
	double r1, er1, r2, er2; 
	size_t n1, n2; 
	h = 1; 
	s = 50.0; 
	cout<<"Scale: "<<s<<"GeV"<<endl;
	cout<<"For (/Pi^++/Pi^-)/2, the average value of z"<<endl; 
	for (p=0; p<6; p++)
	{
		gsl_function gf1; 
		gf1.function = integrand1; 
		gsl_integration_qng(&gf1, 0.05, 1, 1e-10, 1e-10, &r1, &er1, &n1); 
		gsl_function gf2; 
		gf2.function = integrand2; 
		gsl_integration_qng(&gf2, 0.05, 1, 1e-10, 1e-10, &r2, &er2, &n2); 
		switch(p) 
		{
			case 0: 
				cout<<"gluon: "; 
				break; 
			case 1: 
				cout<<"u: "; 
				break; 
			case 2: 
				cout<<"d: "; 
				break; 
			case 3: 
				cout<<"s: "; 
				break; 
			case 4: 
				cout<<"c: "; 
				break; 
			case 5: 
				cout<<"b: "; 
				break; 
			default: 
				cout << "Unknown parton in lhapdf(). Aborting...\n";
				exit(1);
				break;
		}
		cout<<r2/r1<<endl; 
	}
	double r, er, R=0.0; 
	size_t n;  
	for (h=1; h<7; h++)
	{	
		gsl_function gf;
		gf.function = integrand2; 
		gsl_integration_qng(&gf, 0.05, 1, 1e-10, 1e-10, &r, &er, &n); 
		if (h==5) R+=r;
		else R+=2*r; 
	}
	cout<<R<<endl; 

	z^2D-z plot to see the difference between Q=p_h and Q=p_c
	double z, p_c; 
	h = 1;  
	p_c=20.0; 
	for (z=0.1; z<0.8; z+=1e-05)
	{
		s = p_c; 
		fout3<<z<<" ";
		for (p=-5; p<6; p++)
			fout3<<z*z*kkp_ff(h, p, z, s)<<" "; 
		fout3<<endl; 
	}
	p_c = 50.0; 
	for (z=0.1; z<=0.8; z+=1e-05)
	{
		s = p_c; 
		fout4<<z<<" "; 
		for (p=-5; p<6; p++)
			fout4<<z*z*kkp_ff(h, p, z, s)<<" "; 
		fout4<<endl; 
	}
	return 0;
	
	/*h=1; s; 
	double x=0.5; 
	for (s=1.5; s<10.0; s+=0.01)
	{
	for (p=-5; p<6; p++)
		cout<<kkp_ff(h, p, x, s)<<" ";
	cout<<endl; 
	}
	
}*/
