#include <stdio.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_integration.h"
#include "pdf_and_ff.cpp"
ofstream fout("cross_section_pi0_2.txt"); 
using namespace std; 
//#include "params.h"

double lambda=0.165, y=0.0, sqrt_S=200.0, hbar=6.582119e-16, const_c=2.99792e+08; 
int hadron=1, Nf=5; 

//alpha_s
double alpha_s(double q, double lambda, int nf) {

	// Declaration of functions //

	// Declaration of variables //
	double res;

	// Initialization of variables //

	// Core //

	//From arXiv:0802.4364, eq 3.2
	res=12.0*M_PI/( (33.0-2.0*nf) * log( pow(q,2) / pow(lambda,2) ) );
	
	return res;


}


//parton-parton cross-section
//could probably be optimised to reduce the number of if/then/else required
int two_two_type(int parton_a, int parton_b, int parton_c, int parton_d) {

	// Declaration of functions //

	// Declaration of variables //
	int type;

	// Initialisation of variables //

	type = 0;
	// Core //

	//Processes with no gluons (1,2,3,4)
	if ((parton_a != 0)&&(parton_b != 0)&&(parton_c != 0)&&(parton_d != 0)) {

		//type 1
		if ((parton_a == parton_c)&&(parton_b == parton_d)&&(abs(parton_a) != abs(parton_b))) {
			type=1;
		}
		else {

		//type -1
		if ((parton_a == parton_d)&&(parton_b == parton_c)&&(abs(parton_a) != abs(parton_b))) {
			type=-1;
		}
		else {
			
		//type 2
		if ((parton_a == parton_b)&&(parton_b == parton_c)&&(parton_c == parton_d)) {
			type=2;
		}
		else {

		//type 3
		if ((parton_a == -1*parton_b)&&(parton_c == -1*parton_d)&&(abs(parton_a) != abs(parton_c))) {
			type=3;
		}
		else {

		//type 4
		if ((parton_a == -1*parton_b)&&(parton_c == -1*parton_d)&&(parton_a == parton_c)) {
			type=4;
		}
		else {

		if ((parton_a == -1*parton_b)&&(parton_c == -1*parton_d)&&(parton_a == parton_d)) {
			type=-4;
		}

		}

		}

		}

		}
		
		}

	}
	else {

		//type 5
		if (((parton_a == 0)&&(parton_c == 0)&&(parton_b == parton_d)&&(parton_b != 0))||((parton_b == 0)&&(parton_d == 0)&&(parton_a == parton_c)&&(parton_a != 0))) {
			type=5;
		}
		else {
	
		//type -5
		if (((parton_a == 0)&&(parton_d == 0)&&(parton_b == parton_c)&&(parton_b != 0))||((parton_b == 0)&&(parton_c == 0)&&(parton_a == parton_d)&&(parton_a != 0))) {
			type=-5;
		}
		else { 
	
		//type 6
		if ((parton_a == -1*parton_b)&&(parton_c == 0)&&(parton_d == 0)&&(parton_a != 0)) {
			type=6;
		}
		else {
	
		//type 7
		if ((parton_c == -1*parton_d)&&(parton_a == 0)&&(parton_b == 0)&&(parton_c != 0)) {
			type=7;
		}
		else {
	
		//type 8
		if ((parton_a == 0)&&(parton_b == 0)&&(parton_c == 0)&&(parton_d == 0)) {
			type=8;
		}

		}
	
		}
	
		}
	
		}
	
	}

	return type;

}



//type
double twotwocs(/*void * pars, */int parton_a, int parton_b, int parton_c, int parton_d, double s, double t, double u, double scale) {

	/* Declaration of functions */
	double charge(int, int, int, int, int);

	/* Declaration of variables */
	double res;
	double sh, th, uh;
	double shs, ths, uhs;
	double valpha, valphas, vcharge;
	//params plist;
	int type;
	
	type = two_two_type(parton_a, parton_b, parton_c, parton_d); 

	/* Initialisation of variables */
	//plist= *(params *) pars;

	sh=s;
	if (/*plist.*/type < 0) {
		uh=t;
		th=u;
		type=-1*/*plist.*/type;
	}
	else {
		th=t;
		uh=u;
		type=/*plist.*/type;
	}

	shs=sh*sh;
	ths=th*th;
	uhs=uh*uh;

	valphas=alpha_s(scale, /*plist.*/lambda, /*plist.*/Nf);

	/* Core */

	//From Owens, table I
	switch (type) {

		case 1:
			res = valphas*valphas*(4.0/9.0* ( shs+uhs ) /ths);
			//group of statements 1;
			break;
		case 2:
			//group of statements 2;
			res = valphas*valphas*(4.0/9.0* ( ( shs+uhs ) /ths + ( shs+ths ) /uhs )-8.0/27.0*shs/ ( th*uh ));
			break;

		case 3:
			res = valphas*valphas*(4.0/9.0* ( ths+uhs ) /shs);
			break;

		case 4:
			res = valphas*valphas*4.0/9.0*( ( ( shs+uhs ) /ths+ ( uhs+ths ) /shs )-8.0/27.0*uhs/ ( th*sh ) );
			break;

		case 5:
			res = valphas*valphas*(-4.0/9.0* ( sh/uh+uh/sh ) + ( shs+uhs ) /ths);
			break;

		case 6:
			res = valphas*valphas*(32.0/27.0* ( th/uh+uh/th )-8.0/3.0* ( ths+uhs ) /shs);
			break;

		case 7:
			res = valphas*valphas*(1.0/6.0* ( th/uh+uh/th )-3.0/8.0* ( ths+uhs ) /shs);
			break;

		case 8:
			res = valphas*valphas*(9.0/2.0* ( 3.0-th*uh/shs-sh*uh/ths-sh*th/uhs ));
			break;
		
		default:
			std::cout << "Unknown type in twotwocs(): " << /*plist.*/type << " Aborting.../n";
			exit(1);
			break;

	}



	res=M_PI*res/shs;

	return res;

}

struct int_params
{
	int a, b, c, d; 
	double p_T, xa; 
}; 

double inner_integrand(double x_b, void * params)
{
	struct int_params *p = (struct int_params *)params; 
	double x_T = 2*(p->p_T)/sqrt_S; 
	double x_a = p->xa; 
	double z_c = x_T/(2.0*x_b)*exp(-1.0*y)+x_T/(2.0*x_a)*exp(y); 
	double s, t, u; 
	s = x_a*x_b*pow(sqrt_S, 2); 
	t = -1.0*x_a*p->p_T/z_c*sqrt_S*exp(-1.0*y); 
	u = -1.0*x_b*p->p_T/z_c*sqrt_S*exp(y); 
	return cteq5_pdf(p->a, x_a, p->p_T)*cteq5_pdf(p->b, x_b, p->p_T)*kkp_ff(hadron, p->c, z_c, p->p_T)/M_PI/z_c*twotwocs(p->a, p->b, p->c, p->d, s, t, u, p->p_T); 
}

double outer_integrand(double x_a, void * params)
{
	struct int_params *p = (struct int_params *)params; 
	struct int_params outer_set; 
	double x_bmin, x_T, r, er; 
	size_t n; 
	outer_set.xa = x_a; 
	outer_set.a = p->a; 
	outer_set.b = p->b; 
	outer_set.c = p->c; 
	outer_set.d = p->d; 
	outer_set.p_T = p->p_T; 
	x_T = 2*(p->p_T)/sqrt_S; 
	gsl_function gfi; 
	gfi.function = &inner_integrand; 
	gfi.params = &outer_set; 
	x_bmin = x_a*x_T*exp(-1.0*y)/(2*x_a-x_T*exp(y)); 
	gsl_integration_qng(&gfi, x_bmin, 1, 1e-5, 1e-5, &r, &er, &n); 
	return r; 
}
int main(void)
{
	//int parton_a, parton_b, parton_c, parton_d; 
	double x_T, x_amin, pT, R, r, er;
	size_t n; 
	struct int_params all_set; 
	gsl_function gfo; 
	gfo.function = &outer_integrand; 
	gfo.params = &all_set; 

	for (pT=4.0; pT<=15.0; pT+=0.1)
	{
		R = 0.0; 
		all_set.p_T = pT; 
		x_T = 2*(pT)/sqrt_S; 
		for (all_set.a=-5; all_set.a<=5; all_set.a++)
			for (all_set.b=-5; all_set.b<=5; all_set.b++)
				for (all_set.c=-5; all_set.c<=5; all_set.c++)
					for (all_set.d=-5; all_set.d<=5; all_set.d++) 
						if (two_two_type(all_set.a, all_set.b, all_set.c, all_set.d))
						{
							x_amin =  x_T*exp(y)/(2-x_T*exp(-1.0*y)); 
							gsl_integration_qng(&gfo, x_amin, 1, 1e-5, 1e-5, &r, &er, &n); 
							R += r; 
						}
		fout<<pT<<" "<<R*pow(hbar, 2)*pow(const_c, 2)*1e+13<<endl; 
	}
	return 0; 
}

