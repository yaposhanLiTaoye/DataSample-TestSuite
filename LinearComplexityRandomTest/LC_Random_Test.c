#include "LC_Random_Test.h"
#include "LinearComplexity.h"
#include "my_struct.h"
#include "cephes.h"
#include <math.h>
static const int M = 500;

/*-------------------one order ju-----------------------*/
double T11[13] = { 480, 441, 423, 403, 390, 379,369, 360, 351, 341, 327, 317, 299 };//have limit in the end
double fenwei11[13] = { 0.00993521, 0.0496111, 0.0985057, 0.196668, 0.292549, 0.393766,0.498258, 0.597295, 0.694731
, 0.792997, 0.899212, 0.94834, 0.989752 };
double T12[13] = { 482, 443, 425, 404, 391, 380, 371, 362, 352, 342, 328, 318, 300 };//no limit in the end
double fenwei12[13] = { 0.00985309, 0.0489237, 0.0969601, 0.199704, 0.296058, 0.397407, 0.490929, 0.589356, 0.697183
, 0.794679, 0.899932, 0.948615, 0.989736 };
/*-------------------double order ju-----------------------*/
double T21[13] = { 1044.5, 825.5, 738.5, 653.5, 601.5, 562.5, 529.5, 500.5, 471.5, 440.5, 402.5, 375.5, 331.5 };//have limit in the end
double fenwei21[13] = { 0.00999305, 0.0496213, 0.0998859, 0.199246, 0.299574, 0.399401, 0.499772, 0.596826, 0.696678
, 0.797886, 0.899286, 0.948936, 0.989596 };
double T22[13] = { 1051.5, 830.5, 743.5, 657.5, 605.5, 566.5, 533.5, 503.5, 473.5, 443.5, 405.5, 377.5, 333.5 };//no limit in the end
double fenwei22[13] = { 0.00997701, 0.0497821, 0.099758, 0.199665, 0.299349, 0.398346, 0.497813, 0.597438, 0.699957, 0.797106, 0.89811, 0.949325, 0.989571 };
/*-------------------jump times---------------------*/
int  J_T[13] = { 143, 138, 135, 131, 129, 127, 125, 123, 121, 118, 115, 112, 107 };
double j_fenwei[13] = { 0.01057644, 0.04752037, 0.09877715, 0.21730577, 0.29879308, 0.39185401, 0.49156801, 0.59181476, 0.68637532, 0.80670430, 0.89367072, 0.94814558, 0.98818247 };



double NIST_Test(unsigned char *s, int n)
{
	int       N, K = 6;
	double    p_value, nu1, nu2, nu3, nu4, nu5, nu6, nu0, chi2;
	double    pi[7] = { 0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };

	nu0 = 0.00;
	nu1 = 0.00;
	nu2 = 0.00;
	nu3 = 0.00;
	nu4 = 0.00;
	nu5 = 0.00;
	nu6 = 0.00;
	N = n / M;
	#pragma omp parallel
	{
		#pragma omp for reduction(+:nu0,nu1,nu2,nu3,nu4,nu5,nu6)
		for (int i = 0; i < N; i++) {
			
			int L = byteLC(s , M*i, M);
			double mean = M / 2.0 + (9.0 -1) / 36.0 - (M / 3.0 + 2.0 / 9.0) / pow(2, M);
			
			double T_ = (L - mean) + 2.0 / 9.0;

			if (T_ <= -2.5)
				nu0+=1;
			else if (T_ > -2.5 && T_ <= -1.5)
				nu1+=1;
			else if (T_ > -1.5 && T_ <= -0.5)
				nu2+=1;
			else if (T_ > -0.5 && T_ <= 0.5)
				nu3+=1;
			else if (T_ > 0.5 && T_ <= 1.5)
				nu4+=1;
			else if (T_ > 1.5 && T_ <= 2.5)
				nu5+=1;
			else
				nu6+=1;
		}
	}
	chi2 = 0.00;
	chi2 += pow(nu0 - N*pi[0], 2) / (N*pi[0]);
	chi2 += pow(nu1 - N*pi[1], 2) / (N*pi[1]);
	chi2 += pow(nu2 - N*pi[2], 2) / (N*pi[2]);
	chi2 += pow(nu3 - N*pi[3], 2) / (N*pi[3]);
	chi2 += pow(nu4 - N*pi[4], 2) / (N*pi[4]);
	chi2 += pow(nu5 - N*pi[5], 2) / (N*pi[5]);
	chi2 += pow(nu6 - N*pi[6], 2) / (N*pi[6]);
	p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
	return p_value;
}
double Japan_Test(unsigned char*s, int n, int fw) {
	double p_value;
	int N = n / M, N1 = 0, N2 = 0;
	#pragma omp parallel
	{
		#pragma omp for reduction(+: N1, N2)
		for (int i = 0; i < N; i++) {
			if (byteLC(s , M*i, M) != M / 2)continue;
			double  A = byteA(s , M*i, M);
			if (A > T11[fw])N2 += 1;
			N1+=1;
		}
	}
	
	
	double p0 = N1*1.0/N;
	if(Jueduizhi(p0 - 0.5) > 1.5/sqrt(N)){
		
		return 0;
	}
	//cout << N << "," << N1 << "," << N1 << endl;
	//if (N1*Phi < 5)cout << "Too little samples!" << endl;
	double p1 = (N2*1.0) / N1;
	double z1 = (p1 - fenwei11[fw]) / sqrt(fenwei11[fw]*(1 - fenwei11[fw]) / N1);
	p_value = cephes_erfc(Jueduizhi(z1) / sqrt(2));
	
	
	return p_value;

}

double Modified_Japan_Test(unsigned char*s, int n, int fw) {
	double p_value;
	int N = n / M, N1 = 0, N2 = 0;
	#pragma omp parallel
	{
		#pragma omp for reduction(+: N1, N2)
		for (int i = 0; i < N; i++) {
			if (byteLC(s , M*i, M) != M / 2)continue;
			double  A = byteA(s , M*i, M);
			if (A > T11[fw])N2 += 1;
			N1+=1;
		}
	}
	
	
	double p1 = (N2*1.0) / N1;
	double z1 = (p1 - fenwei11[fw]) / sqrt(fenwei11[fw]*(1 - fenwei11[fw]) / N1);
	p_value = cephes_erfc(Jueduizhi(z1) / sqrt(2));
	
	
	return p_value;
}
double2 Our_Test(unsigned char *s, int n, int fw) {
	double2 p_value;
	int N = n / M;
	int N11 = 0, N12=0;
	
	#pragma omp parallel
	{
		#pragma omp for reduction(+:N11,N12)
		for (int i = 0; i < N; i++) {
			double2 D = byteAB(s , M*i, M);
			if (D.d1 > T12[fw])  N11+=1;
			if (D.d2 > T22[fw])  N12+=1;
		}
	}
	
	
	//if (N1*Phi < 5)cout << "Too little samples!"<<endl;
	//cout << N << "," << N1 << endl;
	double p1 = N11*1.0 / N;
	double z1 = (p1 - fenwei12[fw]) / sqrt(fenwei12[fw]*(1 - fenwei12[fw]) / N);
	p_value.d1 = cephes_erfc(Jueduizhi(z1) / sqrt(2));
	
	double p2 = N12*1.0 / N;
	double z2 = (p2 - fenwei22[fw]) / sqrt(fenwei22[fw]*(1 - fenwei22[fw]) / N);
	p_value.d2 = cephes_erfc(Jueduizhi(z2) / sqrt(2));
	return p_value;
}

double Jump_Test(unsigned char *s, int n, int fw) {
	double p_value;
	int N = n / M;
	int N1 = 0;
	
	#pragma omp parallel
	{
		#pragma omp for reduction(+:N1)
		for (int i = 0; i < N; i++) {
			int J = byteJumps(s , M*i, M);
			if (J > J_T[fw])  N1++;
		}
	}
	//if (N1*Phi < 5)cout << "Too little samples!"<<endl;
	//cout << N << "," << N1 << endl;
	double p1 = N1*1.0 / N;
	double z1 = (p1 - j_fenwei[fw]) / sqrt(j_fenwei[fw]*(1 - j_fenwei[fw]) / N);
	p_value = cephes_erfc(Jueduizhi(z1) / sqrt(2));
	
	return p_value;
}
