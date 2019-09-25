#ifndef _LinearComplexity_H_
#define _LinearComplexity_H_
#include "my_struct.h"

int  NumOneBits(unsigned i);
double Jueduizhi(double a);


unsigned char* BoolArrayToUcharArray(char s[], int start, int n);

char * UcharArrayToBoolArray(unsigned char s[], int start, int n);


int bitLC(char *s, int n);
int intLC(char *s, int n);
int byteLC(unsigned char* S, int start, int n);

double bitA(char *s, int n);
double intA(char *s, int n);
double byteA(unsigned char* S, int start, int n);

double bitB(char *s, int n);
double intB(char *s, int n);
double byteB(unsigned char* S, int start, int n);

double2 bitAB(char *s, int n);
double2 intAB(char *s, int n);
double2 byteAB(unsigned char* S, int start, int n);

int byteJumps(unsigned char*S, int start, int n);
#endif 
