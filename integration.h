#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <random>
#include <iomanip>
#include <thread>
#include "timer.h"

using namespace std;

#define R 1
#define RUN_NUM 10
#define MM 1000000
#define MAX_N 1000000000 //max unsigned long: 4.3*(10^9)
#define MAX_ATTEMPT 10
#define MAX_THREAD 8

typedef struct{

	double N;
	double answer;
	double error;
	double upper_bound;
	double lower_bound;
	double run_time;

} RESULT;

typedef struct{

	//input
	int d;
	double N;

	//output
	double p; //ratio of hypersphere in hypercube

} MC;

typedef struct{

	//input
	int d;
	double k_1, K_1;

	//output
	double p_1, n; //ratio of I in X
	double I, O, X;

} CB;

class Integration{

public:
	Integration(int d, double tol);
	~Integration(){};
	static double even_fact(double r, int k, double pi)
		{return k == 0 ? 1 : (pi*r*r/(double)k)*even_fact(r, k-1, pi);}
	static double odd_fact(double r, int k, double pi)
		{return k == 0 ? 1 : (2.0*pi*r*r/(2*k-1))*odd_fact(r, k-1, pi);}
	static double Ball_Volume(double r, int d, double pi){

		int k = d/2;
		if(d%2 == 0) return even_fact(r, k, pi);
		return (2.0*r/d)*odd_fact(r, k, pi);
	}
	static double K2N(double K, int d){return pow(K, (double)d);}
	static double N2K(double N, int d){return pow(N, 1.0/(double)d);}

protected:
	void Pass_Result(RESULT *r, double N);
	void Pass_PI_Result(RESULT *r, double N);
	static double Calculate_PI_1(double answer, int d);
	static void   Time_Estimator(double N, double unit_N);

	int d; //dimension
	double tol;
	double C; //hypercube volume
	double answer, error, upper_bound, lower_bound, run_time;
	int attempt_num;
};

class Monte_Carlo : public Integration{

public:
	Monte_Carlo(int d, int M, double tol, double t);
	~Monte_Carlo(){};
	void Generate_Result(RESULT *r);
	void Generate_Result(RESULT *r, double N);
	void Calculate_PI(RESULT *r, double N);

private:
	void Calculate_Volume(double N);
	void Init_MC(MC *mc, double N);
	static void Thread_Work(void *arg);
	static double Is_Inside(double v[]);

	int M; //M: # of runs
	double t; //t value for the given risk of being wrong
};

class Cube_based : public Integration{

public:
	Cube_based(int d, double tol);
	~Cube_based(){};
	void Generate_Result(RESULT *r);
	void Generate_Result(RESULT *r, double N);
	void Calculate_PI(RESULT *r, double N);

private:
	void Calculate_Volume(double N);
	void Init_CB(CB *cb, double K);
	void Sum_CB(CB *cb, double &n);
	static void Thread_Work(void *arg);
	static double Is_Inside(double node[], int d);
	static void Init_DFS(double **parent, double *node, double *path, CB *cb);
	static void Free_DFS(double **parent, double *node, double *path, int d);
	static double Count_Cubes(int &l, double **parent, double *node, double *path, CB *cb);
	static void Copy_Node(double dest[], double src[], int j, int d);
	static void Update_Node(double node[], int k, double x, int d, double K_1);

	double I, O, X, p_1;
};
