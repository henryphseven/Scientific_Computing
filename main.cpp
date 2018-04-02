#include "integration.h"

using namespace std;

int col_num = 9;

map<double,double> tmap; //t-distribution table
map<double,double>::iterator it;

void build_tmap(void){

	//assume M = 10 => degree of freedom = 10 - 1 = 9
	//source 1: http://www.statisticshowto.com/tables/t-distribution-table/
	//source 2: https://www.medcalc.org/manual/t-distribution.php
	tmap.insert(pair<double,double>(0.2,1.383));
	tmap.insert(pair<double,double>(0.1,1.833));
	tmap.insert(pair<double,double>(0.05,2.262));
	tmap.insert(pair<double,double>(0.02,2.821));
	tmap.insert(pair<double,double>(0.01,3.250));
	tmap.insert(pair<double,double>(0.005,3.690));
	tmap.insert(pair<double,double>(0.002,4.297));
	tmap.insert(pair<double,double>(0.001,4.781));
}

void Auto_Calculate(int part, double N, int D);
void Insert_Tupe(double tuple[], int d, double exact_answer, RESULT *r_mc, RESULT *r_cb);
void Print_Tupe(double tuple[]);
void Result_Tuple(double tuple[], RESULT *r);
void Print_Result(double tuple[]);

int main(int argc, char *argv[])
{
	char input[20];
	char exit_project[] = "exit";
	int part;

    build_tmap();

	while(1){

		do{

			printf("Please specify Part Number: (or input exit to quit)\n");
			printf("1 2 3 4(unused) 5(auto run Part 2) 6(auto run Part 3)\n");
			cin.getline(input, sizeof(input));
			if(strcmp(input, exit_project) == 0) return 0;
			part = atoi(input);

		} while(part < 1 || part > 6);

		if(part == 4) continue;

		double N = 0, tol = 0, alpha = 0, t = 0;
		int precision = 0, M = 0;

		if(part == 1){

			M = RUN_NUM;

			do{

				printf("\nPlease specify precision in terms of digits (e.g., 4 means tol = 0.001): (or input exit to quit)\n");
				cin.getline(input, sizeof(input));
				if(strcmp(input, exit_project) == 0) return 0;
				precision = atoi(input);

			} while(precision <= 0);

			tol = pow(10.0, -(double)(precision-1));

			do{

				printf("\nPlease specify alpha (%%) for Monte Carlo Method: (or input exit to quit)\n");
				for(it = tmap.begin(); it != tmap.end(); ++it) printf("%g ", 100*(it->first));
				printf("\n");
				cin.getline(input, sizeof(input));
				if(strcmp(input, exit_project) == 0) return 0;
				alpha = atof(input)/100;
				it = tmap.find(alpha);

			} while(it == tmap.end());

			t = it->second;
		}
		else{

			//you are allowed to use only 1 sample for Part 2 + 3
			M = 1;

			do{

				printf("\nPlease specify N (input 0 for default value of 1 million): (or input exit to quit)\n");
				cin.getline(input, sizeof(input));
				if(strcmp(input, exit_project) == 0) return 0;
				N = atof(input);
				if(N == 0) N = MM;

			} while(N <= 0);
		}

		//auto run all d from 1 to D
		if(part == 5 || part == 6){

			int D;

			do{

				printf("\nPlease specify max dimension (D): (or input exit to quit)\n");
				cin.getline(input, sizeof(input));
				if(strcmp(input, exit_project) == 0) return 0;
				D = atoi(input);

			} while(D <= 0);

			Auto_Calculate(part, N, D);
			continue;
		}

		int d;

		do{

			printf("\nPlease specify dimension (d): (or input exit to quit)\n");
			cin.getline(input, sizeof(input));
			if(strcmp(input, exit_project) == 0) return 0;
			d = atoi(input);

		} while(d <= 0);

		RESULT r_mc, r_cb;
		double exact_answer, diffs_mc, diffs_cb;

		switch(part){

		case 1:
			printf("\nPart 1: Given precision = %d and alpha = %g%%, find %d-ball volume\n", precision, alpha*100, d);
			break;

		case 2:
			printf("\nPart 2: Given N = %e, find %d-ball volume\n", N, d);
			break;

		case 3:
			printf("\nPart 3: Given N = %e and d = %d, find pi\n", N, d);
			break;
		}

		if(part == 3) exact_answer = M_PI;
		else exact_answer = Integration::Ball_Volume((double)R, d, M_PI);

		Monte_Carlo mc(d, M, tol, t);

		switch(part){

		case 1:
			mc.Generate_Result(&r_mc);
			break;

		case 2:
			mc.Generate_Result(&r_mc, N);
			break;

		case 3:
			mc.Calculate_PI(&r_mc, N);
			break;
		}

		printf("\nThe result of Monte Carlo:\n");
		printf("The answer is %.3lf(%.3lf) when N = %e and M = %d\n", r_mc.answer, r_mc.error, r_mc.N, M);
		if(part == 1) printf("The true value is within [%.3lf, %.3lf] at %.1lf%% confidence level\n",
				r_mc.lower_bound, r_mc.upper_bound, (1-alpha)*100);
		printf("The number of seconds of run time: %.0lf\n", r_mc.run_time);

		diffs_mc = fabs(r_mc.answer - exact_answer);
		printf("The correct answer is %.3lf and the diffs is %e\n", exact_answer, diffs_mc);

		int cube_based, ch_part = 0;

		if(part == 1){

			printf("\nDo you want to continue with Cube-based? (May need to wait for a lifetime): (or input exit to quit)\n");
			printf("Yes: input 1, No: input 0\n");
			cin.getline(input, sizeof(input));
			if(strcmp(input, exit_project) == 0) return 0;
			cube_based = atoi(input);

			if(cube_based == 0) continue;

			do{

				printf("\nPlease specify precision in terms of digits (e.g., 4 means tol = 0.001)");
				printf("\nor input 0 if you would like to specify N instead: (or input exit to quit)\n");
				cin.getline(input, sizeof(input));
				if(strcmp(input, exit_project) == 0) return 0;
				precision = atoi(input);

			} while(precision < 0);

			if(precision == 0){

				do{

					printf("\nPlease specify N (input 0 for default value of 1 million): (or input exit to quit)\n");
					cin.getline(input, sizeof(input));
					if(strcmp(input, exit_project) == 0) return 0;
					N = atof(input);
					if(N == 0) N = MM;

				} while(N <= 0);

				part = 2;
				ch_part = 1;
				tol = 0;
			}
			else tol = pow(10.0, -(double)(precision-1));
		}

		Cube_based cb(d, tol);

		switch(part){

		case 1:
			cb.Generate_Result(&r_cb);
			break;

		case 2:
			cb.Generate_Result(&r_cb, N);
			break;

		case 3:
			cb.Calculate_PI(&r_cb, N);
			break;
		}

		if(ch_part) part = 1;

		printf("\nThe result of Cube-based:\n");
		printf("The answer is %.3lf(%.3lf) when N = %e and K = %.0lf\n", r_cb.answer, r_cb.error, r_cb.N, Integration::N2K(r_cb.N, d));
		printf("The true value is guaranteed to be within [%.3lf, %.3lf]\n",
				r_cb.lower_bound, r_cb.upper_bound);
		printf("The number of seconds of run time: %.0lf\n", r_cb.run_time);

		diffs_cb = fabs(r_cb.answer - exact_answer);
		printf("The correct answer is %.3lf and the diffs is %e\n\n", exact_answer, diffs_cb);

		printf("Comparison:\n");

		printf("[Value] ");
		printf("Monte Carlo: %lf, Cube-based: %lf, diffs: %e\n", r_mc.answer, r_cb.answer, fabs(r_mc.answer - r_cb.answer));

		printf("[Run Time] ");
		if(r_mc.run_time == r_cb.run_time)
			printf("Monte Carlo (%.0lf) and Cube-based (%.0lf) are equally efficient\n", r_mc.run_time, r_cb.run_time);
		else {
			if(r_mc.run_time < r_cb.run_time)
				printf("Monte Carlo (%.0lf) is more efficient than Cube-based (%.0lf)", r_mc.run_time, r_cb.run_time);
			else
				printf("Cube-based (%.0lf) is more efficient than Monte Carlo (%.0lf)", r_cb.run_time, r_mc.run_time);
			printf(" by %.0lf second(s)\n", fabs(r_mc.run_time - r_cb.run_time));
		}

		printf("[Accuracy] ");
		if(diffs_mc == diffs_cb)
			printf("Monte Carlo (%.0lf) and Cube-based (%.0lf) are equally accurate\n", diffs_mc, diffs_cb);
		else{
			if(diffs_mc < diffs_cb)
				printf("Monte Carlo (%e) is more accurate than Cube-based (%e)", diffs_mc, diffs_cb);
			else
				printf("Cube-based (%e) is more accurate than Monte Carlo (%e)", diffs_cb, diffs_mc);
			printf(" by %e\n", fabs(diffs_mc - diffs_cb));
		}

		printf("\n");

		//prepare summary table for this particular d
		double table[col_num];
		double mc_result[5], cb_result[5];

		Insert_Tupe(table, d, exact_answer, &r_mc, &r_cb);
		if(part == 1){

			Result_Tuple(mc_result, &r_mc);
			Result_Tuple(cb_result, &r_cb);
		}

		printf("When d = %d:\n", d);
		printf("d,Exact_Answer,MC_Answer,CB_Answer,Answer_Dffs,");
		printf("MC_Diffs,CB_Diffs,Accuracy_Winner,Accuracy_Diffs,");
		if(part == 1){

			printf("MC_Time,CB_Time,Efficiency_Winner,Time_Diffs,");
			printf("MC_Answer,MC_Error,MC_LBound,MC_UBound,MC_N,MC_M,MC_Time,");
			printf("CB_Answer,CB_Error,CB_LBound,CB_UBound,CB_N,CB_K,CB_Time\n");
		}

		Print_Tupe(table);
		if(part == 1){

			printf(",%.0lf,%.0lf,",r_mc.run_time,r_cb.run_time);
			if(r_mc.run_time == r_cb.run_time) printf("PAR,");
			else if(r_mc.run_time < r_cb.run_time) printf("MC,");
			else  printf("CB,");
			printf("%.0lf,",fabs(r_mc.run_time - r_cb.run_time));
			Print_Result(mc_result);
			printf(",%d,%.0lf,", M, r_mc.run_time);
			Print_Result(cb_result);
			printf(",%.0lf,%.0lf\n", Integration::N2K(r_cb.N, d), r_cb.run_time);
		}

		printf("\n");
	} //while(1)

	return 0;
}

void Auto_Calculate(int part, double N, int D){

	if(!(part == 5 || part == 6)) return;

	int i, d;
	double exact_answer;
	RESULT r_mc, r_cb;

	double **table;
	table = (double**)malloc(sizeof(double*)*D);
	for(i = 0; i < D; i++) table[i] = (double*)malloc(sizeof(double)*col_num);

	for(d = 1; d <= D; d++){

		if(part == 5) exact_answer = Integration::Ball_Volume((double)R, d, M_PI);
		else if(part == 6) exact_answer = M_PI;

		Monte_Carlo mc(d, 1, 0, 0);
		Cube_based cb(d, 0);

		if(part == 5){

			mc.Generate_Result(&r_mc, N);
			cb.Generate_Result(&r_cb, N);
		}
		else if(part == 6){

			mc.Calculate_PI(&r_mc, N);
			cb.Calculate_PI(&r_cb, N);
		}

		Insert_Tupe(table[d-1], d, exact_answer, &r_mc, &r_cb);
	}

	printf("\nWhen N = %e:\n", N);
	printf("d,Exact_Answer,MC_Answer,CB_Answer,Answer_Dffs,");
	printf("MC_Diffs,CB_Diffs,Accuracy_Winner,Accuracy_Diffs\n");

	for(d = 1; d <= D; d++){

		Print_Tupe(table[d-1]);
		printf("\n");
	}

	for(i = 0; i < D; i++) free(table[i]);
	free(table);
}

void Insert_Tupe(double tuple[], int d, double exact_answer, RESULT *r_mc, RESULT *r_cb){

	tuple[0] = d;
	tuple[1] = exact_answer;
	tuple[2] = r_mc->answer;
	tuple[3] = r_cb->answer;
	tuple[4] = fabs(r_mc->answer - r_cb->answer);
	tuple[5] = fabs(exact_answer - r_mc->answer);
	tuple[6] = fabs(exact_answer - r_cb->answer);
	if(tuple[5] == tuple[6]) tuple[7] = 0.5;
	else tuple[7] = tuple[5] < tuple[6] ? 0 : 1;
	tuple[8] = fabs(tuple[5] - tuple[6]);
}

void Print_Tupe(double tuple[]){

	printf("%d,", (int)tuple[0]);
	printf("%lf,", tuple[1]);
	printf("%lf,", tuple[2]);
	printf("%lf,", tuple[3]);
	printf("%e,", tuple[4]);
	printf("%e,", tuple[5]);
	printf("%e,", tuple[6]);
	if(tuple[7] == 0.5) printf("PAR,");
	else if(tuple[7] == 0)	printf("MC,");
	else printf("CB,");
	printf("%e", tuple[8]);
}

void Result_Tuple(double tuple[], RESULT *r){

	tuple[0] = r->answer;
	tuple[1] = r->error;
	tuple[2] = r->lower_bound;
	tuple[3] = r->upper_bound;
	tuple[4] = r->N;
}
void Print_Result(double tuple[]){

	printf("%lf,", tuple[0]);
	printf("%e,", tuple[1]);
	printf("%lf,", tuple[2]);
	printf("%lf,", tuple[3]);
	printf("%e", tuple[4]);
}

