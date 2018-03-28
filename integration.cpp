#include "integration.h"

double unit_N_mc = pow(10.0,6.0); //# of operations a computer can handle in a second
double unit_N_cb = pow(10.0,7.0);

Integration::Integration(int d, double tol){

	this->d = d; //dimension
	this->tol = tol; //tolerance
	C = pow(2.0*R, (double)d);
	answer = 0;
	error = 0;
	upper_bound = 0;
	lower_bound = 0;
	run_time = 0;
	attempt_num = 0;
}

void Integration::Pass_Result(RESULT *r, double N){

	r->N = N;
	r->answer = this->answer;
	r->error = this->error;
	r->upper_bound = this->upper_bound;
	r->lower_bound = this->lower_bound;
	r->run_time = this->run_time;
}

double Integration::Calculate_PI_1(double answer, int d){

	if(answer <= 0) return 0;

	int k = d/2;
	double c = Ball_Volume(R, d, 1);
	return pow(answer/c, 1.0/(double)k);
}

void Integration::Time_Estimator(double N, double unit_N){

	printf("Estimated time: ");

	double time = N/unit_N; //initial unit: seconds
	if(time < 60) {
		printf("%.0lf second(s)\n", time);
		return;
	}

	time = time/60;
	if(time < 60) {
		printf("%.0lf minute(s)\n", time);
		return;
	}

	time = time/60;
	if(time < 24) {
		printf("%.0lf hour(s)\n", time);
		return;
	}

	time = time/24;
	printf("%.0lf day(s)\n", time);
	return;
}

void Monte_Carlo::Generate_Result(RESULT *r, double N){

	printf("\nMonte Carlo starts to generate result...\n");
	start_timer();
	Calculate_Volume(N);
	run_time = elapsed_time(); //the functions return units of seconds
	Pass_Result(r, N);
}

void Cube_based::Generate_Result(RESULT *r, double N){

	printf("\nCube-based starts to generate result...\n");
	start_timer();
	Calculate_Volume(N);
	run_time = elapsed_time(); //the functions return units of seconds
	Pass_Result(r, N);
}

void Monte_Carlo::Calculate_PI(RESULT *r, double N){

	printf("\nMonte Carlo starts to generate result...\n");
	start_timer();
	Calculate_Volume(N);
	run_time = elapsed_time(); //the functions return units of seconds
	Pass_PI_Result(r, N);
}

void Cube_based::Calculate_PI(RESULT *r, double N){

	printf("\nCube-based starts to generate result...\n");
	start_timer();
	Calculate_Volume(N);
	run_time = elapsed_time(); //the functions return units of seconds
	Pass_PI_Result(r, N);
}

void Integration::Pass_PI_Result(RESULT *r, double N){

	r->N = N;
	r->answer = Calculate_PI_1(this->answer, d);

	double pi_u = Calculate_PI_1(this->upper_bound, d);
	double pi_l = Calculate_PI_1(this->lower_bound, d);

	if(pi_u >= pi_l) {

		r->upper_bound = pi_u;
		r->lower_bound = pi_l;
	}
	else{

		r->upper_bound = pi_l;
		r->lower_bound = pi_u;
	}

	r->error = 0.5*(r->upper_bound - r->lower_bound);
	r->run_time = this->run_time;
}

Monte_Carlo::Monte_Carlo(int d, int M, double tol, double t) : Integration(d, tol){

	this->M = M;
	this->t = t;
}

void Monte_Carlo::Generate_Result(RESULT *r){

	printf("\nMonte Carlo starts to generate result...\n");
	start_timer(); //initiate the timer

	double N, N2, p, p_error, z = 2.58;

	//trial run to estimate required N given d and tol
	N = 2*ceil(1.0/(tol*2));
	MC *mc;
	mc = new MC;
	Init_MC(mc, N);
	Thread_Work((void *)mc);
	p = mc->p;
	delete mc;
	p_error = z*sqrt(p*(1-p))/sqrt(N);

	do{
		attempt_num++;

		printf("[Attempt %d] ", attempt_num);
		printf("p = %.3lf(%.3lf), ", p, p_error);

		//choose p as close as to 0.5 to obtain a conservative estimate for N
		if(p + p_error < 0.5) p = p + p_error;
		else if(p - p_error > 0.5) p = p - p_error;
		else p = 0.5;

		N2 = 2*ceil(t*t*C*C*p*(1-p)/(M*tol*tol*2));

		//if last attempt fails, then gradually increase N
		if(N2 > N) N = N2;
		else if(attempt_num > 1 && error > tol) N = 2*ceil(pow(error/tol, 2)*N/2); //error/tol > 1, so N must increase
		else N = 10*N;

		printf("N = %e\n", N);

		Calculate_Volume(N);

	    p = answer/C; //update probability
		p_error = z*sqrt(p*(1-p))/sqrt(N*M);
		//increase confidence level to get a more conservative N

	}while((answer <= 0 || error >= tol) && attempt_num < MAX_ATTEMPT);

	run_time = elapsed_time(); //the functions return units of seconds
	Pass_Result(r, N);
}

void Monte_Carlo::Init_MC(MC *mc, double N){

	mc->d = d;
	mc->N = N;
	mc->p = 0;
}

void Monte_Carlo::Calculate_Volume(double N){

	Time_Estimator(M*N/2, unit_N_mc);

	double t1 = elapsed_time();

	int j, q, q_1, n;
	double average = 0, square_sum = 0, std_error = 0;
	double S;

	//build structures used to pass arguments to thread
	MC *mc[MAX_THREAD];

	for(q = 0; q < MAX_THREAD; q++){

		mc[q] = new MC;
		Init_MC(mc[q], N);
	}

	thread th[MAX_THREAD];

    for(q_1 = 0; q_1 < M; q_1 = q_1 + MAX_THREAD){

    	if(M - q_1 < MAX_THREAD) n = M - q_1;
    	else n =  MAX_THREAD;

		for(q = 0; q < n; q++) th[q] = std::thread(Thread_Work, (void *)mc[q]);

		for(q = 0; q < n; q++) {

			th[q].join();
			S = (mc[q]->p)*C;
			//printf("S: %.3lf, p: %.3lf, C: %.3lf\n", S, mc[q]->p, C);
			j = q_1 + q;
			square_sum = square_sum + ((double)j/(j+1))*(S - average)*(S - average);
			average = average + (S - average)/(j+1);
		}
    }

	for(q = 0; q < MAX_THREAD; q++) delete mc[q];

    if(M > 1) std_error = sqrt(square_sum)/sqrt(M*(M-1)); //standard error
    else if(M == 1) std_error = 0; //only one run => no standard deviation

    answer = average;
    error = t*std_error;
	upper_bound = answer + error;
	lower_bound = max(answer - error, 0.0);

    printf("d = %d, N = %e, M = %d: answer = %.3lf(%.3lf), interval = [%.3lf, %.3lf]\n"
   			, d, N, M, answer, error, lower_bound, upper_bound);

	double t2 = elapsed_time();

	unit_N_mc = (M*N/2)/(t2-t1);
}

void Monte_Carlo::Thread_Work(void *arg){

	MC *mc = (MC *)arg;
	int d = mc->d;
	double N = mc->N;
	double N_1 = N/2;

	double I = 0;
	double x;

    double *v, *av;
    v = new double[d+1];
    av = new double[d+1]; //antithetic sample
    //0th element is used to store distance from the origin (squared)

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-R, R);

	for(double i = 0; i < N_1; i++) //N: number of points
	{
		//generate a random vector
		v[0] = 0, av[0] = 0;

		for(int k = 1; k <= d; k++){ //d: number of dimensions

			  x = dis(gen);

			  //The method of antithetic variates is to use antithetic pairs of samples.
			  v[k] = x;
			  v[0] = v[0] + x*x;
			  av[k] = -x;
			  av[0] = v[0];
		}

		I = I + Is_Inside(v) + Is_Inside(av);

		if(i > 0 && i < 4.3*pow(10.0, 9.0)){ //maximum unsigned long integer

			if((unsigned long)i%100000000 == 0)
			printf("[Random Check (per thread)] i = %e (%%%.0lf, elapsed time: %.0lf)\n"
					, i, 100*i/N_1, elapsed_time());
		}
	}

	delete v;
	delete av;

	mc->p = I/N;
}

double Monte_Carlo::Is_Inside(double v[]){

	if(v[0] > R*R) return 0;
	return 1;
}

Cube_based::Cube_based(int d, double tol) : Integration(d, tol){

	I = 0, O = 0, X = 0, p_1 = 0;
}

void Cube_based::Generate_Result(RESULT *r){

	printf("\nCube-based starts to generate result...\n");
	start_timer(); //initiate the timer

	//expand from [-R -R ... -R] to [R R ... R]
	//gradually increment by h
	//explore each hypercube [x1 x2 ... xd]

	double p, p_error; //p: percentage of X / surface of the hypercube
	double K_1, K, N;

	K_1 = ceil(1.0/(tol*2));
	K = 2*K_1;

	do{
		attempt_num++;
		N = K2N(K, d);

		printf("[Attempt %d] ", attempt_num);
		printf("K = %.0lf, N = %e\n", K, N);

		Calculate_Volume(N);

		if((answer > 0 && error < tol) || attempt_num >= MAX_ATTEMPT) break;

		p = X/(pow(K, (double)d-1)*d*2);
		//percentage of X / surface of the hypercube; increase when K becomes large
		//observation: initial guess underestimates p

		p_error = 1.0/K; //just heuristic

		K_1 = ceil((1 + p_error)*p*d*pow(2.0*R, (double)d)/(tol*2)); //formula verified
		K = 2*K_1; //formula verified

	}while(1);

	run_time = elapsed_time(); //the functions return units of seconds
	Pass_Result(r, N);
}

void Cube_based::Calculate_Volume(double N){

	Time_Estimator(N, unit_N_cb);
	double t1 = elapsed_time();

	double  K = 2*ceil(N2K(N, d)/2);
	double	h = 2.0*R/K; //h: grid size
	double  n;

	//initialize calculation
	I = 0, O = 0, X = 0, p_1 = 0, n = 0;

	if(d == 1){

		CB *cb;
		cb = new CB;
		Init_CB(cb, K);
		Thread_Work((void *)cb);
		Sum_CB(cb, n);
		delete cb;
	}
	else{

		CB *cb[MAX_THREAD];

		int q;
		for(q = 0; q < MAX_THREAD; q++) {

			cb[q] = new CB;
			Init_CB(cb[q], K);
		}

		int check_point;
		check_point = ceil(MAX_N/(MAX_THREAD*pow(K, (double)d - 1)));
		check_point = check_point*MAX_THREAD;

		thread th[MAX_THREAD];

		double t = 0;
		while(t < K){

			for(q = 0; q < MAX_THREAD; q++) {

				cb[q]->k_1 = t + q;
				th[q] = std::thread(Thread_Work, (void *)cb[q]);
			}

			for(q = 0; q < MAX_THREAD; q++) {

				th[q].join();
			}

			if(t > 0 && (int)t%check_point == 0)
				printf("[Random Check (all threads)] k = %.0lf (%%%.0lf, elapsed time: %.0lf)\n"
						, t, 100*t/K, elapsed_time());

			t = t + MAX_THREAD;
		}

		I = 0, O = 0, X = 0, p_1 = 0, n = 0;

		for(q = 0; q < MAX_THREAD; q++) {

			Sum_CB(cb[q], n);
			free(cb[q]);
		}
	} //else d == 1

	if(n != K2N(K, d))
		printf("bug: N = %e is not equal to K^d = %e!\n", n, K2N(K, d));

	if(I + O + X != n) {

		printf("bug: total count should be the same!\n");
		printf("I: %.0lf, O: %.0lf, X: %.0lf, N: %e\n", I, O, X, n);
	}

	double h_volume = pow(h, d);
	double I_volume = I*h_volume;
	double X_volume = X*h_volume;

	answer = I_volume + p_1*h_volume;
	error = 0.5*X_volume;
	upper_bound = I_volume + X_volume;
	lower_bound = I_volume;

	printf("d = %d, N = %e, K = %.0lf: answer = %.3lf(%.3lf), interval = [%.3lf, %.3lf]\n"
			, d, N, K, answer, error, lower_bound, upper_bound);

	double t2 = elapsed_time(); //the functions return units of seconds
	unit_N_cb = N/(t2 - t1);
}

void Cube_based::Sum_CB(CB *cb, double &n){

	I = I + cb->I;
	O = O + cb->O;
	X = X + cb->X;
	p_1 = p_1 + cb->p_1;
	n = n + cb->n;
}

void Cube_based::Init_CB(CB *cb, double K){

	cb->d = d;
	cb->k_1 = 0;
	cb->K_1 = K/2;

	cb->I = 0;
	cb->O = 0;
	cb->X = 0;
	cb->p_1 = 0;
	cb->n = 0;
}

void Cube_based::Thread_Work(void *arg){

	CB *cb = (CB *)arg;
	int d = cb->d;
	double k_1 = cb->k_1, K_1 = cb->K_1;
	double K = 2.0*K_1;

	if(k_1 > K - 1) return;

	int l_1 = 2; //starting level of DFS tree
	if(d == 1) l_1 = 1;

	//construct DFS tree
	double *path, *node;
	double **parent;

	int j;

	path = new double[d+1]; //record frontier in each d, value = 0, 1, ... , K-1
	node = new double[d+2];
	//struct of a node (small cube):
	//node[0]: distance of the nearest edge
	//node[k], k = 1, 2, ... , d: (x1, x2, ... , xk, ... , xd-1, xd)
	//node[d+1]: distance of the farthermost edge

	parent = (double**)malloc(sizeof(double*)*(d+1)); //record parent node in each d
	//in fact parent[d] will not be used, still claim d+1 for convenience

	for(j = 0; j <= d; j++) {

		parent[j] = (double*)malloc(sizeof(double)*(d+2));
	}

	Init_DFS(parent, node, path, cb);

	double N_1 = pow(K, (double)(d - l_1 + 1)); //number of cubes
    double n = 0;
	int l = l_1;

	while(1){

		//return after finish exploring a leaf node
		n = n + Count_Cubes(l, parent, node, path, cb); //explored a node

		if(n >= N_1) {

			cb->n = cb->n + n;
			Free_DFS(parent, node, path, d);

			return;
		}

		while(path[l] >= K-1){ //find the node which is not yet explored

			l = l - 1;
		}

		path[l] = path[l] + 1; //explore other sibling nodes
	}
}

void Cube_based::Init_DFS(double **parent, double *node, double *path, CB *cb){

	int j, k; //j: level, k: kth element in the vector

	for(j = 0; j <= cb->d; j++) {

		path[j] = 0;
	}

	for(k = 0; k <= cb->d+1; k++){

		node[k] = 0;
	}

	for(j = 0; j <= cb->d; j++) {

		for(k = 0; k <= cb->d+1; k++){

			parent[j][k] = 0;
		}

		if(cb->d > 1) Update_Node(parent[j], 1, cb->k_1, cb->d, cb->K_1);
	}
}

void Cube_based::Free_DFS(double **parent, double *node, double *path, int d){

	delete[] path;
	delete[] node;

	int j;
	for(j = 0; j <= d; j++){

		free(parent[j]);
	}

	free(parent);
}

double Cube_based::Count_Cubes(int &l, double **parent, double *node, double *path, CB *cb){

	int d = cb->d;
	double K_1 = cb->K_1;

	Copy_Node(node, parent[l-1], l-1, d); //inherit parent

	//try another value
	Update_Node(node, l, path[l], d, K_1);
	if(node[0] >= node[d+1]) printf("bug: farthermost edge's distance should be longer!\n");

	if(l < d){ //not yet reach the leaf node

		Copy_Node(parent[l], node, l, d);
		l = l + 1;
		path[l] = 0; //explore the left child first
		return Count_Cubes(l, parent, node, path, cb);
	}
	else{ //l == d, has reached the leaf node

		//classify nodes
		double c = Is_Inside(node, d);

		if(c == 1) cb->I = cb->I + 1;
		else if(c == 0) cb->O = cb->O + 1;
		else {

			cb->X = cb->X + 1;
			cb->p_1 = cb->p_1 + c;
		}

		return 1;
	}
}

double Cube_based::Is_Inside(double node[], int d){

	//classify nodes
	if(node[0] > R*R) return 0; //outside
	else if(node[d+1] <= R*R) return 1; //inside
	else return (R - sqrt(node[0]))/(sqrt(node[d+1]) - sqrt(node[0])); //intersect interface
}

void Cube_based::Copy_Node(double dest[], double src[], int j, int d){

	for(int k = 0; k <= j; k++) dest[k] = src[k];
	dest[d+1] = src[d+1];
}

void Cube_based::Update_Node(double node[], int k, double x, int d, double K_1){

	double x_1, h = (double)R/K_1;

	node[k] = x;
	x_1 = -R + x*h;

	if(x < K_1){ //K_1: middle point => x_1 < 0

		node[0] = node[0] + (x_1 + h)*(x_1 + h);
		node[d+1] = node[d+1] + x_1*x_1;
	}
	else{

		node[0] = node[0] + x_1*x_1;
		node[d+1] = node[d+1] + (x_1 + h)*(x_1 + h);
	}
}
