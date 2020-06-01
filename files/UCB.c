#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define min(a,b,c) (a<b?(a<c?a:c):(b<c?b:c))
#define EPSILON 1e-6
#define PI 3.141592653589793

struct double_pair
{
	double zero;
	double one;
};

struct double_int
{
    double value;
    int index;
};

const int K = 14;
const int M = 3;
const int T = 1e6;
const int MC_iterations = 10;
const double p = 0.49;
const double q = 0.45;
double variance_proxy;
double MAX_theta;

// state reward function
int f(int x)
{
	return(2*x-1);
}

double g(int t)
{
	return(log(t) + 3*log(log(t)));
}

//  PF-eigenvalue
double rho(double theta)
{
	double square = pow(q*exp(theta)-(1-p)*exp(-theta),2);
	return(((1-p)*exp(-theta) + q*exp(theta) + sqrt(square + 4*p*(1-q)))/2);
}

// log-PF-eigenvalue
double log_rho(double theta)
{
    	return(log(rho(theta)));
}

// log-PF-eigenvalue derivative
double log_rho_prime(double theta)
{
	double square = pow(q*exp(theta)-(1-p)*exp(-theta),2);
    	return((q*exp(theta)-(1-p)*exp(-theta))/sqrt(square + 4*p*(1-q)));
}

// inverse-log-PF-eigenvalue derivative
double inverse_log_rho_prime(double mu)
{
	double mu_square = pow(mu,2);
	return(log((sqrt(p*(1-q))*mu + sqrt(p*(1-q)*mu_square + q*(1-p)*(1-mu_square)))/(q*sqrt(1-mu_square))));
}

// theta transition probability matrix
struct double_pair theta_transitions(double theta)
{
	double r = rho(theta);
	double p_theta = 1 - (1-p)*exp(-theta)/r;
	double q_theta = q*exp(theta)/r;
	struct double_pair P_theta = {p_theta, q_theta};
	return(P_theta);
}

// KL-divergence rate
double KL(double theta1, double theta2)
{
	return(log_rho(theta2) - log_rho(theta1) - log_rho_prime(theta1)*(theta2 - theta1));
}

// invert the KL-divergence using the bisection method
double inverse_KL(double theta1, double y)
{
	double a = theta1, b = MAX_theta;

	if (a >= b)
	{
		printf("This shouldn't happen! Essentially theta1 infinity was given in the KL.\n");
		return(a);
	};

	double c = a;
	while ((b-a) >= EPSILON)
	{
		c = (a+b)/2;
		double div = KL(theta1,c);
		if (div < y)
			a = c;
		else if (div > y)
			b = c;
		else
			break;	
	};
	return(c);
}


double uniform()
{
	return((double) rand() / (double) ((unsigned)RAND_MAX + 1));
}

struct double_pair standard_normal_pair()
{
	double u1 = uniform(), u2 = uniform();
	double r = sqrt(-2*log(u1)), theta = 2*PI*u2;
	double z1 = r*cos(theta), z2 = r*sin(theta);
	struct double_pair z = {z1, z2};
	return(z);
}

// comparison functions for qsort
int cmp_double_descending(const void * a, const void * b)
{
	return (*(double*)a > *(double*)b) ? -1 : (*(double*)a < *(double*)b) ? 1:0 ;
}

int cmp_double_int_descending(const void *a, const void *b)
{
    	struct double_int *a1 = (struct double_int *) a;
    	struct double_int *b1 = (struct double_int *) b;
    	return ((*a1).value > (*b1).value) ? -1 : ((*a1).value < (*b1).value) ? 1:0;
}

// UCB index
double UCB(int s, int n, int t)
{
	double avg = ((double) s)/n;
	return(avg + sqrt(4*variance_proxy*log(t)/n));
}

// KL-UCB index
double KL_UCB(int s, int n, int t)
{
	double avg = ((double) s)/n;
	if (avg >= 1)
		return(1);
	else if (avg <= -1)
		avg = -nextafter(1,0);

	double theta1 = inverse_log_rho_prime(avg);
	double theta2 = inverse_KL(theta1, g(t)/n);
	double KL_UCB = log_rho_prime(theta2);
	if (KL_UCB < avg)
		printf("KL-UCB smaller than avg! This shouldn't happen.%lf, %lf, %lf, %lf, %lf\n", avg, KL_UCB, theta2, theta1, g(t)/n);
	return(KL_UCB);
}

void index_policy(double Ps[K][2], double best_M_rewards, double regrets[T], double (*index_calculator)(int, int, int))
{
	int states[K] = {0};
	int times_pulled[K] = {0};
	int sums[K] = {0};
	struct double_int indices[K];

	// initialize the K arms, by pulling each one M times
	for (int a = 0; a < K; a++)
	{
		for (int n = 0; n < M; n++)
		{
			int current_state = states[a];
			int new_state = (uniform() < Ps[a][current_state]) ? 1 : 0;
			states[a] = new_state;
			times_pulled[a] += 1;
			int reward = f(new_state);
			sums[a] += reward;
		};
		regrets[K-1] += best_M_rewards - sums[a];
		// printf("%d, %d, %d\n", a, sums[a], times_pulled[a]);
	};

	for (int t = K; t < T; t++)
	{
		// calculate indices for all K arms
		for (int a = 0; a < K; a++)
		{
			indices[a].value = index_calculator(sums[a], times_pulled[a], t);
			indices[a].index = a;
		};
		// sort arms based on their indices
		qsort(indices, K, sizeof(struct double_int), cmp_double_int_descending);
		// play top-M arms
		//for (int i = 0; i < K; i++)
		//{
		// 	printf("%d, %lf || ", indices[i].index, indices[i].value);
		// };
		// printf("\n");
		double round_reward = 0;
		for (int i = 0; i < M; i++)
		{
			int a = indices[i].index;
			// printf("%d ", a);
			int current_state = states[a];
			int new_state = (uniform() < Ps[a][current_state]) ? 1 : 0;
			states[a] = new_state;
			times_pulled[a] += 1;
			int reward = f(new_state);
			sums[a] += reward;
			round_reward += reward;
		};
		regrets[t] += best_M_rewards - round_reward;
		// printf("\n%lf\n", regrets[t]);
	};

}

int main(void)
{

	double thetas[K];
	double mus[K];
	double Ps[K][2];
	
	srand(100);
	for (int a = 0; a < K; a += 2)
	{
		struct double_pair z = standard_normal_pair();
		thetas[a] = z.zero/4;
		thetas[a+1] = z.one/4;
	}
	qsort(thetas, K, sizeof(double), cmp_double_descending);
	double m = 1;
	for (int a = 0; a < K; a++)
	{
	     	double theta = thetas[a];
	      	mus[a] = log_rho_prime(theta);
		struct double_pair P_theta =  theta_transitions(theta);
		Ps[a][0] = P_theta.zero;
		Ps[a][1] = P_theta.one;
		m = min(m, Ps[a][0], 1-Ps[a][1]);
		// printf("%lf, %lf, %lf, %lf\n", thetas[a], mus[a], Ps[a][0], Ps[a][1]);
	};


	// variance_proxy = 1/pow(m,2);
	variance_proxy = 1;
	
	double best_M_rewards = 0;
	for (int a = 0; a < M; a++)
	{
		best_M_rewards += mus[a];
	};

	double optimal_constant = 0;
	for (int b = M; b < K; b++)
	{
		optimal_constant += (mus[M-1] - mus[b])/KL(thetas[b], thetas[M-1]);
	};
	
	// maximum value of theta is the one that makes mu almost one
	MAX_theta = inverse_log_rho_prime(nextafter(1,0));

	srand(time(NULL));
	
	static double regrets[T] = {0};
	for (int i = 0; i < MC_iterations; i++)
	{
		// printf("Round %d\n", i);
		// index_policy(Ps, best_M_rewards, regrets, UCB);
		index_policy(Ps, best_M_rewards, regrets, KL_UCB);
	};
	static double cum_regrets[T] = {0};
	cum_regrets[K-1] = regrets[K-1]/MC_iterations;
	for (int t = K; t < T; t++)
	{
		cum_regrets[t] = cum_regrets[t-1] + regrets[t]/MC_iterations;
	};
	printf("%lf, %lf, %lf\n", regrets[K-1], regrets[T-2], regrets[T-1]);
	printf("%lf, %lf, %lf\n", cum_regrets[K-1], cum_regrets[T-2], cum_regrets[T-1]);
	printf("%lf, %lf, %lf\n", optimal_constant, log(T), optimal_constant*log(T));

	return(0);
}

