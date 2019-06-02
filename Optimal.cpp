// ComputerSimulation_report.cpp : 此檔案包含 'main' 函式。程式會於該處開始執行及結束執行。
//

//#include "pch.h"
#include <iostream>
#include <random>
#include <cmath>
#include <time.h>
#include <iostream>
#include <math.h>
//#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
//#include <stdio.h>
#include <fstream>

//using namespace cv;
using namespace std;

struct member_nodes {
	int UI = 0;
	int US = 0;
	int PV = 0;
	int DV = 0;
	double fraction_UI = 0.0;
	double fraction_US = 0.0;
	double fraction_PV = 0.0;
	double fraction_DV = 0.0;
};

struct prob_every_state {
	double p_ui = 0; //t+1
	double p_us = 0; //t+1
	double p_pv = 0; //t+1
	double p_dv = 0; //t+1
};

struct parametrers_every_state {
	double p_ui = 0; //t+1
	double p_us = 0; //t+1
	double p_pv = 0; //t+1
	double p_dv = 0; //t+1
	double q = 0;
	double m = 0;
};

double lambdaBS = 200; //density*area
double length = 500;
int run_time = 13;
int patch_insert_time = 4;
int number_of_nodes = 1000;
int US_state = 1, UI_state = 2, PV_state = 3, DV_state = 4;
int Average_degree_patch = 4; // 5比較fit
int t = 0;
// patch from central
double m = 0.0;
double delta_v = 0.2; //0.2 0.6 1 for fig4.
//pv->dv
double delta_p = 0.0;
double beta_v = 0.6; //0.6比較吻合
//us ui -> pv
double beta_p = 0.0;
fstream file;

double prob_not_infected(int **adjacency_matrix_patch, prob_every_state* prob_each_node, double beta_v, int i, int BS_number);
double prob_not_patched(int **adjacency_matrix_patch, prob_every_state* prob_each_node, double beta_p, int i, int BS_number);
double prob_UI(prob_every_state* prob_each_node, double q, double r, double delta_v, double m, int i);
double prob_US(prob_every_state* prob_each_node, double q, double r, double delta_v, double m, int i);
double prob_PV(prob_every_state* prob_each_node, double q, double r, double delta_p, double m, int i);
double prob_DV(prob_every_state* prob_each_node, double delta_p, int i);

double lambdaUS(double**lambda, double***prob_detail, int i, int t);
double lambdaUI(double**lambda, double***prob_detail, int**adjacency_matrix_patch, int i, int t, double delta_v, double beta_v);
double lambdaPV(double**lambda, double***prob_detail, int**adjacency_matrix_patch, int i, int t);
double lambdaDV(double**lambda, double***prob_detail, int i);
double lambdaQ(double**lambda, double***prob_detail, int i, int t);
double lambdaR(double**lambda, double***prob_detail, int i, int t, double delta_v);
double M(double**lambda, double***prob_detail, int i, int t, double delta_v);
double DeltaP(double**lambda, double***prob_detail, int i, int t);

int Determine_state_for_time_t(int* state_matrix, prob_every_state* prob_each_node, int i, int BS_number);
void upgrading_state_prob(prob_every_state* temp, prob_every_state* prob_each_node, int BS_number);
void count_members_time_t(int* state_matrix, member_nodes* total_nodes_distribution, int t, int BS_number);

int main()
{
	//圖
	//Mat img(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
	//撒點
	default_random_engine generatorBS(time(NULL));
	poisson_distribution<int> distributionBS(lambdaBS); //lambdaBS = density*Area

	default_random_engine putBS(time(NULL));
	uniform_real_distribution<double> unif(-1 * length, length);

	int BS_number = distributionBS(generatorBS);
	cout << BS_number << endl;

	//double *BS_x, *BS_y;
	//BS_x = new double[BS_number];
	//BS_y = new double[BS_number];
	//position matrix
	//Point* BS = new Point[BS_number];
	/**
	for (int i = 0; i < BS_number; i++)
	{
		BS[i].x = unif(putBS) + length;
		BS[i].y = unif(putBS) + length;
		//cout << BS[i].x << "　"<< BS[i].y << endl;
		circle(img, BS[i], 3, Scalar(0, 0, 0), 0);
	}
	*/
	//adjacency matrix
	int** adjacency_matrix_patch = new int*[BS_number];
	for (int i = 0; i < BS_number; i++) {
		adjacency_matrix_patch[i] = new int[BS_number];
	}

	for (int i = 0; i < BS_number; ++i) {
		for (int j = 0; j < BS_number; ++j) {
			adjacency_matrix_patch[i][j] = 0;
		}
	}
	//機率亂數
	srand((unsigned)time(NULL));
	//add edges via average degree
	for (int i = 0; i < (BS_number*Average_degree_patch) / 2; i++) {
		int num1 = rand() % BS_number;
		int num2 = rand() % BS_number;
		adjacency_matrix_patch[num1][num2] = 1;
		//line(img, BS[num1], BS[num2], Scalar(0, 0, 0));
	}
	//prob in 4 state for each node at time t 和準備一個RE
	prob_every_state* prob_each_node = new prob_every_state[BS_number];
	for (int i = 0; i < BS_number; ++i) {
		prob_each_node[i].p_us = 1.0;
		prob_each_node[i].p_ui = 0;
		prob_each_node[i].p_pv = 0;
		prob_each_node[i].p_dv = 0;
	}

	prob_every_state* prob_each_node_RE = new prob_every_state[BS_number];
	for (int i = 0; i < BS_number; ++i) {
		prob_each_node_RE[i].p_us = 1.0;
		prob_each_node_RE[i].p_ui = 0;
		prob_each_node_RE[i].p_pv = 0;
		prob_each_node_RE[i].p_dv = 0;
	}

	//存每個t時間的所有參數
	double ***prob_detail = new double**[run_time];
	for (int i = 0; i < run_time; ++i) {
		prob_detail[i] = new double*[BS_number];
		for (int j = 0; j < BS_number; ++j) {
			prob_detail[i][j] = new double[8];//ui us pv dv q r m delta_p
		}
	}
	for (int i = 0; i < run_time; i++)
		for (int j = 0; j < BS_number; j++)
			for (int k = 0; k < 8; k++)
				prob_detail[i][j][k] = 0.0;
	//state at time t 和 準備一個重新執行matrix
	int* state_matrix = new int[BS_number];
	for (int i = 0; i < BS_number; ++i) {
		state_matrix[i] = US_state;
	}

	int* state_matrix_RE = new int[BS_number];
	for (int i = 0; i < BS_number; ++i) {
		state_matrix_RE[i] = US_state;
	}

	//select n = 10 node infected
	int n = 1;
	for (int i = 0; i < n; ++i) {
		int infected_at_time_0 = rand() % BS_number;
		prob_each_node[infected_at_time_0].p_us = 0.8;
		prob_each_node[infected_at_time_0].p_ui = 0.2;
		state_matrix[infected_at_time_0] = UI_state;
		//重新執行matrix
		prob_each_node_RE[infected_at_time_0].p_us = 0.8;
		prob_each_node_RE[infected_at_time_0].p_ui = 0.2;
		state_matrix_RE[infected_at_time_0] = UI_state;
		for (int k = 0; k < BS_number; ++k) {
			adjacency_matrix_patch[k][infected_at_time_0] = 1;
		}
	}
	//total node count 和準備一個重新執行
	member_nodes* total_nodes_distribution = new member_nodes[run_time];
	member_nodes* total_nodes_distribution_RE = new member_nodes[run_time];
	//寫入txt
	//fstream file;
	//file.open("C:\\Users\\User\\Desktop\\Matlab Code\\result.txt", ios::out);
	//fstream file1;
	//file1.open("C:\\Users\\User\\Desktop\\Matlab Code\\result_figure3(c)_centralize.txt", ios::out);
	//simulation
	while (t < run_time) {
		//determine whether changing
		for (int i = 0; i < BS_number; ++i) {
			state_matrix[i] = Determine_state_for_time_t(state_matrix, prob_each_node, i, BS_number);
		}
		//count members
		count_members_time_t(state_matrix, total_nodes_distribution, t, BS_number);
		//cout << "比率 at time " << t << " : " << total_nodes_distribution[t].fraction_US << " " << total_nodes_distribution[t].fraction_UI << " " << total_nodes_distribution[t].fraction_PV << " " << total_nodes_distribution[t].fraction_DV << endl;
		cout << "比率 at time " << t << " : " << total_nodes_distribution[t].fraction_UI << endl;
		//calculate probability for t+1
		prob_every_state* temp = new prob_every_state[BS_number];
		// t = patch_input_time
		//if (t == patch_insert_time) { //mixed
			//m = 0.03;
			//delta_p = 0.4;
			//beta_p = 0.6;
		//}
		//if (t == patch_insert_time) { //centralize
			//m = 0.01;
			//delta_p = 0.0;
			//beta_p = 0.0;
		//}
		/*
		if (t == patch_insert_time) { //decentralize
			m = 0.0;
			delta_p = 0.4;
			beta_p = 0.6;
			int random_patch_at_time_c = rand() % BS_number; // one node patch
			cout << random_patch_at_time_c << endl;
			prob_each_node[random_patch_at_time_c].p_us = 0;
			prob_each_node[random_patch_at_time_c].p_ui = 0;
			prob_each_node[random_patch_at_time_c].p_pv = 1;
			prob_each_node[random_patch_at_time_c].p_dv = 0;
		}
		*/
		for (int i = 0; i < BS_number; ++i) {
			if (t >= patch_insert_time) {
				m = (double)(0.25 - 0.0) * rand() / (RAND_MAX + 1.0) + 0.0;
				delta_p = (1.0 - 0.0) * rand() / (RAND_MAX + 1.0) + 0.0;
			}
			
			prob_detail[t][i][6] = m;
			prob_detail[t][i][7] = delta_p;
			double q = prob_not_infected(adjacency_matrix_patch, prob_each_node, beta_v, i, BS_number);
			double r = prob_not_patched(adjacency_matrix_patch, prob_each_node, beta_p, i, BS_number);
			prob_detail[t][i][4] = q;
			prob_detail[t][i][5] = r;
			//cout <<  q  << " " << r << endl;
			temp[i].p_us = prob_US(prob_each_node, q, r, delta_v, m, i);
			temp[i].p_ui = prob_UI(prob_each_node, q, r, delta_v, m, i);
			temp[i].p_pv = prob_PV(prob_each_node, q, r, delta_p, m, i);
			temp[i].p_dv = prob_DV(prob_each_node, delta_p, i);
		}
		//store prob at time t
		for (int i = 0; i < BS_number; ++i) {
			prob_detail[t][i][0] = prob_each_node[i].p_us;
			prob_detail[t][i][1] = prob_each_node[i].p_ui;
			prob_detail[t][i][2] = prob_each_node[i].p_pv;
			prob_detail[t][i][3] = prob_each_node[i].p_dv;
		}
		//update probability for t+1
		upgrading_state_prob(temp, prob_each_node, BS_number);
		cout << "機率 at time " << t + 1 << " : " << prob_each_node[0].p_us << " " << prob_each_node[0].p_ui << " " << prob_each_node[0].p_pv << " " << prob_each_node[0].p_dv << endl;

		++t;
	}
	/*
	//寫入 figure 3(b)
	for (int i = 0; i < run_time; ++i) {
		file << total_nodes_distribution[i].fraction_US << " " << total_nodes_distribution[i].fraction_UI << " " << total_nodes_distribution[i].fraction_PV << " " << total_nodes_distribution[i].fraction_DV << endl;
		//file1 << total_nodes_distribution[i].fraction_UI << endl;
	}
	file.close();
	//imshow("Poisson geometric graph", img);
	//waitKey();
	*/

	//Lagrangian Optimization Algo
	double** lambda = new double*[BS_number];//us ui pv dv q r = 0 at time T
	for (int i = 0; i < BS_number; ++i) {
		lambda[i] = new double [6];
	}
	for (int i = 0; i < BS_number; ++i) {
		for (int j = 0; j < 6; ++j) {
			lambda[i][j] = 0;
		}
	}
	for (int i = 0; i < BS_number; ++i) {
		lambda[i][1] = 1; //ui = 1 at time T
	}
	t = run_time - 1;
	while (t > 0) {//19 18 17 16....
		for (int i = 0; i < BS_number; ++i) {
			double lambdaUS_t_1 = 0.0;
			double lambdaUI_t_1 = 0.0;
			double lambdaPV_t_1 = 0.0;
			double lambdaDV_t_1 = 0.0;
			double lambdaQ_t_1 = 0.0;
			double lambdaR_t_1 = 0.0;
			double optimal_m = 0.0;
			double optimal_delta_p = 0.0;

			lambdaUS_t_1 = lambdaUS(lambda, prob_detail, i, t);//count lambdaUS at time t-1
			lambdaDV_t_1 = lambdaDV(lambda, prob_detail, i);//count lambdaR at time t-1
			lambdaQ_t_1 = lambdaQ(lambda, prob_detail, i, t);//count lamdbaUI at time t-1
			lambdaR_t_1 = lambdaR(lambda, prob_detail, i, t, delta_v);//count lambdaPV at time t-1
			lambdaUI_t_1 = lambdaUI(lambda, prob_detail, adjacency_matrix_patch, i, t, delta_v, beta_v);//count lambdaDV at time t-1
			lambdaPV_t_1 = lambdaPV(lambda, prob_detail, adjacency_matrix_patch, i, t);//count lambdaQ at time t-1
			
			optimal_m = M(lambda, prob_detail, i, t, delta_v);
			optimal_delta_p = DeltaP(lambda, prob_detail, i, t);
			prob_detail[t-1][i][6] = optimal_m;
			prob_detail[t-1][i][7] = optimal_delta_p;
			if (i == 0) {
				cout << "t:" << t << endl;
				cout << lambdaUS_t_1 << endl;
				cout << lambdaUI_t_1 << endl;
				cout << lambdaPV_t_1 << endl;
				cout << lambdaDV_t_1 << endl;
				cout << lambdaQ_t_1 << endl;
				cout << lambdaR_t_1 << endl;
				cout << optimal_m << "　" << optimal_delta_p << endl;
			}
			//count lambdaM and M at time t-1
			//count lambdadelta and deltaP at time t-1

			//更新成t-1
			lambda[i][0] = lambdaUS_t_1;
			lambda[i][1] = lambdaUI_t_1;
			lambda[i][2] = lambdaPV_t_1;
			lambda[i][3] = lambdaDV_t_1;
			lambda[i][4] = lambdaQ_t_1;
			lambda[i][5] = lambdaR_t_1;
		}
		t--;
	}
	
	//重新執行
	t = 0;
	while (t < run_time) {
		//determine whether changing
		for (int i = 0; i < BS_number; ++i) {
			state_matrix_RE[i] = Determine_state_for_time_t(state_matrix_RE, prob_each_node_RE, i, BS_number);
		}
		//count members
		count_members_time_t(state_matrix_RE, total_nodes_distribution_RE, t, BS_number);
		//cout << "比率 at time " << t << " : " << total_nodes_distribution_RE[t].fraction_US << " " << total_nodes_distribution_RE[t].fraction_UI << " " << total_nodes_distribution_RE[t].fraction_PV << " " << total_nodes_distribution_RE[t].fraction_DV << endl;
		cout << "比率 at time " << t << " : " << total_nodes_distribution_RE[t].fraction_UI << endl;
		//calculate probability for t+1
		prob_every_state* temp1 = new prob_every_state[BS_number];
		for (int i = 0; i < BS_number; ++i) {
			if (t >= patch_insert_time) {
				m       = prob_detail[t][i][6];
				delta_p = prob_detail[t][i][7];
			}
			double q = prob_not_infected(adjacency_matrix_patch, prob_each_node_RE, beta_v, i, BS_number);
			double r = prob_not_patched(adjacency_matrix_patch, prob_each_node_RE, beta_p, i, BS_number);
			temp1[i].p_us = prob_US(prob_each_node, q, r, delta_v, m, i);
			temp1[i].p_ui = prob_UI(prob_each_node, q, r, delta_v, m, i);
			temp1[i].p_pv = prob_PV(prob_each_node, q, r, delta_p, m, i);
			temp1[i].p_dv = prob_DV(prob_each_node, delta_p, i);
		}
		//update probability for t+1
		upgrading_state_prob(temp1, prob_each_node_RE, BS_number);
		cout << "機率 at time " << t + 1 << " : " << prob_each_node_RE[0].p_us << " " << prob_each_node_RE[0].p_ui << " " << prob_each_node_RE[0].p_pv << " " << prob_each_node_RE[0].p_dv << endl;
		++t;
	}
	
	return 0;
}

double prob_not_infected(int **adjacency_matrix_patch, prob_every_state* prob_each_node, double beta_v, int i, int BS_number) {
	double result = 1.0;
	for (int j = 0; j < BS_number; ++j) {
		double temp = 0.0;
		temp = 1 - adjacency_matrix_patch[i][j] * prob_each_node[j].p_ui*beta_v;
		result *= temp;

	}
	return result;
};

double prob_not_patched(int **adjacency_matrix_patch, prob_every_state* prob_each_node, double beta_p, int i, int BS_number) {
	double result = 1.0;
	for (int j = 0; j < BS_number; ++j) {
		double temp;
		temp = 1 - adjacency_matrix_patch[i][j] * prob_each_node[j].p_pv*beta_p;
		result *= temp;
	}
	return result;
};

double prob_US(prob_every_state* prob_each_node, double q, double r, double delta_v, double m, int i) {
	double result = 1.0;
	result = prob_each_node[i].p_us*(1 - m)*r*q + prob_each_node[i].p_ui*(1 - m)*r*delta_v;
	return result;
};

double prob_UI(prob_every_state* prob_each_node, double q, double r, double delta_v, double m, int i) {
	double result = 1.0;
	result = prob_each_node[i].p_us*(1 - m)*r*(1 - q) + prob_each_node[i].p_ui*(1 - m)*r*(1 - delta_v);
	return result;
};

double prob_PV(prob_every_state* prob_each_node, double q, double r, double delta_p, double m, int i) {
	double result = 1.0;
	result = prob_each_node[i].p_us*(1 - (1 - m)*r) + prob_each_node[i].p_ui*(1 - (1 - m)*r) + prob_each_node[i].p_pv*(1 - delta_p);
	return result;
};

double prob_DV(prob_every_state* prob_each_node, double delta_p, int i) {
	double result = 1.0;
	result = prob_each_node[i].p_dv + prob_each_node[i].p_pv*delta_p;
	return result;
};
//Lagrangian optimize
double lambdaUS(double**lambda, double***prob_detail, int i, int t) {
	double temp = 0.0;
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	a = -1 * lambda[i][0] * (1 - prob_detail[t - 1][i][6])*prob_detail[t - 1][i][5] * prob_detail[t - 1][i][4];
	b = -1 * lambda[i][1] * (1 - prob_detail[t - 1][i][6])*prob_detail[t - 1][i][5] * (1 - prob_detail[t - 1][i][4]);
	c = -1 * lambda[i][2] * (1 - (1 - prob_detail[t-1][i][6])*prob_detail[t - 1][i][5]);
	temp = -1 * (a + b + c);
	return temp;
};

double lambdaUI(double**lambda, double***prob_detail, int**adjacency_matrix_patch, int i, int t, double delta_v, double beta_v) {
	double temp = 0.0;
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;
	double e = 0.0;
	a = -1 * lambda[i][0] * (1 - prob_detail[t - 1][i][6])*prob_detail[t - 1][i][5] * delta_v;
	b = -1 * lambda[i][1] * (1 - prob_detail[t - 1][i][6])*prob_detail[t - 1][i][5] * (1 - delta_v);
	c = -1 * lambda[i][2] * ((1 - (prob_detail[t - 1][i][6]))*prob_detail[t - 1][i][5]);
	//積分
	//d = prob_detail[t][i][4]/(1- adjacency_matrix_patch)
	return temp;
};

double lambdaPV(double**lambda, double***prob_detail, int**adjacency_matrix_patch, int i, int t) {
	double temp = 0.0;
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;
	a = -1 * lambda[i][2] * (1 - prob_detail[t - 1][i][7]);
	b = -1 * lambda[i][3] * prob_detail[t - 1][i][7];
	return temp;
};

double lambdaDV(double**lambda, double***prob_detail, int i) {
	double temp = 0.0;
	temp = lambda[i][3];
	return 0.0;
};

double lambdaQ (double**lambda, double***prob_detail, int i, int t) {
	double temp = 0.0;
	double a = 0.0;
	double b = 0.0;
	a = -1 * lambda[i][0] * prob_detail[t - 1][i][0] * (1 - prob_detail[t - 1][i][6])*prob_detail[t - 1][i][5];
	b = -1 * lambda[i][1] * prob_detail[t - 1][i][0] * (1 - prob_detail[t - 1][i][6])*prob_detail[t - 1][i][5];
	temp = -1 * (a + b);
	return temp;
};

double lambdaR (double**lambda, double***prob_detail, int i, int t, double delta_v) {
	double temp = 0.0;
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;
	double e = 0.0;
	double f = 0.0;
	a = -1 * lambda[i][0] * prob_detail[t - 1][i][0] * (1 - prob_detail[t - 1][i][6])*prob_detail[t - 1][i][4];
	b = -1 * lambda[i][0] * prob_detail[t - 1][i][1] * (1 - prob_detail[t - 1][i][6])*delta_v;
	c = -1 * lambda[i][1] * prob_detail[t - 1][i][0] * (1 - prob_detail[t - 1][i][6])*(1 - prob_detail[t - 1][i][4]);
	d = -1 * lambda[i][1] * prob_detail[t - 1][i][1] * (1 - prob_detail[t - 1][i][6])*(1 - delta_v);
	e = lambda[i][2] * prob_detail[t - 1][i][0] * (1 - prob_detail[t - 1][i][6]);
	f = lambda[i][2] * prob_detail[t - 1][i][1] * (1 - prob_detail[t - 1][i][6]);
	temp = -1 * (a + b + c + d + e + f);
	return temp;
};

double M(double**lambda, double***prob_detail, int i, int t, double delta_v) {
	double temp = 0.0;
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;
	double e = 0.0;
	double f = 0.0;
	double epsilon = 0.0001;
	a = lambda[i][0] * prob_detail[t - 1][i][0] * prob_detail[t-1][i][5] * prob_detail[t - 1][i][4];
	b = lambda[i][0] * prob_detail[t - 1][i][1] * prob_detail[t - 1][i][5] * delta_v;
	c = lambda[i][1] * prob_detail[t - 1][i][0] * prob_detail[t - 1][i][5] * (1 - prob_detail[t - 1][i][4]);
	d = lambda[i][1] * prob_detail[t - 1][i][1] * prob_detail[t - 1][i][5] * (1 - delta_v);
	e = -1 * lambda[i][2] * prob_detail[t - 1][i][0] * prob_detail[t - 1][i][5];
	f = -1 * lambda[i][2] * prob_detail[t - 1][i][1] * prob_detail[t - 1][i][5];
	//epsilon 很小
	temp = -1 * (a + b + c + d + e + f) - epsilon / (prob_detail[t - 1][i][6] - 0) + epsilon / (1 - prob_detail[t - 1][i][6]);//lambdaM
	double g = 0.0;
	double h = 0.0;
	double result = 0.0;
	g = (temp + (a + b + c + d + e + f)-2*epsilon + sqrt(pow(temp + (a + b + c + d + e + f), 2) + 4 * pow(epsilon, 2))) / (2 * (temp + (a + b + c + d + e + f)));
	h = (temp + (a + b + c + d + e + f)-2*epsilon - sqrt(pow(temp + (a + b + c + d + e + f), 2) + 4 * pow(epsilon, 2))) / (2 * (temp + (a + b + c + d + e + f)));
	
	//double test = sqrt(pow(temp + (a + b + c + d + e + f), 2) + 4 * pow(epsilon, 2));
	//double test1 = temp + (a + b + c + d + e + f);
	//cout << prob_detail[t][i][6] << " " << g << " " << h <<  endl;

	if (g < 0 || g > 1) return h;
	if (h < 0 || h > 1) return g;
	if (g > h) result = g;
	else result = h;
	return result;
};

double DeltaP(double**lambda, double***prob_detail, int i, int t) {
	double temp = 0.0;
	double a = 0.0;
	double b = 0.0;
	double epsilon = 0.0001;
	a = lambda[i][2] * prob_detail[t - 1][i][2];
	b = -1 * lambda[i][3] * prob_detail[t - 1][i][2];
	temp = -1 * (a + b) - epsilon / (prob_detail[t - 1][i][7] - 0) + epsilon / (1 - prob_detail[t - 1][i][7]);
	double c = 0.0;
	double d = 0.0;
	double result = 0.0;
	c = (temp + (a + b) - 2 * epsilon + sqrt(pow(temp + (a + b), 2) + 4 * pow(epsilon, 2))) / (2 * (temp + (a + b)));
	d = (temp + (a + b) - 2 * epsilon - sqrt(pow(temp + (a + b), 2) + 4 * pow(epsilon, 2))) / (2 * (temp + (a + b)));
	
	//double test = sqrt(pow(temp + (a + b), 2) + 4 * pow(epsilon, 2));
	//double test1 = temp + (a + b);
	//cout << prob_detail[t][i][7] <<" " << c << " " << d  << endl;
	
	if (c < 0 || c > 1) return d;
	if (d < 0 || d > 1) return c;
	if (c > d) result = c;
	else result = d;
	return result;
};
int Determine_state_for_time_t(int* state_matrix, prob_every_state* prob_each_node, int i, int BS_number) {
	//US->US or UI or PV 
	//UI->US or US or PV
	//PV -> PV or DV
	//DV -> DV
	double r = (double)rand() / (RAND_MAX);
	//cout << r << endl;
	double total_count = 0.0;
	double us_prob, ui_prob, pv_prob, dv_prob;
	double reset_prob = 1.0;
	switch (state_matrix[i]) {
	case 1:
		total_count = prob_each_node[i].p_us + prob_each_node[i].p_ui + prob_each_node[i].p_pv;
		us_prob = reset_prob / total_count * prob_each_node[i].p_us;
		ui_prob = reset_prob / total_count * prob_each_node[i].p_ui;
		pv_prob = reset_prob / total_count * prob_each_node[i].p_pv;
		if (r <= us_prob) return US_state;
		r -= us_prob;
		if (r <= ui_prob) return UI_state;
		r -= ui_prob;
		if (r <= pv_prob) return PV_state;
		return -1;
		break;
	case 2:
		total_count = prob_each_node[i].p_us + prob_each_node[i].p_ui + prob_each_node[i].p_pv;
		us_prob = reset_prob / total_count * prob_each_node[i].p_us;
		ui_prob = reset_prob / total_count * prob_each_node[i].p_ui;
		pv_prob = reset_prob / total_count * prob_each_node[i].p_pv;
		if (r <= us_prob) return US_state;
		r -= us_prob;
		if (r <= ui_prob) return UI_state;
		r -= ui_prob;
		if (r <= pv_prob) return PV_state;
		return -1;
		break;
	case 3:
		total_count = prob_each_node[i].p_pv + prob_each_node[i].p_dv;
		pv_prob = reset_prob / total_count * prob_each_node[i].p_pv;
		dv_prob = reset_prob / total_count * prob_each_node[i].p_dv;
		if (r <= pv_prob) return PV_state;
		r -= pv_prob;
		if (r <= dv_prob) return DV_state;
		return -1;
		break;
	case 4:
		return DV_state;
		break;
	}
	//cout << r << endl;
	/*
	if (r <= prob_each_node[i].p_us) return US_state;
	r -= prob_each_node[i].p_us;
	if (r <= prob_each_node[i].p_ui) return UI_state;
	r -= prob_each_node[i].p_ui;
	if (r <= prob_each_node[i].p_pv) return PV_state;
	r -= prob_each_node[i].p_pv;
	if (r <= prob_each_node[i].p_dv) return DV_state;
	return -1;
	*/

};

void upgrading_state_prob(prob_every_state* temp, prob_every_state* prob_each_node, int BS_number) {
	for (int j = 0; j < BS_number; ++j) { // time t to t+1
		prob_each_node[j].p_us = temp[j].p_us;
		prob_each_node[j].p_ui = temp[j].p_ui;
		prob_each_node[j].p_pv = temp[j].p_pv;
		prob_each_node[j].p_dv = temp[j].p_dv;
	}
};

void count_members_time_t(int* state_matrix, member_nodes* total_nodes_distribution, int t, int BS_number) {
	for (int j = 0; j < BS_number; ++j) {
		switch (state_matrix[j]) {
		case 1: //US
			total_nodes_distribution[t].US += 1;
			break;
		case 2: //UI
			total_nodes_distribution[t].UI += 1;
			break;
		case 3: //PV
			total_nodes_distribution[t].PV += 1;
			break;
		case 4: //DV
			total_nodes_distribution[t].DV += 1;
			break;
		}
	}
	total_nodes_distribution[t].fraction_US = (double)total_nodes_distribution[t].US / (double)BS_number;
	total_nodes_distribution[t].fraction_UI = (double)total_nodes_distribution[t].UI / (double)BS_number;
	total_nodes_distribution[t].fraction_PV = (double)total_nodes_distribution[t].PV / (double)BS_number;
	total_nodes_distribution[t].fraction_DV = (double)total_nodes_distribution[t].DV / (double)BS_number;
};
