#include "pch.h"
#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>

using namespace std;

//This program runs affine curve-shortening flow (ACSF) on an initial curve.
//The curve is given by a circular list of points. The number of points
//stays constant throughout the execution.
//At each step, for each point, the radius of curvature and normal are estimated
//by taking the unique circle through the point and its two neighbors.
//Then all the points are moved in unison for time ct*dmin^(4/3), where dmin
//is the minimum distance between two consecutive points, and ct is a fixed constant.

//Maximum number of points allowed:
const int MAXN = 5000;

//Number of points
int N;

const double pi = acos(-1);

//The circular list of points. The first two points duplicated at the end of the list
//for convenience.
double gammax[MAXN + 2], gammay[MAXN + 2];

//Total time elapsed:
double t;

//Constant for the tangent velocity component (0 = no tangent velocity).
//The tangent velocity is used to try and keep the points from bunching together
//at the sharp ends of the curve. The tangent velocity should have no effect on the
//ACSF evolution.
const int tangent_c = 10;

//Constant used in vec2ctr to determine whether the three points are almost aligned:
const double almost_zero = 1e-15;

//Given three points p1, p2, p3, returns the vector from p2 to the
//circle center.
//If the three points are aligned, or almost aligned, we return false.
bool vec2ctr(double, double, double, double, double, double, double *, double *);

//Curve initialization functions. They set t to 0:
void initroundtriangle();
void initcircle();
void initpear();
void initsquare();
void inithalfcircle();
void initpolygon(int, const double *, const double *);
void initpolygon_uniform(int, const double *, const double *);
void inittriangle();
void initstar();
double spiralf(double);
double spiralfinv(double);
void initspiral();
void initinf();
void initneedle();
void initcamelfish();
void initbowtie();
void initswan();

//Returns the squared distance between point i and point i+1 in gamma
double dsq(int i);

//Finds the minimum distance between two consecutive points in gamma
double find_dmin();

//Takes a step in the flow. The first input parameter is the inverse radius exponent
//(1 for CSF, 1/3 for ACSF). The second input parameter is
//the time step to take.
//We give a tangential component to the velocities of the points.
//This mitigates their tendency to bunch toegher.
//The tangential component is not supposed to affect the flow evolution.
void power_flow_step(double, double);

//These two functions are wrappers for power_flow_step:
void CSF_step(double);
void ACSF_step(double);

//Cheap scientific notation workaround:
void my_out_double(double, ostream &);
//(A better solution would be to output in JSON format; see
//https://mathematica.stackexchange.com/questions/146313/)

//Output the current curve (with first point = last point for convenience)
//in a format suitable for Mathematica:
void output_curve_Mathematica(ostream &);

//Output the current curve (with first point = last point for convenience)
//in a format suitable for C:
void output_curve_C(ostream &);

//Returns the perimeter (Euclidean length) of the curve:
double perimeter();

//Like initandrun(), but decides when to output based on the perimeter
void initandrun_perimeter(void(*initf)(), bool, int, double, double, int, int, string);

bool vec2ctr(double x1, double y1, double x2, double y2, double x3, double y3,
	double *outx, double *outy)
{
	double det = x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1);

	if (abs(det) < almost_zero)
		return false;

	double xctr, yctr;
	xctr = (x1*x1*(y3 - y2) + x2 * x2*(y1 - y3) + x3 * x3*(y2 - y1)
		+ (y1 - y2)*(y2 - y3)*(y3 - y1)) / (2 * det);
	yctr = -(y1*y1*(x3 - x2) + y2 * y2*(x1 - x3) + y3 * y3*(x2 - x1)
		+ (x1 - x2)*(x2 - x3)*(x3 - x1)) / (2 * det);
	*outx = xctr - x2;
	*outy = yctr - y2;
	return true;
}

//This function initializes the curve to the curve given by
//( ((1-sin(t))/2)^2, ((1-sin(t+2))/2)^1.3 )
//for 0 <= t <= 2*pi.
void initroundtriangle()
{
	if (N > MAXN)
	{
		cout << "N too large." << endl;
		return;
	}

	int i;
	for (i = 0; i <= N + 1; i++)
	{
		gammax[i] = pow((1 - sin(i * 2 * pi / N)) / 2, 2);
		gammay[i] = pow((1 - sin(i * 2 * pi / N + 2)) / 2, 1.3);
	}
	t = 0;
}

//Initializes the curve to a circle of diameter 1:
void initcircle()
{
	if (N > MAXN)
	{
		cout << "N too large." << endl;
		return;
	}
	int i;
	for (i = 0; i <= N + 1; i++)
	{
		gammax[i] = cos(i * 2 * pi / N) / 2;
		gammay[i] = sin(i * 2 * pi / N) / 2;
	}
	t = 0;
}

//Initializes the curve to the curve given by
//((1 - sin(t + 2))/2)^1.5, 0.7 ((1 - sin(t))/2 - 0.3 ((1 - Sin(t + 0.6)) / 2) ^ 15)
//for 0 <= t <= 2*pi.
void initpear()
{
	if (N > MAXN)
	{
		cout << "N too large." << endl;
		return;
	}

	int i;
	for (i = 0; i <= N + 1; i++)
	{
		gammax[i] = pow((1 - sin(i * 2 * pi / N + 2)) / 2, 1.5);
		gammay[i] = 0.7 * ((1 - sin(i * 2 * pi / N)) / 2
			- 0.3 * pow((1 - sin(i * 2 * pi / N + 0.6)) / 2, 15));
	}
	t = 0;
}

//Initializes the curve to the unit square.
//This needs N to be a multiple of 4.
void initsquare()
{
	if (N > MAXN)
	{
		cout << "N too large." << endl;
		return;
	}

	int i;
	if (N % 4 != 0)
	{
		cout << "Error: N should be a multiple of 4 for initsquare." << endl;
		return;
	}
	int m = N / 4;
	for (i = 0; i < m; i++)
	{
		gammax[i] = 0;
		gammay[i] = (double)i / m;
		gammax[m + i] = (double)i / m;
		gammay[m + i] = 1;
		gammax[2 * m + i] = 1;
		gammay[2 * m + i] = 1 - (double)i / m;
		gammax[3 * m + i] = 1 - (double)i / m;
		gammay[3 * m + i] = 0;
	}
	gammax[N] = gammax[0];
	gammay[N] = gammay[0];
	gammax[N + 1] = gammax[1];
	gammay[N + 1] = gammay[1];
	t = 0;
}

//Initializes the curve to a half-circle:
void inithalfcircle()
{
	if (N > MAXN)
	{
		cout << "N too large." << endl;
		return;
	}

	int m = N / (1 + pi / 2);
	int i;

	for (i = 0; i < m; i++)
	{
		gammax[i] = (double)i / m;
		gammay[i] = 0;
	}
	for (i = m; i < N; i++)
	{
		gammax[i] = 0.5*(1 + cos(pi*(i - m) / (N - m)));
		gammay[i] = 0.5*sin(pi*(i - m) / (N - m));
	}
	gammax[N] = gammax[0];
	gammay[N] = gammay[0];
	gammax[N + 1] = gammax[1];
	gammay[N + 1] = gammay[1];
	t = 0;
}

//Initializes the curve to a given polygon:
void initpolygon(int num_sides, const double *xx, const double *yy)
{
	if (N > MAXN)
	{
		cout << "N too large." << endl;
		return;
	}

	int imin, imax = 0, side_n, i;

	for (side_n = 0; side_n < num_sides; side_n++)
	{
		imin = imax; //Old max becomes new min
		imax = N * (1 + side_n) / num_sides;
		for (i = imin; i < imax; i++)
		{
			double alpha = (double)(i - imin) / (imax - imin);
			gammax[i] = alpha * xx[(side_n + 1) % num_sides] + (1 - alpha)*xx[side_n];
			gammay[i] = alpha * yy[(side_n + 1) % num_sides] + (1 - alpha)*yy[side_n];
		}
	}
	gammax[N] = gammax[0];
	gammay[N] = gammay[0];
	gammax[N + 1] = gammax[1];
	gammay[N + 1] = gammay[1];
	t = 0;
}

//Initializes the curve to a given polygon, spacing points uniformly:
void initpolygon_uniform(int num_sides, const double *xx, const double *yy)
{
	int i;

	//Find polygon perimeter:
	double per = 0, deltax, deltay;
	for (i = 0; i < num_sides; i++)
	{
		deltax = xx[(i + 1) % num_sides] - xx[i];
		deltay = yy[(i + 1) % num_sides] - yy[i];
		per += sqrt(deltax*deltax + deltay * deltay);
	}

	double persofar = 0; //perimeter seen so far
	int edgemin, edgemax = 0;
	for (i = 0; i < num_sides; i++)
	{
		//Set new edgemin to old edgemax:
		edgemin = edgemax;

		//Set new edgemax:
		if (i == num_sides - 1)
			edgemax = N;
		else
		{
			//Update persofar:
			deltax = xx[(i + 1) % num_sides] - xx[i];
			deltay = yy[(i + 1) % num_sides] - yy[i];
			persofar += sqrt(deltax*deltax + deltay * deltay);
			edgemax = N * persofar / per;
		}

		//Populate gammax[j], gammay[j] for edgemin <= j < edgemax
		//with points uniformly from i-th polygon edge e_i, so that
		//endpoints of e_i are at j==edgemin and j==edgemax.
		for (int j = edgemin; j < edgemax; j++)
		{
			double alpha = (double)(j - edgemin) / (edgemax - edgemin);
			gammax[j] = alpha * xx[(i + 1) % num_sides] + (1 - alpha) * xx[i];
			gammay[j] = alpha * yy[(i + 1) % num_sides] + (1 - alpha) * yy[i];
		}
	}

	gammax[N] = gammax[0];
	gammay[N] = gammay[0];
	gammax[N + 1] = gammax[1];
	gammay[N + 1] = gammay[1];
	t = 0;
}


//Initializes the curve to a triangle with vertices (0,0), (1,3/4), (2/5,1):
void inittriangle()
{
	const double xx[3] = { 0, 1, 2. / 5 }, yy[3] = { 0, 3. / 4, 1 };
	initpolygon(3, xx, yy);
}

//Initializes the curve to a star:
void initstar()
{
	double xx[10], yy[10];
	int i;
	for (i = 0; i < 5; i++)
	{
		xx[2 * i] = cos(2 * pi*i / 5);
		yy[2 * i] = sin(2 * pi*i / 5);
		xx[2 * i + 1] = 4 * cos((2 * i + 1)*pi / 5);
		yy[2 * i + 1] = 4 * sin((2 * i + 1)*pi / 5);
	}

	initpolygon(10, xx, yy);
}

//Two functions used by initspiral:
double spiralf(double x)
{
	return sqrt(x + 1) - 1;
}

double spiralfinv(double x)
{
	return (1 + x)*(1 + x) - 1;
}

//Initializes the curve to a spiral:
void initspiral()
{
	int N1 = 0.558*N, N2 = 0.43*N, N3 = N - N1 - N2 - 1, i;
	double a;

	for (i = 0, a = spiralfinv(7 * pi); i < N1; i++, a -= spiralfinv(7 * pi) / N1)
	{
		gammax[i] = -spiralf(a)*cos(spiralf(a));
		gammay[i] = -spiralf(a)*sin(spiralf(a));
	}
	for (i = N1, a = 0; i <= N1 + N2; i++, a += spiralfinv(6 * pi) / N2)
	{
		gammax[i] = spiralf(a)*cos(spiralf(a));
		gammay[i] = spiralf(a)*sin(spiralf(a));
	}
	for (i = N1 + N2 + 1, a = pi * (1 - 1 / N3); i < N; i++, a -= pi / N3)
	{
		gammax[i] = pi * (6.5 + 0.5*cos(a));
		gammay[i] = pi * 0.5*sin(a);
	}
	gammax[N] = gammax[0];
	gammay[N] = gammay[0];
	gammax[N + 1] = gammax[1];
	gammay[N + 1] = gammay[1];
	t = 0;
}

//Initializes the curve to a self-intersecting "infinity" with two unequal lobes:
void initinf()
{
	int N1 = 0.75*N, N2 = N - N1, i;

	for (i = 0; i < N1; i++)
	{
		gammax[i] = cos(2 * pi*i / N1);
		gammay[i] = sin(2 * pi*i / N1)*sin(pi*i / N1);
	}
	for (i = N1; i < N; i++)
	{
		gammax[i] = (4 - cos(2 * pi*(i - N1) / N2)) / 3;
		gammay[i] = sin(2 * pi*(i - N1) / N2)*sin(pi*(i - N1) / N2) / 3;
	}
	gammax[N] = gammax[0];
	gammay[N] = gammay[0];
	gammax[N + 1] = gammax[1];
	gammay[N + 1] = gammay[1];
	t = 0;
}

//Initializes the curve to a square with a needle:
void initneedle()
{
	double xx[7] = { 0,1,1,0.5,0.5,0.5,0 };
	double yy[7] = { 0,0,1,1,0.2,1,1 };

	initpolygon(7, xx, yy);
}

//Initializes the curve to a polygon with a self-intersection and an inflection point
//disappear soon:
void initcamelfish()
{
	double xx[8] = { 0.,0.16,0.4,0.64,0.94,1.,0.56,0.52 };
	double yy[8] = { 0., 0.81, 0.45, 1., 0.3, 0.45, 0.07, 0.13 };
	initpolygon_uniform(8, xx, yy);
}

//Initializes the curve to a bowtie:
void initbowtie()
{
	double xx[4] = { 1,1,-1,-1 };
	double yy[4] = { 1,-1,1,-1 };
	initpolygon_uniform(4, xx, yy);
}

//Initializes the curve to the curve given by
//(sin(x) + 3 sin(x + .5)^2, cos(x) + 2 ((sin(x + .5) + 1)/2)^7)
//for 0 <= t <= 2*pi.
void initswan()
{
	if (N > MAXN)
	{
		cout << "N too large." << endl;
		return;
	}

	int i;
	for (i = 0; i <= N + 1; i++)
	{
		double x = 2 * pi*i / N;
		gammax[i] = sin(x) + 3 * pow(sin(x + .5), 2);
		gammay[i] = cos(x) + 2 * pow((sin(x + .5) + 1) / 2, 7);
	}
	t = 0;
}

double dsq(int i)
{
	return (gammax[i + 1] - gammax[i])*(gammax[i + 1] - gammax[i])
		+ (gammay[i + 1] - gammay[i])*(gammay[i + 1] - gammay[i]);
}

double find_dmin()
{
	double min_dsq = dsq(0);
	double curr_dsq;
	int i;

	for (i = 1; i <= N - 1; i++)
	{
		curr_dsq = dsq(i);
		if (curr_dsq < min_dsq)
			min_dsq = curr_dsq;
	}
	return sqrt(min_dsq);
}

void power_flow_step(double exp, double dt)
{
	if (N > MAXN)
	{
		cout << "N too large." << endl;
		return;
	}

	double newgammax[MAXN + 2], newgammay[MAXN + 2];
	double vx, vy, rad_sq;
	int i;

	//Recall that point 0 = point N, and point 1 = point N+1.
	for (i = 1; i <= N; i++)
	{
		bool ok = vec2ctr(gammax[i - 1], gammay[i - 1], gammax[i], gammay[i],
			gammax[i + 1], gammay[i + 1], &vx, &vy);
		rad_sq = vx * vx + vy * vy;
		//Give the normal vector the right length:
		vx *= dt / pow(rad_sq, (exp + 1) / 2);
		vy *= dt / pow(rad_sq, (exp + 1) / 2);
		//We treat separately the case where the three points are aligned:
		if (!ok)
		{
			vx = 0;
			vy = 0;
		}
		//Tangent vector:
		double wx, wy;
		wx = vy;
		wy = -vx;
		//Invert the tangent vector if necessary, so it points towards point i+1:
		if (gammax[i] * gammay[i + 1] - gammax[i + 1] * gammay[i]
			- gammax[i - 1] * gammay[i + 1] + gammax[i + 1] * gammay[i - 1]
			+ gammax[i - 1] * gammay[i] - gammax[i] * gammay[i - 1] < 0)
		{
			wx = -wx;
			wy = -wy;
		}
		//Give the tangent vector a length that depends on the ratio between the distances
		//between point i and points i-1 and i+1. (If this ratio is 1 then the length
		//will be 0.)
		double tan_fac = tangent_c * log(dsq(i) / dsq(i - 1));
		wx *= tan_fac;
		wy *= tan_fac;
		newgammax[i] = gammax[i] + vx + wx;
		newgammay[i] = gammay[i] + vy + wy;
	}
	//Copy new array into old:
	for (i = 1; i <= N; i++)
	{
		gammax[i] = newgammax[i];
		gammay[i] = newgammay[i];
	}
	//Copy points 0 and N+1:
	gammax[0] = gammax[N];
	gammay[0] = gammay[N];
	gammax[N + 1] = gammax[1];
	gammay[N + 1] = gammay[1];
	//Update global time:
	t += dt;
}

void CSF_step(double dt)
{
	power_flow_step(1, dt);
}

void ACSF_step(double dt)
{
	power_flow_step(1. / 3, dt);
}

void my_out_double(double d, ostream &s)
{
	int exp;
	double fr;
	fr = frexp(d, &exp);
	s << fr << "*2^" << exp;
}

void output_curve_Mathematica(ostream &s)
{
	bool first = true;
	int i;

	s << "{";
	for (i = 0; i <= N; i++)
	{
		if (!first)
			s << ", ";
		first = false;
		s << "{";
		my_out_double(gammax[i], s);
		s << ", ";
		my_out_double(gammay[i], s);
		s << "}";
	}
	s << "}";
}

void output_curve_C(ostream &s)
{
	int i;

	s << N + 1 << endl;
	for (i = 0; i <= N; i++)
	{
		s << gammax[i] << " " << gammay[i] << " ";
	}
	s << endl;
}

double perimeter()
{
	double p = 0;
	int i;
	for (i = 0; i < N; i++)
	{
		double deltax = gammax[i + 1] - gammax[i], deltay = gammay[i + 1] - gammay[i];
		p += sqrt(deltax * deltax + deltay * deltay);
	}
	return p;
}

//If, for example, we want to print when the perimeter reaches 95%, 90%, and 85%, then we should
//set nparts to 20 (because 100% / 20 = 5%) and printparts to 3.
void initandrun_perimeter(void(*initf)(), bool affine, int desired_N, double ct,
	double min_dt, int nparts, int printparts, string filename)
{
	clock_t runt0 = clock(), runt1;

	ofstream Cfile, mathfile;
	Cfile.open(filename + " C.txt");
	mathfile.open(filename + " math.txt");

	N = desired_N;

	initf();

	output_curve_C(Cfile);

	mathfile << "flowlist = {" << endl;
	output_curve_Mathematica(mathfile);

	int i;
	int p = 1;

	double initialper = perimeter(), currper;

	cout << "Initial perimeter: " << initialper << endl;

	i = 0;
	while (p <= printparts)
	{
		double dt;

		if (affine)
			dt = ct * pow(find_dmin(), 4. / 3);
		else
			dt = ct * pow(find_dmin(), 2);

		if (dt < min_dt)
			dt = min_dt;

		if (affine)
			ACSF_step(dt);
		else
			CSF_step(dt);

		i++;

		currper = perimeter();

		if (currper <= initialper * (nparts - p) / nparts)
		{
			output_curve_C(Cfile);

			mathfile << "," << endl;
			output_curve_Mathematica(mathfile);

			cout << 100 * (nparts - p) / nparts << "% perimeter. t = " << t << endl;
			p++;
		}

		if (i % 100000 == 0)
			cout << "i = " << i << ". t = " << t << ". dt = " << dt
			<< ". perimeter = " << perimeter() << endl;
	}

	Cfile << "0" << endl;
	Cfile.close();

	mathfile << "};" << endl;
	mathfile.close();

	cout << "Final perimeter: " << currper << " i = " << i << endl;

	runt1 = clock();
	cout << "Running time: " << ((double)runt1 - runt0) / CLOCKS_PER_SEC << " s." << endl;
}

int main()
{
	initandrun_perimeter(initcamelfish, true, 1000, 0.0003, 3e-9, 50, 15, "ACSF camelfish hires");
	return 0;
}
