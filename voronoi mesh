#include <cmath>
#include <vector>
#include <list>
#include <utility>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <limits>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <iomanip>

#include<chrono>
#include<cstdlib>
#include <time.h>
using namespace std;


double max(double d1, double d2, double d3, double d4)
{
	if (d1 > d2 && d1 > d3 && d1 > d4) return d1;
	if (d2 > d1 && d2 > d3 && d2 > d4) return d2;
	if (d3 > d1 && d3 > d2 && d3 > d4) return d3;
	if (d4 > d1 && d4 > d2 && d4 > d3) return d4;
}

bool check(size_t n, size_t c, size_t Nvr, size_t Nhr, size_t Nv, double sr, std::vector <size_t> cells, bool* bool1, double* x_cord, double* y_cord, double r2)
{
	size_t counter(0);
	double d1(0.0), d2(0.0), d3(0.0), d4(0.0), d(0.0);
	double yc1(0.0), yc2(0.0), yc3(0.0), yc4(0.0);
	double xc1(0.0), xc2(0.0), xc3(0.0), xc4(0.0);
	size_t* twenty1 = new size_t[20];


	twenty1[0] = n - 1 - 2 * Nv;
	twenty1[1] = n - 0 - 2 * Nv;
	twenty1[2] = n + 1 - 2 * Nv;


	twenty1[3] = n - 2 - Nv;
	twenty1[4] = n - 1 - Nv;
	twenty1[5] = n - 0 - Nv;
	twenty1[6] = n + 1 - Nv;
	twenty1[7] = n + 2 - Nv;


	twenty1[8] = n - 2;
	twenty1[9] = n - 1;
	twenty1[10] = n + 1;
	twenty1[11] = n + 2;


	twenty1[12] = n - 2 + Nv;
	twenty1[13] = n - 1 + Nv;
	twenty1[14] = n + 0 + Nv;
	twenty1[15] = n + 1 + Nv;
	twenty1[16] = n + 2 + Nv;


	twenty1[17] = n - 1 + 2 * Nv;
	twenty1[18] = n + 0 + 2 * Nv;
	twenty1[19] = n + 1 + 2 * Nv;


	xc1 = (sr * size_t(c / (Nvr)));
	xc2 = xc1;
	xc3 = xc1 + sr;
	xc4 = xc3;

	yc1 = (sr * (c % Nhr));
	yc2 = yc1 + sr;
	yc3 = yc1;
	yc4 = yc2;

	if (bool1[n] == true) return false;
	for (size_t i = 0; i < 20; i++)
	{
		//if (bool1[n] == true) continue;
		if (bool1[twenty1[i]] == false)
		{
			counter++;
			continue;
		}

		// the distances between the circle center and the corneres of the cell
		d1 = (x_cord[twenty1[i]] - xc1) * (x_cord[twenty1[i]] - xc1) + (y_cord[twenty1[i]] - yc1) * (y_cord[twenty1[i]] - yc1);
		d2 = (x_cord[twenty1[i]] - xc2) * (x_cord[twenty1[i]] - xc2) + (y_cord[twenty1[i]] - yc2) * (y_cord[twenty1[i]] - yc2);
		d3 = (x_cord[twenty1[i]] - xc3) * (x_cord[twenty1[i]] - xc3) + (y_cord[twenty1[i]] - yc3) * (y_cord[twenty1[i]] - yc3);
		d4 = (x_cord[twenty1[i]] - xc4) * (x_cord[twenty1[i]] - xc4) + (y_cord[twenty1[i]] - yc4) * (y_cord[twenty1[i]] - yc4);
		d = max(d1, d2, d3, d4);

		if (d > r2)
		{
			counter++;
		}
		else
		{
			return false;
		}
	}

	if (counter == 20)
	{

		return true;
	}
	else return false;
}



int main()
{
	//get the start time
	auto start = chrono::steady_clock::now();

	srand(time(NULL));
	double r = 0.01;
	double r2 = r * r;
	double s = r / sqrt(2);
	size_t Nv = ceil(1 / s);
	size_t Nh = ceil(1 / s);
	size_t N = Nv * Nh;

	double* x_cord = new double[N];
	double* y_cord = new double[N];
	bool* bool1 = new bool[N];

	for (size_t i = 0; i < N; i++)
	{
		bool1[i] = false;
		x_cord[i] = 0.0;  //check
		y_cord[i] = 0.0;  //check
	}

	std::vector <size_t> covered;
	std::vector <size_t> cells, cells1;
	int bug = 0;

	double x(0.0), y(0.0), X(0.0), Y(0.0);

	int* twenty = new int[20];
	size_t n = 0;
	size_t miss = 0;

	while (true)
	{
		n = size_t(((double)rand() / RAND_MAX) * N);
		if (bool1[n] == true) { continue; }
		double x = ((double)rand() / RAND_MAX) * s;
		double y = ((double)rand() / RAND_MAX) * s;

		if (n <= (2 * Nv) && n >= 0) continue;
		if (n <= (N) && n >= (N - 2 * Nv)) continue;
		if ((n % Nv) == 0 || (n % Nv) == 1) continue;
		if ((n % Nv) == (Nv - 1) || (n % Nv) == (Nv - 2)) continue;


		twenty[0] = n - 1 - 2 * Nv;
		twenty[1] = n - 0 - 2 * Nv;
		twenty[2] = n + 1 - 2 * Nv;


		twenty[3] = n - 2 - Nv;
		twenty[4] = n - 1 - Nv;
		twenty[5] = n - 0 - Nv;
		twenty[6] = n + 1 - Nv;
		twenty[7] = n + 2 - Nv;


		twenty[8] = n - 2;
		twenty[9] = n - 1;
		twenty[10] = n + 1;
		twenty[11] = n + 2;


		twenty[12] = n - 2 + Nv;
		twenty[13] = n - 1 + Nv;
		twenty[14] = n + 0 + Nv;
		twenty[15] = n + 1 + Nv;
		twenty[16] = n + 2 + Nv;


		twenty[17] = n - 1 + 2 * Nv;
		twenty[18] = n + 0 + 2 * Nv;
		twenty[19] = n + 1 + 2 * Nv;


		X = (s * size_t(n / Nv)) + x;
		Y = (s * (n % Nv)) + y;

		size_t check = 0;
		for (size_t i = 0; i < 20; i++)
		{
			if (bool1[twenty[i]] == false)
			{
				check++;
				continue;
			}
			if ((((x_cord[twenty[i]] - X) * (x_cord[twenty[i]] - X)) + ((y_cord[twenty[i]] - Y) * (y_cord[twenty[i]] - Y))) < r2)
			{
				miss++;
				break;
			}
			check++;	// the current square has a point but on distance > R
		}

		if (check == 20) // point is valid
		{
			bool1[n] = true;
			x_cord[n] = X;
			y_cord[n] = Y;
			covered.push_back(n);
			miss = 0;
		}

		if (miss > 100) break;
	}

	for (size_t i = 0; i < N; i++)
	{
		if (bool1[i] == false)
		{
			if (i <= (2 * Nv) && i >= 0) continue;
			if (i <= (N) && i >= (N - 2 * Nv)) continue;
			if ((i % Nv) == 0 || (i % Nv) == 1) continue;
			if ((i % Nv) == (Nv - 1) || (i % Nv) == (Nv - 2)) continue;
			cells.push_back(i);

		}
	}


	size_t Nr = 4 * N, Nvr = 2 * Nv, Nhr = 2 * Nh, f = 2, h = 0;
	double sr = 0.5 * s;

	/size_t twenty1 = new size_t[20];*/
	size_t g = 0;
	while (true)
	{
		g++;
		for (size_t i = 0; i < cells.size(); i++)
		{
			cells1.push_back(cells[i]);
		}

		cells.clear();

		for (size_t i = 0; i < cells1.size(); i++)    // full the cells vector
		{
			if (g == 1)
			{
				n = cells1[i];
			}
			else
			{
				n = size_t(size_t(cells1[i] / (Nvr / 2)) / (f / 2)) * Nv + size_t((cells1[i] % (Nvr / 2)) / (f / 2)); 
			}

			size_t c1, c2, c3, c4;

			// get children cells from current parent cells
			c1 = (size_t(cells1[i] / (Nvr / 2)) * Nhr * 2) + (cells1[i] % (Nhr / 2)) * 2;
			c2 = c1 + 1;
			c3 = c1 + Nhr;
			c4 = c3 + 1;

			if (check(n, c1, Nvr, Nhr, Nv, sr, cells, bool1, x_cord, y_cord, r2)) cells.push_back(c1);
			if (check(n, c2, Nvr, Nhr, Nv, sr, cells, bool1, x_cord, y_cord, r2)) cells.push_back(c2);
			if (check(n, c3, Nvr, Nhr, Nv, sr, cells, bool1, x_cord, y_cord, r2)) cells.push_back(c3);
			if (check(n, c4, Nvr, Nhr, Nv, sr, cells, bool1, x_cord, y_cord, r2)) cells.push_back(c4);

		}

		std::cout << cells.size() << std::endl;
		if (cells.size() == 0) { break; }
		for (int i = 0; i < cells.size(); i++)
		{
			n = size_t(size_t(cells[i] / (Nvr)) / (f)) * Nv + size_t((cells[i] % (Nvr)) / (f));
			if (bool1[n] == true) { h++; }
		}
		int* twenty2 = new int[20];
		size_t miss1 = 0;

		// dart throwing on cells vector
		while (true)
		{
			int m, n;
			m = size_t(((double)rand() / RAND_MAX) * (cells.size() - 1));
			n = size_t(size_t(cells[m] / (Nvr)) / (f)) * Nv + size_t((cells[m] % (Nvr)) / (f)); // no of father cell


			if (bool1[n] == true)
			{
				h++;
				cells.erase(cells.begin() + m);
				if (cells.size() == 0) break;
				continue;
			}     

			double x1 = ((double)rand() / RAND_MAX) * sr;
			double y1 = ((double)rand() / RAND_MAX) * sr;

			twenty2[0] = n - 1 - 2 * Nv;
			twenty2[1] = n - 0 - 2 * Nv;
			twenty2[2] = n + 1 - 2 * Nv;


			twenty2[3] = n - 2 - Nv;
			twenty2[4] = n - 1 - Nv;
			twenty2[5] = n - 0 - Nv;
			twenty2[6] = n + 1 - Nv;
			twenty2[7] = n + 2 - Nv;


			twenty2[8] = n - 2;
			twenty2[9] = n - 1;
			twenty2[10] = n + 1;
			twenty2[11] = n + 2;


			twenty2[12] = n - 2 + Nv;
			twenty2[13] = n - 1 + Nv;
			twenty2[14] = n + 0 + Nv;
			twenty2[15] = n + 1 + Nv;
			twenty2[16] = n + 2 + Nv;


			twenty2[17] = n - 1 + 2 * Nv;
			twenty2[18] = n + 0 + 2 * Nv;
			twenty2[19] = n + 1 + 2 * Nv;


			double X1 = (sr * size_t(cells[m] / Nvr)) + x1;
			double Y1 = (sr * (cells[m] % Nvr)) + y1;

			size_t check1 = 0;
			for (size_t i = 0; i < 20; i++)
			{
				if (bool1[twenty2[i]] == false)
				{
					check1++;
					continue;
				}
				if ((((x_cord[twenty2[i]] - X1) * (x_cord[twenty2[i]] - X1)) + ((y_cord[twenty2[i]] - Y1) * (y_cord[twenty2[i]] - Y1))) < r2) //
				{
					miss1++;
					break;
				}
				check1++;
			}

			if (check1 == 20) // point is valid
			{
				bool1[n] = true;
				x_cord[n] = X1;
				y_cord[n] = Y1;
				covered.push_back(n);
				miss1 = 0;
				cells.erase(cells.begin() + m);
			}

			if (miss1 > 100) break;
		}

		Nvr = 2 * Nvr; Nhr = 2 * Nhr; sr = 0.5 * sr; Nr = 4 * Nr; f = 2 * f; 
		cells1.clear();
	}
	std::cout << covered.size() << std::endl;

	// get the end time 
	auto end = chrono::steady_clock::now();

	//find the difference
	double elapsed_time_ns = double(std::chrono::duration_cast <std::chrono::nanoseconds> (end - start).count());

	// output
	std::cout << "elapsed time(s): " << elapsed_time_ns / 1e9 << endl;



    std::string file_name = "myplot.ps";
	std::fstream file(file_name.c_str(), std::ios::out);
	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale% one unit = one inch" << std::endl;
	double xmin(0), xmax(1);
	double ymin(0), ymax(1);
	double Lx(xmax - xmin);
	double Ly(ymax - ymin);
	double scale_x, scale_y, scale;
	double shift_x, shift_y;
	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
		shift_x = 1.0 - xmin * scale;
		shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
	}
	else
	{
		scale = scale_y;
		shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
	shift_y = 1.0 - ymin * scale;
	}
	file << shift_x << " " << shift_y << " translate" << std::endl;
	file << "/Courier findfont" << std::endl;
	file << "0.12 scalefont" << std::endl;
	file << "setfont" << std::endl;
	
	for (int i = 0; i < N; i++)
	{
		
		file << " 0 0 1 setrgbcolor" << std::endl;
		file << "newpath" << std::endl;
		file << x_cord[i] * scale << " " << y_cord[i] * scale << " " << 0.05*r * scale << " " << "0 360 " << "arc" << std::endl;;
		file << "0.0 setlinewidth" << std::endl;
		file << "fill" << std::endl;
		file << "closepath" << std::endl;
		
	}
	

	double mid_x, mid_y, normal_x, normal_y, V_x, V_y;
	double validty_check, corner_mid_x, corner_mid_y, trimming_ratio, c1_x, c1_y, c2_x, c2_y, corneri_cornerj_x, corneri_cornerj_y;
	for (size_t i = 0; i < covered.size(); i++)
	{
		size_t* fourty_four = new size_t[44];
		double* corners_x = new double[12];
		double* corners_y = new double[12];
		double* newcorners_x = new double[12];
		double* newcorners_y = new double[12];
		size_t* neighbor = new size_t[12];
		size_t* newneighbor = new size_t[12];
		
		size_t* right = new size_t[12];
		size_t* left = new size_t[12];
		bool* corners_bool = new bool[12];
		size_t no_corners = 4;
		size_t neighbour, seed;

		corners_x[0] = 0;
		corners_x[1] = 1;
		corners_x[2] = 1;
		corners_x[3] = 0;

		corners_y[0] = 0;
		corners_y[1] = 0;
		corners_y[2] = 1;
		corners_y[3] = 1;

		right[0] = 1;
		right[1] = 2;
		right[2] = 3;
		right[3] = 0;

		left[0] = 3;
		left[1] = 0;
		left[2] = 1;
		left[3] = 2;

		neighbor[0] = 0; //SIZE_MAX;
		neighbor[1] = 0;//SIZE_MAX;
		neighbor[2] = 0; //SIZE_MAX;
		neighbor[3] = 0;//SIZE_MAX;


		seed = covered[i];
		fourty_four[0] = seed - 2 - 3 * Nv;
		fourty_four[1] = seed - 1 - 3 * Nv;
		fourty_four[2] = seed - 0 - 3 * Nv;
		fourty_four[3] = seed + 1 - 3 * Nv;
		fourty_four[4] = seed + 2 - 3 * Nv;


		fourty_four[5]  = seed - 3 - 2 * Nv;
		fourty_four[6]  = seed - 2 - 2 * Nv;
		fourty_four[7]  = seed - 1 - 2 * Nv;
		fourty_four[8]  = seed - 0 - 2 * Nv;
		fourty_four[9]  = seed + 1 - 2 * Nv;
		fourty_four[10] = seed + 2 - 2 * Nv;
		fourty_four[11] = seed + 3 - 2 * Nv;


		fourty_four[12] = seed - 3 - Nv;
		fourty_four[13] = seed - 2 - Nv;
		fourty_four[14] = seed - 1 - Nv;
		fourty_four[15] = seed - 0 - Nv;
		fourty_four[16] = seed + 1 - Nv;
		fourty_four[17] = seed + 2 - Nv;
		fourty_four[18] = seed + 3 - Nv;


		fourty_four[19] = seed - 3;
		fourty_four[20] = seed - 2;
		fourty_four[21] = seed - 1;
		fourty_four[22] = seed + 1;
		fourty_four[23] = seed + 2;
		fourty_four[24] = seed + 3;


		fourty_four[25] = seed - 3 + Nv;
		fourty_four[26] = seed - 2 + Nv;
		fourty_four[27] = seed - 1 + Nv;
		fourty_four[28] = seed + 0 + Nv;
		fourty_four[29] = seed + 1 + Nv;
		fourty_four[30] = seed + 2 + Nv;
		fourty_four[31] = seed + 3 + Nv;


		fourty_four[32] = seed - 3 + 2 * Nv;
		fourty_four[33] = seed - 2 + 2 * Nv;
		fourty_four[34] = seed - 1 + 2 * Nv;
		fourty_four[35] = seed + 0 + 2 * Nv;
		fourty_four[36] = seed + 1 + 2 * Nv;
		fourty_four[37] = seed + 2 + 2 * Nv;
		fourty_four[38] = seed + 3 + 2 * Nv;


		fourty_four[39] = seed - 2 + 3 * Nv;
		fourty_four[40] = seed - 1 + 3 * Nv;
		fourty_four[41] = seed - 0 + 3 * Nv;
		fourty_four[42] = seed + 1 + 3 * Nv;
		fourty_four[43] = seed + 2 + 3 * Nv;


		for (size_t j = 0; j < 44; j++)
		{
			if (bool1[fourty_four[j]] == false) continue;

			for (size_t i = 0; i < no_corners; i++)
			{
				corners_bool[i] = true;
			}

			neighbour = fourty_four[j];

			mid_x = (x_cord[seed] + x_cord[neighbour]) / 2;
			mid_y = (y_cord[seed] + y_cord[neighbour]) / 2;

			normal_x = x_cord[neighbour] - x_cord[seed];
			normal_y = y_cord[neighbour] - y_cord[seed];

			for (size_t k = 0; k < no_corners; k++)
			{
				V_x = corners_x[k] - mid_x;
				V_y = corners_y[k] - mid_y;
				validty_check = V_x * normal_x + V_y * normal_y;
				if (validty_check > 0)
				{
					corners_bool[k] = false;
				}
			}

			int h = 0;
			for (size_t k = 0; k < no_corners; k++)
			{
				if ((corners_bool[k] == true) && (corners_bool[right[k]] == false))
				{
					corner_mid_x = mid_x - corners_x[k];
					corner_mid_y = mid_y - corners_y[k];

					corneri_cornerj_x = corners_x[right[k]] - corners_x[k];
					corneri_cornerj_y = corners_y[right[k]] - corners_y[k];

					trimming_ratio = (corner_mid_x * normal_x + corner_mid_y * normal_y) / (corneri_cornerj_x * normal_x + corneri_cornerj_y * normal_y);
					c1_x = corners_x[k] + trimming_ratio * (corneri_cornerj_x);
					c1_y = corners_y[k] + trimming_ratio * (corneri_cornerj_y);
				}

				if ((corners_bool[k] == true) && (corners_bool[left[k]] == false))
				{
					h = k;
					corner_mid_x = mid_x - corners_x[k];
					corner_mid_y = mid_y - corners_y[k];

					corneri_cornerj_x = corners_x[left[k]] - corners_x[k];
					corneri_cornerj_y = corners_y[left[k]] - corners_y[k];

					trimming_ratio = (corner_mid_x * normal_x + corner_mid_y * normal_y) / (corneri_cornerj_x * normal_x + corneri_cornerj_y * normal_y);
					c2_x = corners_x[k] + trimming_ratio * (corneri_cornerj_x);
					c2_y = corners_y[k] + trimming_ratio * (corneri_cornerj_y);
				}
			}

			
			int counterr(-1);
			for (size_t k = 0; k < no_corners; k++)
			{
				if ((corners_bool[k] == true) && (corners_bool[right[k]] == false))
				{
					counterr++;
					newcorners_x[counterr] = corners_x[k];
					newcorners_y[counterr] = corners_y[k];
					newneighbor[counterr] = neighbor[k];
					counterr++;
					newcorners_x[counterr] = c1_x;
					newcorners_y[counterr] = c1_y;
					newneighbor[counterr] = neighbour;
					counterr++;
					newcorners_x[counterr] = c2_x;
					newcorners_y[counterr] = c2_y;
					newneighbor[counterr] = neighbor[left[h]];
					continue;
				}
				if (corners_bool[k] == true)
				{
					counterr++;
					newcorners_x[counterr] = corners_x[k];
					newcorners_y[counterr] = corners_y[k];
					newneighbor[counterr] = neighbor[k];
				}
			}

			
			no_corners = counterr+1;
			for (size_t k = 0; k < no_corners; k++)
			{
				corners_x[k] = newcorners_x[k];
				corners_y[k] = newcorners_y[k];
			}

			for (size_t k = 0; k < no_corners; k++)
			{
				neighbor[k] = newneighbor[k];
				neighbor[k] = newneighbor[k];
			}

			for (size_t k = 0; k < no_corners; k++)
			{
				if (k == 0)
				{
					right[k] = 1;
					left[k] = no_corners - 1;
					continue;
				}
				if (k == no_corners - 1)
				{
					right[k] = 0;
					left[k] = no_corners - 2;
					continue;
				}
				right[k] = k + 1;
				left[k] = k - 1;
			}
		}

		
		
		for (size_t k = 0; k < no_corners; k++)
		{
			if (neighbor[k] == 0 || neighbor[left[k]] == 0) continue;

			file << "newpath" << std::endl;
			file << x_cord[neighbor[k]] * scale << " " << (y_cord[neighbor[k]]) * scale << " ";
			file << "moveto" << std::endl;
			file << x_cord[neighbor[left[k]]] * scale << " " << (y_cord[neighbor[left[k]]]) * scale << " ";
			file << "lineto" << std::endl;
			file << ".005 setlinewidth" << std::endl; // Inactive seeds
			file << "0 0 0 setrgbcolor" << std::endl;
			file << "stroke" << std::endl;

			file << "newpath" << std::endl;
			file << x_cord[neighbor[k]] * scale << " " << (y_cord[neighbor[k]]) * scale << " ";
			file << "moveto" << std::endl;
			file << x_cord[seed] * scale << " " << (y_cord[seed]) * scale << " ";
			file << "lineto" << std::endl;
			file << ".005 setlinewidth" << std::endl; // Inactive seeds
			file << "0 0 0 setrgbcolor" << std::endl;
			file << "stroke" << std::endl;

			file << "newpath" << std::endl;
			file << x_cord[neighbor[left[k]]] * scale << " " << (y_cord[neighbor[left[k]]]) * scale << " ";
			file << "moveto" << std::endl;
			file << x_cord[seed] * scale << " " << (y_cord[seed]) * scale << " ";
			file << "lineto" << std::endl;
			file << ".005 setlinewidth" << std::endl; // Inactive seeds
			file << "0 0 0 setrgbcolor" << std::endl;
			file << "stroke" << std::endl;
		}
	}
	file << "showpage" << std::endl;
}
