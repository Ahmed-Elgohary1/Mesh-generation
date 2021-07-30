#include<bits/stdc++.h>


using namespace std;
struct seed
{
    double x;
    double y;
    int check_bool1;
};
typedef long long ll;
	size_t counter=0;double radius=0.01;double s = radius / sqrt(2);
	int nc = ceil(1 / s);
	int nr = nc;
	double ymax = nc * s;
	double xmax = ymax;
	ll N = nc*nc;
	ll element=N;
	int g=0;
	size_t false_try(0);
		vector <seed>storing(N);    //it's a container that contain [x,y coord. and bool check]  IT's main function is to facilitate the check of disk free

std::vector<double> px; std::vector<double> py;

int check (ll randN,double x, double y)
{
    vector<int> neighboring_cells= { randN + 1 ,randN + 2,randN - 1,randN - 2 ,randN + nc ,randN + 2*nc,randN - nc,randN - 2*nc ,randN + 1 + nc,randN + 2 + nc,randN + 1 + 2 * nc, randN - 1 + nc,randN - 2 + nc,randN - 1 + 2 * nc ,randN + 1 - nc,randN + 2 - nc,randN + 1 - 2 * nc,randN - 1 - nc,randN - 2 - nc,randN - 1 - 2 * nc };

		if (storing[randN].check_bool1 == 0)
		{
			for (size_t i = 0; i < 20; i++)
			{
				if (storing[neighboring_cells[i]].check_bool1 == 1)
				{
					if (sqrt(((x - storing[neighboring_cells[i]].x) * (x - storing[neighboring_cells[i]].x)) + ((y - storing[neighboring_cells[i]].y) * (y - storing[neighboring_cells[i]].y))) <= radius)
					{
						counter++;
						break;
					}
				}
			}
		}
		else { counter++; }

		if (counter == 0)
		{

			counter = 0;
			false_try = 0;
			return 1;
		}
		else
		{
			false_try++;
			counter = 0;
			return 0;
		}
}


float ConstantShift(ll par,int x)
{
    int coun=0;double step=s/pow(2,g),r=par ,it=g,xadd=0,yadd=0;
    while(it)
    {
    par=ceil((float)r/pow(4,coun));
int dummy=(int)(par-1)%4;
    if(dummy==2||dummy==3)
    {
    xadd+=step;
    }
    if(dummy==2||dummy==1)
    {
    yadd+=step;
    }
    step*=2;
    it--;coun++;
    }
if(x)
{
    return xadd;
}
else
    return yadd;
}




int IsCovered (ll p)    // To check the child cell is covered totally by another disk so reject it from selection pool or not
{
    double step=s/pow(2,g);int z=ceil(p/pow(4,g))-1;
   float x_shift=ConstantShift(p,1), y_shift=ConstantShift(p,0);
    int constant1 = ((z) / nc);
		double x = constant1*s;x+=x_shift;
		double y = (z % nc) * s;y+=y_shift;

    int ch1=check (z, x,y);
    int ch2= check (z, x,y+step);
    int ch3= check (z, x+step,y+step);
    int ch4= check (z, x+step,y);
    if(ch1||ch2||ch3||ch4)
    {
        return 0;
    }

    else
        return 1;
}



int eraser (ll itre,int g,vector<int>&operations) // to delete the whole cell with its children
{
int start,ending=pow(4,g),con=0;
    start=itre-ending;
      if(start<0)
          start =0;
int w= ceil((float)operations[itre]/ending);
    for(start;start<=itre+ending;start++)
    {
        int z= ceil((float)operations[start]/ending);

        if(w==z)
        {
        while(ending--)
        {
        operations.erase(operations.begin()+start);
        con++;
        }
        break;
        }
    }
return con;
}




void generation (vector<int>&operation)
{

   ll n= operation.size(), counter=0;

    while(n>0)
    {
     for(int j=1;j<=4;j++)
       {
         int z=(operation[0]-1)*4+j;
            if(IsCovered(z))
                continue;
         operation.push_back(z);
       }
              counter++;

       eraser(0,0,operation) ;
       n= operation.size();n-=(pow(4,1))*counter;
    }

}

int plot_circle(std::string file_name, std::vector<double> px, std::vector<double> py, double rad);




void dart(vector<int>&oper,double se)
{

    while (element)
	{

		int dummy= rand() % element;// random cell is selected

        int ParentCell =ceil((oper[dummy])/pow(4,g))-1;

		float rand_delta_x = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / se));
		float rand_delta_y = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / se));
		int constant1 = ((ParentCell) / nc);
		double x = constant1*s + rand_delta_x;
		x+=ConstantShift(oper[dummy],1);
		double y = (ParentCell % nc) * s + rand_delta_y; // random x,y is produced with the limits of the picked cell
		y+=ConstantShift(oper[dummy],0);


if(check (ParentCell,x, y))   // disk free seed condition is checked
{

            storing[ParentCell].check_bool1 = 1;
			storing[ParentCell].x = x;
			storing[ParentCell].y = y;
			px.push_back(x);
			py.push_back(y);
			false_try = 0;
             //delete the cell or all children of the cell from the selection pool
            element-= eraser (dummy,g,oper);
}
		if (false_try > .2 *N)
            {

                generation(oper);
                g++;
                se/=2;
                element=oper.size();cout<<element<<endl;
            }

	}


}
int main()
{
	auto start = chrono::steady_clock::now();

	cout << s << std::endl;
	vector <int>oper(N);        //it's a container that contain the position  of cell "selection pool"

for(int i=0;i<N;i++)
{
    oper[i]=i+1;
}
dart(oper,s);           // dart function is responsible for produce bias free and disk free seeds


	auto end = chrono::steady_clock::now();
	cout<<N<<"  "<<oper.size();
plot_circle("myplot.ps", px, py, radius);   // visualize  the disks
	double elapsed_time_ns = double(std::chrono::duration_cast <std::chrono::nanoseconds> (end - start).count());
	std::cout << "elapsed time(s): " << elapsed_time_ns / 1e9 << endl;

	return 0;
}




int plot_circle(std::string file_name, std::vector<double> px, std::vector<double> py, double rad)
{

	std::fstream file(file_name.c_str(), std::ios::out);
	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale% one unit = one inch" << std::endl;

	double xmin(-1), xmax(1);
	double ymin(-1), ymax(1);

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

	for (int i(0); i < px.size(); i++)
	{
		file << "newpath" << std::endl;
		file << px[i] * scale << " " << py[i] * scale << " " << rad * scale << " 0 360 arc" << std::endl;
		file << "closepath" << std::endl;
		file << "1 1 0 setrgbcolor" << std::endl; // yellow filled circle(big area)
		file << "fill" << std::endl;
		file << "0.0 setlinewidth" << std::endl;
	}
	for (int i(0); i < px.size(); i++)
	{
		file << "newpath" << std::endl;
		file << px[i] * scale << " " << py[i] * scale << " " << rad * scale << " 0 360 arc" << std::endl; // big stroke black circle(premiter)
		file << "closepath" << std::endl;
		file << "0 0 0 setrgbcolor" << std::endl;
		file << "stroke" << std::endl;
		file << "0.0 setlinewidth" << std::endl;
	}
	for (int i(0); i < px.size(); i++)
	{
		file << "newpath" << std::endl;
		file << px[i] * scale << " " << py[i] * scale << " " << 0.001 * scale << " 0 360 arc" << std::endl;
		file << "closepath" << std::endl;
		file << "0 0 0 setrgbcolor" << std::endl; // Inactive seeds
		file << "fill" << std::endl; // black filled small dot(small area)
		file << "0.0 setlinewidth" << std::endl;
	}

	file << "newpath" << std::endl; // code to drawing the square
	file << 0 * scale << " " << (0) * scale << " ";
	file << "moveto" << std::endl;
	file << 0 * scale << " " << 1.00409 * scale << " ";
	file << "lineto" << std::endl;
	file << 1.00409 * scale << " " << 1.00409 * scale << " ";
	file << "lineto" << std::endl;
	file << 1.00409 * scale << " " << 0 * scale << " ";
	file << "lineto" << std::endl;
	file << 0 * scale << " " << 0 * scale << " ";
	file << "lineto" << std::endl;
	file << ".005 setlinewidth" << std::endl; // Inactive seeds
	file << "0 0 0 setrgbcolor" << std::endl;
	file << "stroke" << std::endl;
	file << "showpage" << std::endl;
	return 0;
}
