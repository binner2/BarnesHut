// generate_data.cpp : generate random datas
//

#include <iostream>
#include <fstream>

using namespace std;

const int NDIM = 3 ;
int main(int argc, char* argv[])
{
	cout<<"this program generates random datas for a file like this :";
	cout<<endl<<endl;
	cout<<"n	//number of particles"<<endl;
	cout<<"t	//start of simulation time "<<endl;
	cout<<"t_end	//end of simulation time"<<endl;
	cout<<"dt	//time step of simulation"<<endl;
	//cout<<"dt_out	//output time step of simultaion "<<endl;	
	cout<<"mass[1]	x[1]	y[1]	z[1]	vx[1]	vy[1]	vz[1]"<<endl;
	cout<<"mass , x,y,z coordinates and x,y,z velocities of particles.."<<endl;
	cout<<"....."<<endl;
	cout<<endl;

	char filename[25]=" ";
	cout<<"enter filename = " ; cin.get ( filename , 25 );
	
	ofstream outfile( filename );

	int n;
	cout<<"enter number of particle (n) = ";	cin >> n;
	outfile << n <<endl;

	float t;
	cout<<"simulation start time (t) = ";	cin >> t;
	outfile << t <<endl;

	float t_end;
	cout<<"simulation end time (t_end) = ";	cin >> t_end;
	outfile << t_end <<endl;

	float dt;
	cout<<"time step of simulation (dt) = ";	cin >> dt;
	outfile << dt <<endl;

	cout<<endl<<"now generating the random mass , coordinates and velocities..."<<endl;

	/*	generates data 	*/

	for ( int i=0; i<n; i++ )	{
		
		outfile << rand()% 10000 + 5000;	// mass between 5000 - randmax 

		for( int k=0; k < NDIM ; k++ )	{
			outfile << ' ' << 10 * (float)rand() / RAND_MAX;
		}

		for( k=0; k < NDIM ; k++ )	{
			outfile << ' ' << (float) ( rand() % 100) ;
		}

		outfile << endl;
	}
	
	cout<<endl<<endl<< "all datas generated .... ";

	outfile.close();
	return 0;
}



