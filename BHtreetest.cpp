/*-------------------------------------------------------------------------
-- treetest.cpp - Fast Multipole Method N-body implementation for coulomb forces..
-- Burak ÝNNER , 6-7-2004
	01-03-2005
		teta and particleLEAF global variable.
		LEAF has a STL::vector for particle pointers.
		Nodelist is a STL::vector for NODE pointers.
	07-02-2005
		review..
	26-8-2004
	-->TETA parameter NOT defined as static double in stdinc.h 
	instead of static, tetaparameter send as a function argument and
	defined in TREE
---------------------------------------------------------------------------*/

#include "file.h"
#include "tree.h"

double TETA;
int particleLEAF;

int main( int argc, char *argv[] )
{
	int i=0;

	long N = 0;				//number of particle...
	particlePtr partList = NULL;		//particle array..
	double starttime=0.0;
	double endtime=0.0;
	double stepsize=0.0;
	clock_t start,end;
	char message[300]=" ";
	//char timetakes1[200]=" ";
	double timeLoad , timeUpward , timeForce;

/****************** read from filename mass velocities position ************************/
	char filename[25] = "filenamedata1.dat";
	if ( argc == 4 ) {
		strcpy( filename , argv[1] );
		cout<<"reading particle information from file " <<filename<<endl;
		TETA = atof(argv[2]);
		particleLEAF = atoi( argv[3] );

		//generatedata2( filename );
		readfilePart( filename , &N , &partList , &starttime , &endtime , &stepsize);

	}
	else
	{
		cout<<"\nparameter problem..\nbh N TETA pLEAF\nN is the number of particles, pLEAF is number of particles in LEAF\n";
		exit(2);
	}

/************* time ****************************************************************************/
#ifdef showtime    		
	start = clock();
#endif    
/*****************************************************************************************/

	cout <<"Starting a BARNES HUT algorithm with Lepfrog integration for a \n" << N
     << "-body system,from time t =" << starttime <<" to t=" << endtime <<" with time step = " << stepsize
	 << " \nand number of max. particle in a LEAF " << particleLEAF << " TETA="<< TETA << endl;

/************* time ****************************************************************************/
#ifdef showtime    		
	end = clock();
	cout<<"read from file elapsed "<< (end-start) /(float)CLOCKS_PER_SEC<<" seconds.."<< endl;
#endif    
/*****************************************************************************************/

	clock_t onestep_time_start;				//one step time of the simulation time..

	/* construct tree with the simulation parameter */
	tree bhtree( partList , stepsize , N );
	
/****************************************************************************************/

	int counter=0;			//number of total step.

	double dt_out= (endtime - starttime) / 10;	//snapshot output frequency.
	double t_out = starttime + dt_out;			//put snapshot on t_out
	
	while ( endtime > starttime  )
	{
		onestep_time_start = clock();

		strcpy( message , " " );
		
/************* time ****************************************************************************/
#ifdef showtime    		
		start = clock();
#endif    
/*****************************************************************************************/

		bhtree.loadparticles( );

/************* time ****************************************************************************/
#ifdef showtime    		
		end = clock();
		timeLoad = (end-start)/(float)CLOCKS_PER_SEC;
		cout<<"\nload particles takes "<< timeLoad <<" seconds"<< endl;
		start = clock();
#endif    
/*****************************************************************************************/

		bhtree.upwardpass();	//	upward pass .. calculate center of mass for the cubes..

/************* time ****************************************************************************/
#ifdef showtime    		
		end = clock();
		timeUpward = (end-start)/(float)CLOCKS_PER_SEC;
		cout<<"\nupward pass ( without load particles) takes "<< timeUpward <<" seconds"<< endl;
		start = clock();
#endif    
/*****************************************************************************************/

		//bhtree.downwardpass();		// downward pass ... 

/************* time ****************************************************************************/
#ifdef showtime    		
		end = clock();
		//cout<<"\ndownward pass ( without load particles with upwardpass) takes "<< ( end-start) / (float)CLOCKS_PER_SEC<<" seconds"<< endl;
		start = clock();
#endif    
/*****************************************************************************************/

		bhtree.calculateforces();		//	calculate forces with the tree..

/************* time ****************************************************************************/
#ifdef showtime    		
		end = clock();
		timeForce = ( end-start) / (float)CLOCKS_PER_SEC;
		cout<<"\ncalculate forces take ( without load particles.."<< timeForce <<" seconds.."<< endl;
#endif    
/*****************************************************************************************/
		
		bhtree.emptytree( message );

		starttime += stepsize;

		counter++;

/************* time ****************************************************************************/
#ifdef showtime    		
	clock_t onestep_time_end=clock();
	cout<<"\nall takes "<< ( onestep_time_end-onestep_time_start) / (float)CLOCKS_PER_SEC<<" seconds"<< endl;
	char timetakes[300]=" ";
	sprintf(timetakes,"BH;pLEAF;%d;step;%d;teta;%1.2f;timeSort;-;timeLoad;%f;timeupward;%f;timeDownward;;timeforce;%f;all;%f;",particleLEAF, counter , TETA, timeLoad , timeUpward, timeForce , ( onestep_time_end-onestep_time_start) / (float)CLOCKS_PER_SEC );
	strcat( timetakes , message );
	cout<<"number of step = "<<counter;
#endif    
/*****************************************************************************************/

		///////////////////////////////////////////////////////////////////////////////
		// now it write file the forces array..
		putsnapshotfileForce( N , partList , timetakes  );

	}

/************* time ****************************************************************************
#ifdef showtime    		
	clock_t onestep_time_end=clock();
	cout<<"\nall takes "<< ( onestep_time_end-allstart) / (float)CLOCKS_PER_SEC<<" seconds"<< endl;
	char timetakes[300]=" ";
	sprintf(timetakes,"BH;teta;%f;pLEAF;%d;all;%f;",TETA, particleLEAF,  ( onestep_time_end-allstart) / (float)CLOCKS_PER_SEC );

	cout<<"number of step = "<<counter;
#endif    
********************************************************************************************/

	///////////////////////////////////////////////////////////////////////////////
	// now it write file the forces array..
	//putsnapshotfileForce( N , partList , timetakes  );
		
	//cout<<endl<<"end of the BARNES-HUT"<<endl;
	return 0;
}