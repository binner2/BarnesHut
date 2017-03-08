/*		
		23 Aðustos 2004						
		version 2							
		|-----> putsnapshot with a message that contains information
				and the values are in scientific format not in fixed format..
		29 nisan 2004						
		version 1.1							
		|-----> readfilepart() with particle class 
		version 1	
		|-----> file.cpp and file.h	

*/

#include "file.h"							//for vector 


/********************************************************************************************
*    this function reads datas from a file with particle class          ************
**   n   //number of particles"<<endl;
**   t   //start of simulation time "<<endl;
**   t_end   //end of simulation time"<<endl;
**   dt   //time step of simulation"<<endl;
**   //dt_out   //output time step of simultaion "<<endl;   
**   mass[1]   x[1]   y[1]   z[1]   vx[1]   vy[1]   vz[1]"<<endl;
**   mass , x,y,z coordinates and x,y,z velocities of particles.."<<endl;
********************************************************************************************/
int readfilePart( char* filename , long *n , particlePtr *partList , double* starttime , double* endtime , double* stepsize)
{
	long i; int k;

   if ( strcmp( filename , "") == 0 ) {	/*check filename for NULL */
      cout<<endl<<endl<<"FILENAME is NULL .."<<endl<<endl;
      exit(3);
   }

   ifstream infile( filename  ) ;

   if ( !infile ) {
      cout<<filename<<" not found in this directory..."<<endl;
	  exit(3);
      return 1;
   }

   infile >> *n;				//nmumber of particles..
   /* control n > 0 */
   if ( *n <= 0  ) {	cout<<"number of particle="<<*n<<" must be bigger then 0";	exit(4);	}

   infile >> *starttime ;

   infile >> *endtime;
   /* endtime must be greater than starttime */
   if ( *endtime < *starttime ) 
   {	cout<<"simulation endtime="<<*endtime<<" must be bigger than starttime="<<*starttime;exit(5); } 

   infile >> *stepsize ;

   /*   read data mass , position and velocities   */
   double mass;

   vektor pos , vel;         //position of n particle...

   *partList = new particle[ *n ];
   if ( partList == NULL ) {	cout<<" bellek yetersiz ";	exit(6);	}

   for ( i=0; i<*n; i++ )   {
      infile >> mass;
	  if ( mass <= 0 ) {	cout<<"mass["<<i<<"]="<<mass<<" must be bigger than zero";	}
		
      for(  k=0; k < NDIM ; k++ )
         infile >> pos[k];

      for( k=0; k < NDIM ; k++ )   
         infile >> vel[k];

	  particle newpart( mass , pos , vel );
	  (*partList)[i] = newpart;

   }//for
    
   return 0;
}

/********************************************************************************************
**				put POSITION of the system to a filename.
**   x[1]   y[1]   z[1]   
**   x,y,z coordinates of particles.. so we can plot with gnuplot
********************************************************************************************/
int putsnapshotfilePos( long n , particle* partList , const char *mesg , double tetaparam )
{
	long i;

	static unsigned int counter = 1;

	char filename[] = "snapPOS__BH100K-tetaCOUNTER.dat";

	sprintf( filename  , "snapPOS__BH%d-teta%1.2f-pLEAF%d-L%d.dat" , n/1000 , TETA , particleLEAF, counter++ );

	ofstream outfile( filename );

	outfile << mesg << " and Number of particle ="
			<<endl<<n<<endl;
	

	for ( i=0; i<n; i++ )   {
  
		outfile << showpos << partList[i].getpos() << endl;

	}//for

	outfile.close();

	return 0;

}

/********************************************************************************************
**				put FORCES of the system to a filename.
**   x[1]   y[1]   z[1]   
**   x,y,z coordinates of particles.. so we can plot with gnuplot
********************************************************************************************/
int putsnapshotfileForce( long n , particle* partList , const char *mesg  )
{
	long i;

	static unsigned int counter = 1;

	char filename[100] = "snapFORCE_BHCOUNTER.dat";

//%s yerine lota yüzünden %d gelmeli..


	sprintf( filename  , "snapFORCE__BH%dK-teta%1.2f-pLEAF%d-Step-%d.dat" , n/1000 , TETA , particleLEAF, counter++ );

	ofstream outfile( filename );

	outfile.precision(40);
	
	outfile << mesg <<endl<<n<<endl;
	
	for ( i=0; i<n; i++ )   {
  
      outfile << showpos << partList[i].getforce() << endl;

   }//for

	outfile.close();

	return 0;

}
