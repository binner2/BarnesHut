/*		23 Aðustos 2004						
		version 2							
		|-----> putsnapshot with a message that contains information
				and the values are in scientific format not in fixed format..
		29 nisan 2004						
		version 1.1							
		|-----> readfilepart() with particle class 
		version 1	
		|-----> file.cpp and file.h		
*/

#ifndef _file_h
#define _file_h

extern int particleLEAF;
extern double TETA;

#include "vektor.h"							//for vector 
#include "particle.h"


/********************************************************************************************
**	 GENERATE and write datas to a file like this :
**   n   //number of particles"<<endl;
**   t   //start of simulation time "<<endl;
**   t_end   //end of simulation time"<<endl;
**   dt   //time step of simulation"<<endl;
**   //dt_out   //output time step of simultaion "<<endl;   
**   mass[1]   x[1]   y[1]   z[1]   vx[1]   vy[1]   vz[1]"<<endl;
**   mass , x,y,z coordinates and x,y,z velocities of particles.."<<endl;
********************************************************************************************/
int generatedata( char* filename );

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
int readfilePart( char* filename , long *n , particlePtr *partList , double* starttime , double* endtime , double* stepsize);

/********************************************************************************************
**				put Position of the system to a file..
**   x[1]   y[1]   z[1]   
**   x,y,z coordinates of particles.. so we can plot with gnuplot
********************************************************************************************/
int putsnapshotfilePos( long n , particle* particlePointerList , const char* ,double tetaparam  );

/********************************************************************************************
**				put FORCES of the system to a filename.
**   x[1]   y[1]   z[1]   
**   x,y,z coordinates of particles.. so we can plot with gnuplot
********************************************************************************************/
int putsnapshotfileForce( long n , particle* particlePointerList , const char* );

#endif
