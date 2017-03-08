/****************************************************************************
* STDINC.H: standard include inspired from J. Barnes nbody code..
* This file contains the standart includes
* Burak ÝNNER , 5-7-2004                                                          
****************************************************************************/

#ifndef _stdinc_h
#define _stdinc_h

#include <iostream>
#include <fstream>                     //for file I/O
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>

#define NDIM 3                  //Number of dimension for problem.
#define NSUB 1<<NDIM            //number of child nodes.

//#define gravity 6.672e-1
#define gravity 1               //gravity or coulomb constant.

#define showtime 1

#define xdebug 1

using namespace std;

#endif
