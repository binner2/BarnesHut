// particle.cpp: implementation of the particle class.
//
//////////////////////////////////////////////////////////////////////

#include "particle.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

particle::particle()
{
	data = 1.0;
	force = 0;
	//pot = 0;
	vel = 0;
	pos = 0;
	Id = 0;
	parent = NULL;

}

particle::particle( double dat , vektor p , vektor v )
{
	data = dat;
	force = 0;
	//pot = 0;
	vel = v;
	pos = p;
	Id = -1;
	parent = NULL;

}


int particle::calcVelPos ( double dt )
{
	// acc = force / m
	vektor acc = force / data;

	vel += acc * 0.5 * dt;
	pos += vel * dt;
	vel += acc * 0.5 * dt;
	return 0;
}

/*particle::~particle()
{
}
*/