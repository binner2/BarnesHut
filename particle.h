// particle.h: interface for the particle class.
//
//////////////////////////////////////////////////////////////////////


#ifndef _particle_h
#define _particle_h

#include <vector>      //STL vector.
#include "vektor.h"      //3d vector.

class particle;

enum _type { EMPTY=0 , NODE=1 , LEAF=2 };

struct Node
{
	long index;						//index of the NODE array.

	_type Type;						//type of this NODE or LEAF

	vektor geocenter;				//geometric center coordinate of the node

	double rsize;					//length of the cubes one side.

	vektor masscenter;

	double mass;

    vector<particle*> plist;		// list of particles in cell IF LEAF

	int Nparticle;			//number of particle in this node..

	long level;			//level of the tree
	
	Node* child[NSUB];				//link to children

	Node* parent;

};

// Reset The Node all variables to default value.  Type = EMPTY 


class particle  
{
private:
	double data;		//mass or q
	vektor pos;
	vektor vel;
	//double pot;
	vektor force;
	unsigned int Id;
	Node *parent;

public:

	friend class tree;

	particle();
	particle( double data , vektor p , vektor v );

	//virtual ~particle();

	double getdata() { return data; };
	vektor getpos()  { return pos; };
	vektor getvel()  { return vel;};
	//double getpot()  { return pot;};
	vektor getforce(){ return force;};
	unsigned int getId() { return Id; };
	
	void setpos( vektor p ) { pos = p; };
	void setvel( vektor v ) { vel = v; };

	void setforce( vektor f ) 
	{ 
		force = f; 
		//cout<<"particle içindeki force "<<force;
	};

	//void setpot ( double potential ) { pot = potential; };
	void setId( unsigned int Identity ) { Id = Identity; };
	
	int calcVelPos ( double dt ) ;

	void display() { 
		cout<<"parent="<<parent->index;
		cout<<" data="<<data<<" force=";
		force.print();
		cout<<" -pos="; 
		pos.print(); 
	};
};

typedef particle *particlePtr;



#endif 
