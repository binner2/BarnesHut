#ifndef _tree_h
#define _tree_h

#include <vector>
#include "particle.h"
#include "vektor.h"

extern double TETA;
extern int particleLEAF;

class tree
{
private:

	Node *BHroot;					//root of the BH tree...
	particle* treepartList;			//particle array holds the particle information.
	double dt;						//timestep = dt example = 0.001 or 0.1 
	long nlevel;			// # level of the tree
	long Nparticle;

	long currentNode;					//current Node in the array.
	vector<Node> Nodelist;

	long countdirectforce;
	long countparticlecell;

public:
		
	friend class particle;

	/* Tree construction */
	tree( particlePtr partList , double stepsize , long N  );

	//we need this because we can't access BHroot from main()
	void emptytree( char *msg );

	Node* getroot( void ) { return BHroot; };		//get root node of the tree..

	//find the length of one side of the root node that encloses all the particle..
	void findroot ( const int n , particle pos[] , double &rsize  , vektor &root );

	//loadparticles with particle array and number of particle = N
	void loadparticles();

	//add particles to the tree one by one with the index that show the indix of the particle array.
	int addparticle ( long index , particle* part , Node* root) ;

	//find the number of child of the root NODE's where the particle (part ) will be inserted.
	int whichchild( vektor part , Node* root );
	
	//convert the leaf to Node but be careful when LEAF has more then 1 particle.
	void leaf2node( Node* root , int wcld );

	//add a LEAF to the root's wcld child ..
	void addleaf( long index , particle *part , Node* root , int wcld);

	//upwarad pass of the tree.
	void upwardpass( );

	void calCMS( Node* parent ); 
	
	//downward pass of the tree.
	void downwardpass( );

	//convert multipole expansion to local expansion.
	void interact( particle *thisparticle , Node* otherCELL );

	//if cell's well seperated return true..
	bool wellseperate( particle thisparticle , Node othercell );

	void particlecell( particle *part , Node* cell );

	//calculate forces in LEAF cells but near neighbour will be calculated directly.
	void calculateforces();

	//near neighbour will be calculated directly.
	//void calforcenear ( long index , Node* parent ) ;

	//direct force calculation with two particle index and index2 are the index of particle array.
	void directforcalc( particle* part1 , particle *part2 );

	//display whole the tree..
	void displaytree( Node* root ) ;

	//dislpay the Node information.
	void display( Node* root );
	
	//calculate velocities and position..
	void calcVelPos();

	int clearNode ( Node* oldNode );

	void createRepo( long N );

	Node* mynew( void );

	void clearNodelist();

};

#endif

