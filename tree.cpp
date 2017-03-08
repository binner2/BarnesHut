/************************************************************************************
**		tree.cpp -- FAST MULTIPOLE METHOD tree routines..
added createRepo()	initialize the Repository with 3.BlockSize
added mynew()	use a Node array..
***********************************************************************************/
#include "tree.h"

/************************************************************************************
**		clearNode()
**		clear the content of the Node 
************************************************************************************/
int tree::clearNode ( Node* oldNode )
{
	int k;
/************* extra debugging ***********************************************************/
#ifdef xdebug    		
	if ( oldNode == NULL ) { cout<<"clearNode için null root ";return -1;}
#endif    
/*****************************************************************************************/
	
	for ( k=0; k < NSUB ; k ++ ) 
		oldNode->child[ k ] = NULL;
	
	oldNode->geocenter = 0;

	oldNode->index = 0;
	
	oldNode->level = 0;

	oldNode->mass = 0;

	oldNode->masscenter = (0);

	oldNode->Nparticle = 0;

	//clear the pointers to particles in this LEAF (NODE dont have any particle pointer)
	oldNode->plist.clear();

	oldNode->rsize = 0;

	oldNode->Type = EMPTY;	

	oldNode->parent = NULL;
	
	return 0;
}

/************************************************************************************
**		clearRepo, allocate and initialize the Node's array.. 
**		it's a two dimensional array..
**		Nodelist[index1][index2]
************************************************************************************/
void tree::createRepo ( long N )
{
	currentNode=0;
	Nparticle = N;
}

/************************************************************************************
**		mynew() return a free Node from the array.
**		
************************************************************************************/
Node* tree::mynew( void )
{
   if ( currentNode < (long)Nodelist.size() )
   {
	   clearNode( &Nodelist[currentNode] );
      return ( &(Nodelist[currentNode++]) );
   }
   
   Node *returnthisNode = new Node;

   clearNode( returnthisNode );
   
   Nodelist.push_back(*returnthisNode);

   currentNode++;

   return returnthisNode;

}


/************************************************************************************
**		clearNodelist() reset the Node's 
**		
************************************************************************************/
void tree::clearNodelist()
{
	currentNode = 0;
}

/************************************************************************************
**		emptytree release the memory for all NODES 
**		free the tree memory 
************************************************************************************/
void tree::emptytree( char *msg )
{

sprintf(msg , "DirectForce;%d;ParticleCell;%d;Usednode;%d;availableNode;%d", countdirectforce,countparticlecell,currentNode,(int)Nodelist.size() );

	countdirectforce   = 0;

	clearNodelist();
	
	clearNode( BHroot );

	BHroot->index = 0;
	
	BHroot->geocenter = ( 50 );

	BHroot->rsize = 100;
}

//////////////////////////////////////////////////////////////////////
// TREE Construction
//////////////////////////////////////////////////////////////////////
tree::tree( particlePtr partList , double stepsize=0.01 , long N=1 ) 	
{	
	treepartList = partList;		//particle list array 

	dt=stepsize;					//stepsize of simulation needed for force calculation  ??? 

	nlevel = 0;						//number of level of tree

	Nparticle = N;
	
	BHroot = new Node;				//root of the Barnes-Hut tree = BHroot
	clearNode ( BHroot );

	//default values
	BHroot->mass = Nparticle;
	BHroot->geocenter = 0;
	BHroot->rsize = 0;
	//geocenter and rsize will be replaced by findroot function

	createRepo( Nparticle );

	countdirectforce   = 0;
	countparticlecell=0;

	for (int i=0;i<Nparticle;i++)	
		treepartList[i].Id = i;


}


/************************************************************************************
**		load particles into tree , start from BHroot node use the particle array partList
**		and the number of particles N .
************************************************************************************/
void tree::loadparticles( )
{
	long i;

	findroot( Nparticle , treepartList , BHroot->rsize , BHroot->geocenter );

	int j;
	bool isinserted = false;
	Node *root;

	//insert particle into the root cell one by one..

	for ( i=0; i<Nparticle; i++)
	{
		root = BHroot;

		while ( isinserted == false )
		{
			j = addparticle( i , &treepartList[i] , root );
			if ( j == -1 ) isinserted = true;
			else if ( j>=0 && j< NSUB ) root = root->child[j];
			else cout<<"return value of addparticle must be between -1 and 8..";

		}//while

		isinserted = false;

	}
}

/************************************************************************************
**		add particles into tree recursively , start from root node 
**		
************************************************************************************/
int tree::addparticle ( long index , particle *part , Node* root )		
{			
	//particle will place which child of the root ??
	int wcld = whichchild( part->getpos() , root );

	Node *childcell = root->child[ wcld ];

	if ( childcell == NULL || childcell->Type == EMPTY ) 
	{
		addleaf( index  , part , root , wcld );
		root->Nparticle++;

		//leaf's parent is root..
		root->child[wcld]->parent = root;

		return -1;//particle inserted.. added.
	}
	else if ( childcell->Type == LEAF )
	{
	
		//look for # of particle in LEAF
		if ( childcell->Nparticle < particleLEAF )
		{
			childcell->Nparticle++;
			childcell->plist.push_back(part);   //add a pointer to this particle
			
			part->parent = childcell;
			root->Nparticle++;
			return -1;							//particle inserted return ok.
		}
		else //# of particle in LEAF is bigger then NPART
		{
			leaf2node( root , wcld );			//convert root's wcld child to a node.
			root->Nparticle++;
			return wcld;
		}
	}

	else if ( childcell->Type == NODE )
	{
		root->Nparticle++;
		return wcld;
	}
	cout<<"there is a problem in addparticle";
	return -2;
}


/************************************************************************************
**		which child of the root NODE will place the particle with position pos 
**		and return the number of child node
************************************************************************************/
int tree::whichchild( vektor pos , Node* root )
{
#ifdef xdebug 
	if ( root == NULL ) 
	{
		cout<<"In whichchild function root is a NULL value..";
		return -2;
	}
#endif

	int childnumber=0;

	for ( int k=0; k<NDIM; k++ )
	{
		if (pos[k] >= root->geocenter[ k ] )
		{
			childnumber += 1<<k;
/************* extra debugging ****************************************************************************/
#ifdef xdebug 
			if ( pos[k] > (root->geocenter[k] +	root->rsize / 2)   )
			{
				cout<<"sorun var .. "<<pos[k] << " root " << root->geocenter[k]<<" rsize="<<root->rsize<<endl;;
				return -2;
			}

		}
		else if ( pos[k] < ( root->geocenter[k] -	root->rsize / 2)   )
		{
			cout<<"pos root içinde deðil.." << pos[k] << " root " << root->geocenter[k]<<endl;
			return -2;
#endif
/************* extra debugging ****************************************************************************/

		}
	}//for

	if ( childnumber>=0 && childnumber<8 ) 
		return childnumber;
	else 
	{
		cout<<"childnumber'da sorun var = "<<childnumber;
		return -2;
	}
	
	
}

/************************************************************************************
**		add leaf to the NODE root 
**		must to know which child of the root..
************************************************************************************/
void tree::addleaf( long index , particle* part , Node* root , int wcld)
{

/************* extra debugging ****************************************************************************/
#ifdef xdebug 
	// insert a leaf to root's wcld child the particle.	
	if ( root->Type == LEAF )
	{
		cout<<"\naddleaf LEAF node 'a ekleme yapýyor.";
		display( root );
	}
#endif
/************* extra debugging ****************************************************************************/

	root->Type = NODE;					//root was a EMPTY and now it is a NODE

	Node* newLeaf = mynew();

	newLeaf->level = root->level+1;

	newLeaf->Nparticle++;
	newLeaf->plist.push_back(part);

	part->parent = newLeaf;


	newLeaf->rsize = root->rsize / 2.0;
/************* extra debugging ****************************************************************************/
#ifdef xdebug 
	if ( newLeaf->rsize == 0 ) 
		cout<<"add leaf icinde rsize hesaplamasýnda 0 var.. L="<<newLeaf->level<<" index="<<newLeaf->index;
#endif
/************* extra debugging ****************************************************************************/
	
	newLeaf->Type = LEAF;

	/* calculate the center of the leaf cube..	*/
	for ( int k =0; k < NDIM ; k++ ){
         if ( ( wcld >>k) % 2 ) 
			 newLeaf->geocenter[k] = root->geocenter[k] + newLeaf->rsize/2.0;
         else newLeaf->geocenter[k] = root->geocenter[k] - newLeaf->rsize / 2.0;
	}
	
	/* newLeaf ok */ 
	
	root->child[wcld] = newLeaf;

	//number of level of the tree..
	if ( nlevel < newLeaf->level ) nlevel = newLeaf->level; 

}

/************************************************************************************
**		convert leaf to a NODE  root->child[ wlcd ] is a leaf convert to node.
**		convert root nodes wcld'th number of children leaf to a node
**		add the old leaf to the new node as a leaf..
************************************************************************************/
void tree::leaf2node( Node* root , int wcld )		
{
  Node* oldleaf = root->child[ wcld ];
   
   int Nparticle2 = oldleaf->Nparticle;
   
   oldleaf->Nparticle = 0;

   vector<particle*> tempplist = oldleaf->plist;

   oldleaf->Type = NODE;

   vector<particle*>::const_iterator temppart; 

   for (temppart = oldleaf->plist.begin(); temppart != oldleaf->plist.end(); ++temppart) 
   {
      wcld = whichchild( (*temppart)->getpos() , oldleaf );
      addparticle( (*temppart)->getId() , *temppart , oldleaf );

   }

   oldleaf->plist.clear();

}

/*************************************************************************
**      Upward PASS of the tree 										** 
**		Calculate multipole expansion of particles and					**
**		shift multipole expansion to parent node 
**************************************************************************/
void tree::upwardpass( void )
{
	//calculate Center Of Mass recursively start from Barnes-Hut Tree's root Node.
	calCMS( BHroot );
} 

/*************************************************************************
**		Calculate multipole expansion of particles and					**
**		shift multipole expansion to parent NODE						** 
**	if root is a LEAF calculate ME of LEAF parent			 			**
**************************************************************************/
void tree::calCMS( Node* parent ) 
{
	int i,k;
	Node* childNode;

   vektor  cmsNode =(0);         //center of mass of this cube..
   double  summass=0.0;         //sum of the children masses.

	vektor cmsNodeLEAF = (0);
	double summassLEAF = 0.0;

	
	// for all children of parent NODE 
	for( i=0; i<NSUB; i++)				
	{
		// childNode is the child of parent
		childNode = parent->child[i];

		// always control for NULL
		if ( childNode != NULL  && childNode->Type!= EMPTY )
		{

			if ( childNode->Type == LEAF )
			{
				cmsNodeLEAF = (0);
				summassLEAF = 0;

				for ( k=0;k<childNode->Nparticle;k++)
				{
					cmsNodeLEAF += childNode->plist[k]->data * childNode->plist[k]->pos;
					summassLEAF += childNode->plist[k]->data;
				}

				childNode->mass = summassLEAF;
				childNode->masscenter = cmsNodeLEAF / summassLEAF;

				cmsNode += cmsNodeLEAF;
				summass += summassLEAF;

				cmsNodeLEAF = 0; summassLEAF = 0;
			}//if
			

			/*if its a NODE, than recursively calculate all childrens Multipole expansion
			and shift all children multipole expansion's to parent cell */
			else if ( childNode->Type == NODE )
			{
				calCMS( childNode );		//recursive calc.

				cmsNode += childNode->mass * childNode->masscenter;
				summass += childNode->mass;

			}//else if
		}//if NULL
	}//for

   cmsNode = cmsNode / summass;

   parent->masscenter = cmsNode;
   parent->mass = summass;


}



/*************************************************************************
**      Downward PASS of the tree 										** 
**		Calculate local expansion for cells and							**
**		shift local expansion to parent node							**
**************************************************************************/
void tree::downwardpass( void )
{
	//cout<<"in downwardpass nothing to do.."<<endl;
}

/************************************************************************************
**		for calculating forces, first we need to calculate local expansions,
**		than with interaction list top-down traversal tree and calculate children local
**		expansion.
**		
************************************************************************************/
void tree::calculateforces()
{
	int i=0 , k=0;

	//control the ID of all particles.ID is the index of the particle array.
	for (i=0;i<Nparticle;i++)
	{
		treepartList[i].force = 0;
/************* extra debugging ****************************************************************************/
#ifdef xdebug 
		if ( treepartList[i].Id != i ) cout<<"ID lerde sorun var. i = " << i << endl;
#endif
/************* extra debugging ****************************************************************************/
	}

/************* time ****************************************************************************/
#ifdef showtime    		
	clock_t start = clock();
#endif    
/*****************************************************************************************/


	double potential=0;
	vektor force=0;

	//for all particles look for the local expansion and near forces contribution..
	for ( i=0; i<Nparticle;i++)
	{
		//root cannot interact any particle..
		for(k=0;k<NSUB;k++)
			if ( BHroot->child[k] != NULL && BHroot->child[k]->Type!= EMPTY )
				interact( &(treepartList[i]) , BHroot->child[k] );
	}//for i.

	////////////////////////////////////////////////////////////////////////
	// now calculate the near neighbour cells particle directly.
	//////////////////////////////////////////////////////////////////////
	
/************* time ****************************************************************************/
#ifdef showtime    		
	clock_t end = clock();
	cout<<"\ncalculate force take"<< ( end-start) / (float)CLOCKS_PER_SEC<<" seconds.."<< endl;

	/*when all forces calculated then calculate acceleration, velocities and positions.*/

	start = clock();
#endif    
/*****************************************************************************************/

	calcVelPos();

/************* time ****************************************************************************/
#ifdef showtime    		
	end = clock();
	cout<<"\ncalculate velocities and position from force take"<< ( end-start) / (float)CLOCKS_PER_SEC<<" seconds.."<< endl;
#endif    
/*****************************************************************************************/

/************* extra debugging ****************************************************************************/
#ifdef xdebug 
		cout<<"direct force ="<<countdirectforce<<"usage of node ="<<currentNode<<" available node ="<<(long)Nodelist.size()<<endl;
#endif
/************* extra debugging ****************************************************************************/

	//displaytree( BHroot );
}


/************************************************************************************
**	wellseperate ??
*************************************************************************************/
bool tree::wellseperate( particle thisparticle , Node othercell )
{
	double r = dist ( thisparticle.pos , othercell.masscenter );
    r = sqrt( r );	

	if ( othercell.rsize / r <= TETA )
	   return true;
   else 
	   return false;
}


/************************************************************************************
**		interact this cell with the thiscell's parent's near neighbours
**	neparent = is the near neighbour of thiscell 's parent so look for all children of
**  neparent and interact with or add near neighbours
************************************************************************************/
void tree::interact( particle *thisparticle , Node* otherCELL )
{
	int i=0;
	Node *childNode;
	//double rkare , rkup,r;
	if ( wellseperate( *thisparticle , *otherCELL ) == true ) 
	{
		if ( otherCELL->Type == NODE )
		{
			particlecell( thisparticle , otherCELL );
		}
		else if ( otherCELL->Type == LEAF ) 
		{
			particlecell( thisparticle , otherCELL );
		}//if LEAF
		else if ( otherCELL->Type == EMPTY ) 
			cout<<"there is a problem in interact because otherCELL is empty"<<endl;

	}
	else		//not wellseperated..
	{
		if ( otherCELL->Type == NODE )
		{
			for( i=0; i<NSUB; i++)            
			{
				childNode = otherCELL->child[i];
				if ( childNode != NULL  && childNode->Type!= EMPTY )
				{
					interact( thisparticle , childNode );				
				}//if !NULL
			}//for i
		}//if
		else if ( otherCELL->Type == LEAF ) 
		{
			for( i=0;i<otherCELL->Nparticle;i++)
				directforcalc( thisparticle , otherCELL->plist[i] );

		}//if LEAF
		else if ( otherCELL->Type == EMPTY ) 
			cout<<"there is a problem in interact because otherCELL is empty"<<endl;

	}//else
}


/************************************************************************************
**      calculate forces between particle and cell ..
**      ...
************************************************************************************/
void tree::particlecell( particle *part , Node* cell )
{
	countparticlecell++;
   double rkare = dist( part->pos , cell->masscenter );
   double rkup = rkare * sqrt( rkare );
   vektor r = part->pos - cell->masscenter;

   part->force += -gravity * part->data * cell->mass / rkup * r;
}

/************************************************************************************
**		calculate forces between two particle ..
**		particle particle method or Direct method...
************************************************************************************/
void tree::directforcalc( particle *part1 , particle *part2 )
{
	/*control that index is not equal index2*/
	if ( part1->Id == part2->Id ) return;

	vektor r = part1->getpos() - part2->getpos();

	double rkare = dist( part1->getpos() , part2->getpos() );

	double rkup = rkare * sqrt( rkare );

/************* extra debugging ****************************************************************************/
#ifdef xdebug    		 
	/*not Necessary but want to control collision */
	if ( rkare == 0 ) 
	{
		cout<<endl<<"fatal error distance between two particle is 0 " ;		cout<<endl;
		part1->getpos().print();		part2->getpos().print() ;
		return;
	}
#endif
/*******************************************************************************************************/

	vektor force = part1->getforce();

	//				      M1 . M2 
	// (vektor ) F = G . --------- . (vektor) r
	//			           rkup

	force += -gravity * part1->getdata() * part2->getdata() / rkup * r;

	part1->setforce( force );

	countdirectforce++;
}


/************************************************************************************
**		calculate velocities and positions..
**		...
************************************************************************************/
void tree::calcVelPos( )
{
	/*all treeAcc is updated.*/
	for ( int i=0; i<Nparticle ; i++)
	{
		treepartList[i].calcVelPos ( dt ) ;
	}//for

	//displaytree( BHroot );

}


/*************************************************************************
**      Find the max and min coordinate for the root cell            ** 
**      calculates the max length of one side of the root cells.      **
**  root cell coordinates vektor root , and length of one side is rsize
**************************************************************************/

void tree::findroot ( const int n , particle pos[] , double &rsize  , vektor &root )
{
   vektor maxpos , minpos , partpos;
   int i , k;

   minpos = maxpos = pos[0].getpos();            //overloaded operator..
   
   // find the min. and max koordinate of the system
   for ( i=1; i<n; i++ )   {
      partpos = pos[i].getpos();
      for ( k=0;k<NDIM;k++)
      {
         if ( maxpos[k] < partpos[k] ) maxpos[k] = partpos[k];
         if ( minpos[k] > partpos[k] ) minpos[k] = partpos[k];
      }
   }

   /* needed to calculate the length of the root cubes one side ( rsize )  
   rsize is the max. distance between the max position and min position.
   but be careful to control all dimension one by one independent.
   */
   
   double maxdist[ NDIM ];

   for ( k=0;k< NDIM ; k++){
      maxdist[k] = maxpos[k] - minpos[k];
		//root cell coordinates.
      root[k] = minpos[k] + ( maxdist[k] / 2 );
	  //root[k] = minpos[k] + maxpos[k] / 2;
   }
   
   /* find the rsize , max distance between the max. and min.*/

   rsize = maxdist[0];
   for ( k=1; k<NDIM ; k++) 
      if ( rsize < maxdist[k] ) 
         rsize = maxdist[k];

 	  if ( rsize == ceil( rsize ) )
		  rsize++;
	  else 
		  rsize = ceil( rsize );

/************* extra debugging ****************************************************************************/
#ifdef xdebug           
   cout<<"\ncoordinate of center of the root cell= " << root
      <<endl<<"and length of one side =" << rsize << endl
      <<"region =";
   for ( k=0; k<NDIM; k++)
      cout<<"("<<root[k]-rsize/2<<","<<root[k]+rsize/2<<")";
   cout<<endl;
#endif    
/*****************************************************************************************/

   

}


/************************************************************************************
**		display root NODE or LEAF 's all information 
**
*************************************************************************************/
void tree::display( Node* root )
{
	cout.setf( ios::right );
	cout.setf( ios::fixed );
	cout.setf( ios::showpoint );
	cout.setf( ios::showpos );
	cout<< setw(4) ;
	cout<< setprecision( 2 ) ;
	
	cout<<" Id="<<root->index;
	cout<<" l="<<root->level;
	cout<<" m="<<root->mass;
	cout<<" Npar="<<root->Nparticle;
	cout<<" geo="<<root->geocenter;
	cout<<" rsize="<<root->rsize;
	cout<<"\tcms="<<root->masscenter;
	if ( root->Type == NODE )
		cout<<"\tType=NODE"<<endl;
	else if ( root->Type == EMPTY ) 
		cout<<"\tType=EMPTY "<<endl;
	else if ( root->Type == LEAF ) 
	{
		cout<<"\tType=LEAF"<<endl;
		for ( int i=0; i<root->Nparticle; i++ )
		{
			
			cout<<i+1<<". particle ID= "<<root->plist[i]->getId()<<" - ";
			root->plist[i]->display();
			cout<<endl;
		}
	}
	//cout<<endl;

}

/************************************************************************************
**		display all tree start from root NODE toward LEAVES
**
*************************************************************************************/
void tree::displaytree( Node* root ) 
{
	int i;
	if ( root == NULL ) {
		cout<<"NULL geldi .. " <<endl;
		return;	
	}

	if ( root->Type == NODE ) 
	{
		display( root );

		for( i=0; i<NSUB; i++)				{
			if ( root->child[i] == NULL || root->child[i]->Type == EMPTY ) 
				cout<<" index = "<<i<<" level="<<root->level+1<<" is AN EMPTY or NULL NODE"<<endl;
			else //if ( root->child[i] != NULL )
			{
				displaytree( root->child[i] );
			}//if 

		}//for
	}//if
	
	// if root is a LEAF
	else if ( root->Type == LEAF ) 
		display( root );

	else if ( root->Type == EMPTY ) 
		cout<<" empty node = "<<root->index << " " << root->level;

  
	//cout<<"\nlevel of the tree = " << nlevel;

}
