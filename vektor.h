/*       version 1.1
         added dist() function 
         //scaler distance between two vector
*/

#ifndef _vektor_h
#define _vektor_h

#include "stdinc.h"

class vektor 
{
private:
   double element[3];
public:
   vektor( double c=0 ) { element[0] = element[1] = element[2] = c; }   // constructor

   vektor( double x, double y, double z ) 
   {   element[0] = x; element[1] = y; element[2] = z;   }   //3 arg. constructor.

   
//  []: the return type is declared as a reference (&), so that it can be used
//  on the left-hand side of an asignment, as well as on the right-hand side,
//  i.e.  v[1] = 3.14  and  x = v[2]  are both allowed and work as expected.
   double & operator [] (int i)  {   return element[i];   }

    inline void print() {cout << element[0] << " " << element[1] << " "
              << element[2] << " ";}

//   Unary -

        vektor operator - ()
       {return vektor(-element[0], -element[1], -element[2]);}

//   Dot product.

    	double operator * (const vektor& b)
       {return element[0]*b.element[0]
        + element[1]*b.element[1]
        + element[2]*b.element[2];}

//   Outer product.
	/*
        vektor operator ^ (const vektor &b)
       {return vektor(element[1]*b.element[2] - element[2]*b.element[1],
         element[2]*b.element[0] - element[0]*b.element[2],
         element[0]*b.element[1] - element[1]*b.element[0]);}
	*/
//   vektortor +, -

        vektor operator + (const vektor &b)
       {return vektor(element[0]+b.element[0],
         element[1]+b.element[1],
         element[2]+b.element[2]);}
        vektor operator - (const vektor &b)
       {return vektor(element[0]-b.element[0],
         element[1]-b.element[1],
         element[2]-b.element[2]);}

        friend vektor operator + (double, const vektor & );
        friend vektor operator + (const vektor &, double);

//   Scalar *, /

        friend vektor operator * (double, const vektor & );
        friend vektor operator * (const vektor &, double);
        friend vektor operator / (const vektor &, double);

//   Bool ==, !=: apparently can't do bit-by-bit compare because the
//   member data are private.  Hence these need to be friends in order
//   to apply these operators.  Almost never used, however...

   friend bool operator == (const vektor &, const vektor &);
   friend bool operator != (const vektor &, const vektor &);

//      Post operations (again...) -- the GCC3 december 2002 sage:
//      this is a some interesting way to resolve the ambiguity of
//      resolving  v*s vs. s*v and the confusion of having a v*v with
//      a non-explicit vektor constructor from a scalar.....

   vektor operator * (const double b)
     {return vektor(element[0]*b, element[1]*b, element[2]*b); }

   vektor operator + (const double b)
     {return vektor(element[0]+b, element[1]+b, element[2]+b); }

//   vektortor +=, -=, *=, /=

        vektor& operator += (const vektor& b)
       {element[0] += b.element[0];       
        element[1] += b.element[1];
        element[2] += b.element[2];
        return *this;}

   vektor& operator -= (const vektor& b)
       {element[0] -= b.element[0];
        element[1] -= b.element[1];
        element[2] -= b.element[2];
        return *this;}

   vektor& operator *= (const double b)
       {element[0] *= b;
        element[1] *= b;
        element[2] *= b;
        return *this;}

   vektor& operator /= (const double b)
       {element[0] /= b;
        element[1] /= b;
        element[2] /= b;
        return *this;}

//      Input / Output

   friend ostream& operator << (ostream & , const vektor & );
   friend istream& operator >> (istream & , vektor & );
};

inline  ostream & operator << (ostream & s, const vektor & v)
       {return s << v.element[0] << "  " << v.element[1]
            << "  " << v.element[2];}

inline  istream & operator >> (istream & s, vektor & v)
       {s >> v.element[0] >> v.element[1] >> v.element[2];
        return s;}

inline  double square(vektor v) {return v*v;}
//inline  double abs(vektor v)    {return sqrt(v*v);}

// Another measure of vektortor magnitude; less work than abs():

//inline  double abs1(vektor v)   {return abs(v[0]) + abs(v[1]) + abs(v[2]);}

inline  vektor operator + (double b, const vektor & v)
       {return vektor(b+v.element[0],
         b+v.element[1],
         b+v.element[2]);}

inline  vektor operator + (const vektor & v, double b)
       {return vektor(b+v.element[0],
         b+v.element[1],
         b+v.element[2]);}

inline  vektor operator * (double b, const vektor & v)
       {return vektor(b*v.element[0],
         b*v.element[1],
         b*v.element[2]);}

inline  vektor operator * (const vektor & v, double b)
       {return vektor(b*v.element[0],
         b*v.element[1],
         b*v.element[2]);}

inline  vektor operator / (const vektor & v, double b)
       {return vektor(v.element[0]/b,
         v.element[1]/b,
         v.element[2]/b);}

inline  bool operator == (const vektor& u, const vektor& v)
       {return (u.element[0] == v.element[0]
           && u.element[1] == v.element[1]
           && u.element[2] == v.element[2]);}

inline  bool operator != (const vektor& u, const vektor& v)
       {return !(u==v);}

/*************************************************************************
**      scaler distance between two vector                         ** 
**                                                       **
**************************************************************************/
static double dist( vektor vec1 , vektor vec2 )
{   
   double distance=0.0;
   for ( int k=0;k < NDIM ; k++ )
      distance += (vec1[k]-vec2[k]) * (vec1[k]-vec2[k]);

   return distance;
}

#endif
