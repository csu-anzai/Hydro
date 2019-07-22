#ifndef LAVA_VEC3_H
#define LAVA_VEC3_H

#include <cstdlib>
#include <ostream>
#include <cstdarg>
#include <iostream>
#include <string.h>

#include "VecTN.h"


/** \class Vec3
    \brief 3 component Vecor class
**/

class Vec3 : public VecTN<double,3>
{
    // protected:
        // double values[3];
    
    public:
        
        /** \brief Default constructor
            
            Sets all values to 0.
        **/
        Vec3()
        {
            this->init_mem();
        }
    
        
        /** \brief Constructor
            
            Sets values to (a,b,c).
        **/
        Vec3(double a, double b, double c)
        {
            this->init_mem();
            this->values[0] = a;
            this->values[1] = b;
            this->values[2] = c;
        }
        
        Vec3(const VecTN<double,3>& rhs)
        {
            this->init_mem();
            // std::cout<<"In copy constructor"<<std::endl;
            for (int i = 0; i < 3; i++)
            {
                // this->values[i] = rhs.get(i);
                rhs.copyValue(i,this->values[i]);
            }
            // return *this;
        }


        Vec3(const VecTN<double,3>* rhs)
        {
            this->init_mem();
            // std::cout<<"In copy constructor"<<std::endl;
            for (int i = 0; i < 3; i++)
            {
                // this->values[i] = rhs.get(i);
                rhs->copyValue(i,this->values[i]);
            }
            // return *this;
        }


        // // Vec3(int nToSet, ... );
        Vec3& operator=(const VecTN<double,3>& rhs)
        {
            //this->init_mem();
            // std::cout<<"In assignment operator (reference)"<<std::endl;
            for (int i = 0; i < 3; i++)
            {
                // this->values[i] = rhs.get(i);
                rhs.copyValue(i,this->values[i]);
            }
            return *this;
            
        }

        Vec3& operator=(const VecTN<double,3>* rhs)
        {
            //this->init_mem();
            // std::cout<<"In assignment operator (pointer)"<<std::endl;
            for (int i = 0; i < 3; i++)
            {
                // this->values[i] = rhs.get(i);
                rhs->copyValue(i,this->values[i]);
            }
            return *this;
            
        }
        
};

#endif
