#ifndef LAVA_VecTN_H
#define LAVA_VecTN_H

#include <cstdlib>
#include <ostream>
#include <cstdarg>
#include <iostream>
#include <string.h>
#include <cmath>


/** \class VecTN
    \brief 3 component Vector class
**/
template <class T, int N>
class VecTN
{
    protected:
        T values[N];

    protected:
        void init_mem();
        // {
        //     memset(&values[0], 0, N*sizeof(T));
        // }

    public:
        
        /** \brief Default constructor
            
            Sets all values to 0.
        **/
        VecTN();
    
        
        /** \brief Constructor
            
            Sets values to (a,b,c).
        **/
        // VecTN(T a, T b, T c)
        // {
        //     values[0] = a;
        //     values[1] = b;
        //     values[2] = c;
        // }
        
        VecTN(int nToSet, ... );
        
        /** \brief Returns value at index
            
            Overloaded [] operator returns value at index i
            @param i Index of Vecor
        **/
        T& operator[](const int i);

        T& get(const int i);
        void copyValue(const int i, T& retVal) const;
      
        /** \brief Prints Vecor values
            
            Values are written tab delimited with an endl.
            @param out ostream to write to
        **/
        void print(std::ostream& out);
        
        /** \brief Dot product
        **/
        VecTN<T,N> operator*(VecTN<T,N> b);
        
        /** \brief Normalizes the Vecor
        **/
        VecTN<T,N>& normalize();

        double norm();
        
        
        /** \brief Adds 2 Vecors
        **/
        VecTN<T,N> operator+(VecTN<T,N> b);
        VecTN<T,N> operator+(VecTN<T,N>* b);
        
        
        /** \brief Subtracts 2 Vecors
        **/
        VecTN<T,N> operator-(VecTN<T,N> b);
        
        
        /** \brief Scalar multiplication
        **/
        friend VecTN<T,N> operator*(T a, VecTN<T,N> V)
        {
            VecTN<T,N> res;
            for (int i = 0; i < N; i++)
            {
                res[i] = a*V[i];
            }
            return res;
        }

        /** \brief Scalar multiplication
        **/
        friend VecTN<T,N> operator*(VecTN<T,N> V, T a)
        {
            VecTN<T,N> res;
            for (int i = 0; i < N; i++)
            {
                res[i] = a*V[i];
            }
            return res;
        }
            
        /** \brief Scalar division
        **/
        friend VecTN<T,N> operator/(T a, VecTN<T,N> V)
        {
            VecTN<T,N> res;
            for (int i = 0; i < N; i++)
            {
                res[i] = V[i]/a;
            }
            return res;
        }

        /** \brief Scalar division
        **/
        friend VecTN<T,N> operator/(VecTN<T,N> V, T a)
        {
            VecTN<T,N> res;
            for (int i = 0; i < N; i++)
            {
                res[i] = V[i]/a;
            }
            return res;
        }
            
        
    };

template <class T, int N> 
void VecTN<T,N>::init_mem()
{
    memset(&values[0], 0, N*sizeof(T));
}

template <class T, int N> 
VecTN<T,N>::VecTN()
{
    init_mem();
}


template <class T, int N> 
VecTN<T,N>::VecTN(int nToSet, ... )
{
    init_mem();

    va_list arguments;

    va_start(arguments,nToSet);
    for (int i = 0; i < nToSet; i++)
    {
        values[i] = va_arg(arguments, T);
    }
    va_end(arguments);

}

/** \brief Returns value at index
    
    Overloaded [] operator returns value at index i
    @param i Index of Vecor
**/
template <class T, int N> 
T& VecTN<T,N>::operator[](const int i)
{
    if(i<N)
    {
        return values[i];
    }
    else
    {
        std::cerr<<"[VecTN->operator[]] Accessor value of "<<i<<" is greater than max value "<<N-1<<std::endl;
        return values[0];
        // exit(-1);
    }
}

template <class T, int N> 
T& VecTN<T,N>::get(const int i)
{
    if(i<N)
    {
        return values[i];
    }
    else
    {
        std::cerr<<"[VecTN->operator[]] Accessor value of "<<i<<" is greater than max value "<<N-1<<std::endl;
        return values[0];
        // exit(-1);
    }
}

template <class T, int N> 
void VecTN<T,N>::copyValue(const int i, T& retVal) const
{
    if(i<N)
    {
        retVal = values[i];
        retVal = T(retVal);
    }
    else
    {
        std::cerr<<"[VecTN->operator[]] Accessor value of "<<i<<" is greater than max value "<<N-1<<std::endl;
        retVal = values[0];
        retVal = T(retVal);
        // exit(-1);
    }
}



/** \brief Prints Vecor values
    
    Values are written tab delimited with an endl.
    @param out ostream to write to
**/
template <class T, int N> 
void VecTN<T,N>::print(std::ostream& out)
{
    if(this==NULL)
    {
        std::cerr<<"[VecTN<T,N> | print] The vector points to NULL!"<<std::endl;
        return;
    }
    for (int i = 0; i < N; i++)
    {
        out<<values[i]<<"\t";
    }
    out<<std::endl;
}

/** \brief Dot product
**/
template <class T, int N> 
VecTN<T,N> VecTN<T,N>::operator*(VecTN<T,N> b)
{
    VecTN<T,N> res;
    for (int i = 0; i < N; i++)
    {
        res[i] = values[i]*b[i];
    }
    return res;
}

/** \brief Normalizes the Vecor
**/
template <class T, int N> 
VecTN<T,N>& VecTN<T,N>::normalize()
{
    T norm = this->norm();
    
    for (int i = 0; i < N; i++)
    {
        this->values[i] /= norm;
    }
    
    return *this;
}

template <class T, int N> 
double VecTN<T,N>::norm()
{
    T norm = 0;
    
    for (int i = 0; i < N; i++)
    {
        norm += this->values[i]*this->values[i];
    }
    
    norm = sqrt(norm);
    return norm;
}

/** \brief Adds 2 Vecors
**/
template <class T, int N> 
VecTN<T,N> VecTN<T,N>::operator+(VecTN<T,N> b)
{
    VecTN<T,N> res;
    for (int i = 0; i < N; i++)
    {
        res[i] = values[i]+b[i];
    }
    return res;
}

template <class T, int N> 
VecTN<T,N> VecTN<T,N>::operator+(VecTN<T,N>* b)
{
    VecTN<T,N> res;
    for (int i = 0; i < N; i++)
    {
        res[i] = this->values[i]+(*b)[i];
    }
    return res;
}


/** \brief Subtracts 2 Vecors
**/
template <class T, int N>
VecTN<T,N> VecTN<T,N>::operator-(VecTN<T,N> b)
{
    VecTN<T,N> res;
    for (int i = 0; i < N; i++)
    {
        res[i] = values[i]-b[i];
    }
    return res;
}



#endif
