#include "ArrayOperations.h"
#include <iostream>

using namespace std; 

std::vector<const char*> VecStringToVecChar(const std::vector<std::string> &vars)
{
    std::vector<const char*> varsChar;
    std::transform(vars.begin(), vars.end(), back_inserter(varsChar), convert);
    
    return varsChar;
}

const char *convert(const std::string & s)
{
   return s.c_str(); 
}

std::vector<std::string> CharToVecOfString(const char** charArray, int numStrings)
{
    std::vector<std::string> vecOfStrings;
    for (int i=0; i<numStrings; i++)
    {
        if (charArray[i] != NULL)
            vecOfStrings.push_back(charArray[i]);
    }

    return vecOfStrings;
}
