#ifndef LAVA_UTIL_SORTOPS_H
#define LAVA_UTIL_SORTOPS_H

#include <vector>
#include <algorithm>




namespace Util
{

	template<class T, class S>
	bool pairAscendingComp( const std::pair<T, S>& i, const std::pair<T, S>& j ) {
	    if( i.first < j.first ) return true;
	    if( j.first < i.first ) return false;
	    return i.second < j.second;
	}

	template<class T, class S>
	bool pairDescendingComp( const std::pair<T, S>& i, const std::pair<T, S>& j ) {
	    if( i.first < j.first ) return false;
	    if( j.first < i.first ) return true;
	    return j.second < i.second;
	}

	template<class T, class S>
	void sortVectorOfPairs(std::vector<std::pair<T,S> >& toBeSorted, bool doAscending = true)
	{
		if(doAscending)
		{
			std::stable_sort(toBeSorted.begin(), toBeSorted.end(), pairAscendingComp<T,S>);
		}
		else
		{

			std::stable_sort(toBeSorted.begin(), toBeSorted.end(), pairDescendingComp<T,S>);
		}

	}




};

#endif