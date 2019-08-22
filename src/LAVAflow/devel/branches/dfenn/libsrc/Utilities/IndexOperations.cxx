#include "IndexOperations.h"

namespace Util
{

int sub2ind(int i, int j, int k, int nI, int nJ)
{

	return i + nI*(j + nJ*k);
}


}; // end namespace