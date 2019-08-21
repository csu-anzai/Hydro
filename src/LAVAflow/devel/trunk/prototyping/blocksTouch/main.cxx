// Final proposal: combine adjacent rectangles, 
// if they are 'flush': almost touching

#include <iostream>
using namespace std;

struct R
{
    int x1,y1,x2,y2;
    int height() const { return y2-y1; }
    int width() const  { return y2-y1; }

    void normalize() 
    { 
        if (x1>x2) std::swap(x1,x2);
        if (y1>y2) std::swap(y1,y2);
    }

    /*
     * adjacent: return whether two rectangles
     * are adjacent; the tolerance in pixels
     * allow you to specify the gap:
     *    tolerance = 0: require at least one pixel overlap
     *    tolerance = 1: accepts 'flush' adjacent neighbours
     * Negative tolerance require more overlap;
     * tolerance > 1 allows gaps between rects;
     */
    bool adjacent(R const& other, int tolerance=1) const
    {
        return !( (other.x1 - x2) > tolerance
               || (x1 - other.x2) > tolerance
               || (other.y1 - y2) > tolerance
               || (y1 - other.y2) > tolerance);
    }
};


std::ostream& operator<<(std::ostream &os, R const& r)
{
    return os << '(' << r.x1 << "," << r.y1 << ")-(" << r.x2 << ',' << r.y2 << ')';
}

int main()
{
    const int tolerance = 1;
    std::cout << "sample from original question" << std::endl;
    R a = { 0, 0, 3, 1 }; /* a.normalize(); */
    R b = { 0, 2, 2, 3}; /* b.normalize(); */

    std::cout << "a: " << a << "\t b: " << b << std::endl;
    std::cout << boolalpha << "are adjacent? "<< a.adjacent(b,tolerance)<<endl;
    
}