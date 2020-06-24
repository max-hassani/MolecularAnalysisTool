#include "p_vec_ten.hh"
typedef p_vec<> Vec;
class particle{
public:
 Vec pos, posU, vel;
 int ix,iy, iz;
 int ind;
 double lx, ly, lz;
 particle(int n);
 particle();
 ~particle();
 void apply_pbc();
 void unwrapped();
}
