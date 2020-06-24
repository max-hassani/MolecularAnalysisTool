#include "particle.h"

particle::particle(int n, double l1, double l2, double l3):
ix(0), iy(0),iz(0),pos(0.,0.,0.),posU(0.,0.,0.),vel(0.,0.,0.),lx(0.),ly(0.),lz(0.)
{
  ind = n
}
particle::particle():
ix(0), iy(0),iz(0),pos(0.,0.,0.),posU(0.,0.,0.),vel(0.,0.,0.),ind(0)
{
  //particle is initialized without any index
}
particle::~particle(){

}
void apply_pbc(){

}
void unwrapped(){

}
