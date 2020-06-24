#include"chain.h"

chain::chain(int nr1):nrAdded(0), ind(NULL), nrBeads(0){
  /*if(nr1 == 0){
    cout << "fatal error: the number of beads are not defined apriory" << endl;
    exit(1000);
  }*/

  //nrBeads = nr1;
  //ind = (int*) malloc(nrBeads*sizeof(int));
  cen.set_to(0.,0.,0.);
  rG.set_to(0.,0.,0.);
  for(int i = 0; i < nrBeads; i++){
    ind[i] = -1;
    //con[i]=-1;
  }
  //nrAdded = 0;
  indSt= -1;
  indEnd=-1;
}
chain::chain():nrBeads(1),cen(0.,0.,0.),rG(0.,0.,0.),ind(NULL),indSt(0),indEnd(0),nrAdded(0)
{
  //when there is no argument is passed only the default values are set
}
chain::~chain(){
  free(ind);
  //free(con);
}
void chain::clear(){
  for(int i = 0; i < nrBeads ; i++){
    ind[i]= -1;
    //con[i] = -1;
  }
}
Vec chain::calCenMass(Vec *pos, double *masses, int *typen){
  int ii;
  double denom = 0.;
  for(int i = 0; i < nrBeads; i++){
    ii = ind[i];
    cen += pos[ii]*masses[typen[ii]];
    denom += masses[typen[ii]];
  }
  if(nrBeads != 0){
      cen = cen * (1./denom);
  }
  return cen;
}

void chain::calRg(Vec *pos, double *masses, int *typen){
  calCenMass(pos, masses, typen);
  for(int i = 0; i < nrBeads; i++){
    rG.x += pow((pos[ind[i]].x - cen.x),2.);
    rG.y += pow((pos[ind[i]].y - cen.y),2.);
    rG.z += pow((pos[ind[i]].z - cen.z),2.);
  }
  rG = rG *(1./(double)nrBeads);
}
void chain::allocInd(int n){
  ind = (int*) malloc(n*sizeof(int));
  nrBeads = n;
}
void chain::addOne(int par_id){
  ind[nrAdded]=par_id;
  nrAdded++;
}
Vec chain::calEndToEnd(Vec *pos){
  Vec R;
  R.reset();
  for(int i = 1 ; i < nrBeads; i++){
    R += pos[ind[i]] - pos[ind[i-1]] ;
  }
  return R;
}
