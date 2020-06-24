#include<iostream>
#include<cmath>
#include"p_vec_ten.hh"

using namespace std;

typedef p_vec<> Vec;

class chain{
public:
	//int num;
	int nrAdded;
	int *ind;//indices of particles making up the chain

	//[MaxChainBeads]
	int nrBeads;
  //int chId; // chain ID
	Vec cen; // position of the center of mass
	Vec rG; // radius of gyration
  int indSt, indEnd;// starting and ending indices.
	chain();
	chain(int nr1); // defaulf constructor
	~chain(); // defaulf deconstructor
	void clear(); // clear the data
	Vec calCenMass(Vec *pos, double *masses, int *typen); // calculate the center of mass of the chain
	void calRg(Vec *pos, double *masses, int *typen); // calculate the radius of gyration for a chain
	void allocInd(int n);
	void addOne(int par_id);
	Vec calEndToEnd(Vec *pos);
};
