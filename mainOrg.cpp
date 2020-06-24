#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <time.h>
#include "p_vec_ten.hh"

#define SQ(x)     ((x) * (x))
#define CORREL_STR
void * memset(void *dest, int c, size_t count);
typedef p_vec<> Vec;
using namespace std;

int NRPART;
int NRBONDS;
int NRChains;

double LX;
double LY;
double LZ;
double LX_0;
double LY_0;
double LZ_0;
double Z_BULK;

int t_first;
int num_times;
int t_step;
int t_delta;


double intvl;
double gdot;

double LXINV;
double LYINV;
double LZINV;
double LXINV_0;
double LYINV_0;
double LZINV_0;
double XMIN;
double YMIN;
double ZMIN;
double XMAX;
double YMAX;
double ZMAX;
double XMIN_0;
double YMIN_0;
double ZMIN_0;
double XMAX_0;
double YMAX_0;
double ZMAX_0;

const int MaxChainBeads = 300;


class chain{
public:
	//int num;
	int nr;
	int ind[MaxChainBeads],con[MaxChainBeads];
  int typ;
  int indSt, indEnd;
	chain(){
    //ind[MaxChainBeads]; //= (int *) malloc(MaxChainBeads * sizeof(int));
    //con[MaxChainBeads]; //(int *) malloc(MaxChainBeads * sizeof(int));
    for(int i = 0; i < MaxChainBeads; i++){
      ind[i] = -1;
      con[i]=-1;
    }
    nr = 0;
    indSt= -1; indEnd=-1;
	}
  ~chain(){
    free(ind);
    free(con);
  }
  int calSize(){
    int j = 0;
    nr = 0;
    while(ind[j] >= 0){
      nr++;
    }
    return(nr);
  }

	void clear(){
		for(int i = 0; i < MaxChainBeads ; i++){
      ind[i]= -1;
      con[i] = -1;
    }
	}
};

void loadData0(char *fname, Vec* data, Vec* vel, int *typen, int* chID, int **bonds, int *bondTyp);
void SetNrBondsPerAtom(int **bonds, int **parBondIDs, int *parBondsNr, int **BondedTo);
//void FormChains(chain *chList, int* chID, int **bonds, int *bondsTyp,int **BondedTo, int *parBondsNr);
void FormChains(chain *chList, int* chID);
void MaxOfChID(int *val, int N);
void CalChainRG(Vec *pos, chain *chList, Vec *chCM, Vec* chRG);
void loadData(int t, Vec* data, int *typen, Vec* posU, Vec* vel);
void CalChainCM(Vec *pos, chain *chList, Vec *chCM);
Vec timeAvg(Vec* val, int N);
void apply_pbc(int time, Vec& mol, Vec& out);

int main(int argc, char **argv){
	clock_t start, end;
	start = clock();
	if(argc != 9){
		cout << "not sufficient input; please refer to the main function " << endl;
		return 1;
	}
	LX = atof(argv[1]);
	NRPART = atoi(argv[2]);
	NRBONDS = atoi(argv[3]);
	NRChains = atoi(argv[4]);
	num_times = atoi(argv[5]);
	t_first = atoi(argv[6]);
	intvl = atof(argv[7]);
	t_step = atoi(argv[8]);
	LZ=LX;LY=LX;
	LXINV = 1./LX;
	LYINV = 1./LY;
	LZINV = 1./LZ;
	XMIN = -LX/2.;
	YMIN = -LY/2.;
	ZMIN = -LZ/2.;
	XMAX = LX/2.;
	YMAX = LY/2.;
	ZMAX = LZ/2.;

	cout << "MIN and MAXs:"<< endl;
	cout << XMIN << "\t" << XMAX << endl;
	cout << YMIN << "\t" << YMAX << endl;
	cout << ZMIN << "\t" << ZMAX << endl;

	cout << "allocating the arrayes" << endl;
	char fname[30];
	Vec *pos,*vel,*posU,*chRG,*chCM;
	int *chID, **bond, *bondTyp, *typen;
	pos = (Vec*) malloc(NRPART*sizeof(Vec));
	posU = (Vec*) malloc(NRPART*sizeof(Vec));
	vel = (Vec*) malloc(NRPART*sizeof(Vec));

	chCM = (Vec*) malloc(NRChains*sizeof(Vec));
	chRG = (Vec*) malloc(NRChains*sizeof(Vec));

	chID = (int*) malloc(NRPART*sizeof(int));
	typen = (int*) malloc(NRPART*sizeof(int));
	bond = (int**) malloc(NRBONDS*sizeof(int*));
	bondTyp = (int*) malloc(NRBONDS*sizeof(int));
	for(int i = 0; i < NRBONDS; i++){
		bond[i] = (int*) malloc(2*sizeof(int));
	}
	cout << "loading the initial data" << endl;
	sprintf(fname, "/home/users/hassasp9/lammps_runs/Polymer/2ndTry/N_300__M_1600.init");
	loadData0(fname, pos, vel, typen, chID, bond,bondTyp);
	cout << "initial data loaded" << endl;
	MaxOfChID(chID, NRPART);

	chain *chList;
	chList = (chain*) malloc(NRChains*sizeof(chain));

 	for(int i = 0; i < NRChains; i++){
		chList[i].nr = 0;
		chList[i].indSt= -1; chList[i].indEnd=-1;
	}

	Vec *posRep, *posCh;
	posRep = (Vec*) malloc(NRPART*sizeof(Vec));
	posCh = (Vec*) malloc(300*sizeof(Vec));
 	FormChains(chList, chID);
        cout << "the chains are formed" << endl;
        ofstream oFileN;
        for(int t = 0; t < num_times; t++){
                int tCurrent = t_first + (t)*t_step;
                cout << "loading data at: " << tCurrent << endl;
                loadData(tCurrent,pos, typen, posU,vel);
								for(int ip = 0; ip < NRPART; ip++){
									apply_pbc(tCurrent, pos[ip],posRep[ip]);
									posRep[ip].x = posRep[ip].x + 1.5 * LX_0;
								}
                //cout << "loading data at: " << tCurrent << endl;
                for(int i = 0 ; i < NRChains; i++){
                  if(i%500 == 0){
                          sprintf(fname,"dump_ch%d.%d",i,tCurrent);
                          oFileN.open(fname);
													oFileN << "ITEM: TIMESTEP" << endl;
													oFileN << tCurrent << endl;
													oFileN << "ITEM: NUMBER OF ATOMS" << endl;
													oFileN << NRPART+300 << endl;
													oFileN << "ITEM: BOX BOUNDS pp pp pp" << endl;
													oFileN << XMIN_0 << " " << XMAX_0<< endl;
													oFileN << YMIN_0<< " " << YMAX_0<< endl;
													oFileN << ZMIN_0<< " " << ZMAX_0<< endl;
													oFileN << "ITEM: ATOMS id type xu yu zu" << endl;
													for(int ip = 0; ip < NRPART; ip++){
                                  oFileN << ip+1 << " " << typen[ip] << " " << posRep[ip] << endl;
                          }
													// Hassani: here we add the extra chain
                          for(int j = 0; j < MaxChainBeads; j++){
																	apply_pbc(tCurrent, pos[chList[i].ind[j]],posCh[j]);
                                  oFileN << j+1+NRPART << " " << typen[chList[i].ind[j]] << " " << posCh[j] << endl;
                          }
                          oFileN.close();
                  }
                }
        }
	/*ofstream CM_file("CM-t.dat");
	ofstream RG_file("RG-t.dat");
	Vec Rg_t , Cm_t;
	ofstream testChID("chID-test.dat");
	for(int i = 0; i < NRPART; i++){
		testChID << i+1 << " " << chID[i]+1 << endl;
	}
	testChID.close();
	for(int t = 0; t < num_times; t++){
		double tCurrent = t_first + (t)*t_step;
		loadData(tCurrent,pos, typen, posU,vel);
		CalChainCM(pos, chList, chCM);
		CalChainRG(pos, chList, chCM,chRG);
		Rg_t = timeAvg(chRG, NRChains);
		Cm_t = timeAvg(chCM, NRChains);
		RG_file << tCurrent << "\t" << Rg_t << endl;
		CM_file << tCurrent << "\t" << Cm_t << endl;
	}
	CM_file.close();
	RG_file.close();*/
	free(posRep);free(posCh);
	end = clock();
	cout << "Time required for execution:"<< (double)(end-start)/CLOCKS_PER_SEC<< " seconds." << endl;
	cout << "DONE." << endl;
	return 0;
}
void apply_pbc(int time, Vec& mol, Vec& out)
{
	out = mol;
	out.x = mol.x - ((int) floor((mol.x - XMIN_0) / LX_0)) * LX_0;
	out.z = mol.z - ((int) floor((mol.z - ZMIN_0) / LZ_0))* LZ_0;
	out.y = mol.y - ((int) floor((mol.y - YMIN_0) / LY_0) ) * LY_0;
	if ( out.x < XMIN_0 ) out.x += LX_0 ;
  	if ( out.y < YMIN_0 ) out.y += LY_0;
	if ( out.z < ZMIN_0 ) out.z += LZ_0 ;
	//return(out);
}
void loadData(int t, Vec* data, int *typen, Vec* posU, Vec* vel)
{
	int i, iniLine= 5, iniLine1 = 3,tempInd;
	char fname[200];
	// Hassani + Fathollah:
//	sprintf(fname, "../dump-p.%d",t);
	sprintf(fname, "/home/users/varnifwh/SMP-Strain200-full-cycle/dump.%d",t);
	ifstream inputFile(fname);
	//ofstream testInp("testInput.dat");
	if (!inputFile)
	{
		cout << "ERROR: no input file: " << fname << endl;
		exit(1);
	}
	string temp;
	double tempD;
	for(i = 0; i < iniLine; i++){
		getline(inputFile,temp);
		//testInp << temp << endl;
	}
	inputFile >> XMIN_0 >> XMAX_0;
	//testInp << XMIN_0 << "\t" << XMAX_0<< endl;
	inputFile >> YMIN_0 >> YMAX_0;
	//testInp << YMIN_0 << "\t" << YMAX_0<< endl;
	inputFile >> ZMIN_0 >> ZMAX_0;
	//testInp << ZMIN_0 << "\t" << ZMAX_0<< endl;


	LX_0 = XMAX_0 - XMIN_0;
	LY_0 = YMAX_0 - YMIN_0;
	LZ_0 = ZMAX_0 - ZMIN_0;
	LXINV_0 = 1. / LX_0;
	LYINV_0 = 1. / LY_0;
	LZINV_0 = 1. / LZ_0;
	cout << "LiS:\n" << LX_0 << "\n" << LY_0 << "\n" << LZ_0 << endl;
	getline(inputFile,temp);
	//testInp << temp << endl;
	getline(inputFile,temp);
	//testInp << temp << endl;
	//getline(inputFile,temp);
	//testInp << temp << endl;
	inputFile.precision(5);
	for(i = 0; i < NRPART; i++){
		inputFile >> tempInd;
		inputFile >> typen[tempInd - 1] >> data[tempInd - 1].x >> data[tempInd - 1].y >> data[tempInd - 1].z;
		posU[tempInd - 1] = data[tempInd - 1];
		//testInp << tempInd << "\t" << typen[tempInd - 1] << "\t" << data[tempInd - 1] << "\t" << vel[tempInd - 1] << endl;
		//apply_pbc(t,data[tempInd - 1]);
	}
	//testInp.close();
	inputFile.close();
	/*for(int i = 0; i < iniLine1; i++){
		getline(inputFile,temp);
	}
	for(int i = 0; i < NRBONDS; i++){
		inputFile >> tempInd;
		inputFile >> bondTyp[tempInd-1] >> bonds[tempInd-1][0] >> bonds[tempInd-1][1];
	}*/
}
void loadData0(char *fname, Vec* data, Vec* vel, int *typen, int* chID, int **bonds, int *bondTyp)
{
	cout << "in the load function" << endl;
	int i, iniLine= 31, iniLine1 = 4,tempInd;

	//sprintf(fname, "d.%d",t);
	//ofstream test("testInp.dat");
	ifstream inputFile(fname);
	if (!inputFile)
	{
		cout << "ERROR: no input file: " << fname << endl;
		exit(1);
	}
	string temp;
	double tempD;
	cout << "loading the header" << endl;
	for(i = 0; i < iniLine; i++){
		getline(inputFile,temp);
		//test << temp << endl;
	}
  int nochTmp, boxInd;
	inputFile.precision(5);
	cout << "loading the particles" << endl;
	for(i = 0; i < NRPART; i++){
		inputFile >> tempInd;
		inputFile >> nochTmp >> typen[tempInd - 1] >> data[tempInd - 1].x >> data[tempInd - 1].y >> data[tempInd - 1].z >> boxInd >> boxInd >> boxInd;
    chID[tempInd - 1] = nochTmp - 1;
    if(chID[tempInd - 1] >= NRChains){
      cout << "chain id is out of limit" << endl;
      exit(13);
    }
		//test << tempInd << " "<< chID[tempInd - 1] << " " << typen[tempInd - 1] << " " << data[tempInd - 1].x <<" " << data[tempInd - 1].y << " "<< data[tempInd - 1].z << endl;
		//posU[tempInd - 1] = data[tempInd - 1];
		//apply_pbc(t,data[tempInd - 1]);
	}
	/*for(int i = 0; i < iniLine1; i++){
		getline(inputFile,temp);
		//test << temp << endl;
	}
	cout << "loading the velocities" << endl;
  for(int i = 0; i < NRPART; i++){
    inputFile >> tempInd;
		inputFile >> vel[tempInd - 1].x >> vel[tempInd - 1].y >> vel[tempInd - 1].z;
		//test << tempInd << " " << vel[tempInd - 1] << endl;
  }
  for(int i = 0; i < iniLine1; i++){
    getline(inputFile,temp);
		//test << temp << endl;
  }
	cout << "loading the bonds" << endl;
	for(int i = 0; i < NRBONDS; i++){
		inputFile >> tempInd;
		if(tempInd < 1 || tempInd > NRBONDS){
			cout << "the number of bonds does not match the expectations" << endl;
			cout << tempInd << endl;
			exit(15);
 		}
		inputFile >> bondTyp[tempInd-1] >> bonds[tempInd-1][0] >> bonds[tempInd-1][1];
		//test << tempInd <<" " << bondTyp[tempInd-1] << " "<< bonds[tempInd-1][0] << " "<< bonds[tempInd-1][1] << endl;
	}*/
  inputFile.close();
	//test.close();
}
void SetNrBondsPerAtom(int **bonds, int **parBondIDs, int *parBondsNr, int **BondedTo){

  for(int i = 0; i < NRPART; i++){
    parBondsNr[i]=0;
    parBondIDs[i][0] = -1; parBondIDs[i][1] = -1;
    BondedTo[i][0] = -1;BondedTo[i][1] = -1;
  }
  for(int i = 0 ; i < NRBONDS; i++){
    for(int j = 0; j < 2; j++){
      int ind =bonds[i][j];
      int ind2;
      if(j == 0){
        ind2 = bonds[i][j+1];
      }else{
        ind2 = bonds[i][j-1];
      }
      int nr = parBondsNr[ind];
      if(nr > 1){
        cout << "number of bonds per particle exceeds the expected value" << endl;
        exit(10);
      }
      parBondIDs[ind][nr] = i;
      BondedTo[ind][nr] = ind2;
      parBondsNr[ind] ++;
    }
  }
}
void FormChains(chain *chList, int* chID){
  int bondPar1,bondPar2,bondPar3;
  for(int i = 0; i < NRPART; i++){
		int n = chList[chID[i]].nr;
    chList[chID[i]].ind[n] = i;
    chList[chID[i]].nr ++;
  }
  /*for(int i = 0; i < NRChains; i++){
    for(int j = 0; j < chList[i].nr; j++){
      int ind = chList[i].ind[j];
      bondPar1 = BondedTo[ind][0];
      bondPar2 = BondedTo[ind][1];
      if(bondPar1 != -1 && bondPar2 != -1){
        if(chID[bondPar1] != chID[bondPar2]){
          if(chList[chID[ind]].indSt == -1){
            chList[chID[ind]].indSt = ind;
          }else if(chList[chID[ind]].indEnd == -1){
            chList[chID[ind]].indEnd = ind;
          }else{
            cout << "the end/start particle could not be set " << endl;
            exit(12);
          }
        }
      }else if(bondPar1 != -1 || bondPar2 != -1){
        if(chList[chID[ind]].indSt == -1){
          chList[chID[ind]].indSt = ind;
        }else if(chList[chID[ind]].indEnd == -1){
          chList[chID[ind]].indEnd = ind;
        }else{
          cout << "the end/start particle could not be set " << endl;
          exit(12);
        }
      }else{
        cout << "bonding particles could not be found! strange !!!!" << endl;
        exit(13);
      }
    }
  }*/
}
void MaxOfChID(int *val, int N){
  int max=val[0];
  for(int i = 1; i < N; i++){
    if(val[i] > max){
      max = val[i];
    }
  }
  cout << "max of chain ids are: "<< max << endl;
  if(max >= NRChains){
    cout << "the maximum chain id is larger than expected" << endl;
    exit(11);
  }
}
void CalChainCM(Vec *pos, chain *chList, Vec *chCM){
	//ofstream test("testCH-CM.dat");
	for(int i = 0; i < NRChains; i++){
		chCM[i].set_to(0.,0.,0.);
		int N = chList[i].nr;
		if(i == 0){
			for(int j = 0; j < N; j++){
				int ind = chList[i].ind[j];
				if(ind < 0 || ind > NRPART){
					cout << "ind error" << endl;
					exit(15);
				}
				chCM[i] = chCM[i] + pos[ind];
			}
			chCM[i] = chCM[i] * (1./(double)N);
			//test << i << "\t"<< chCM[i] << endl;
		}else{
		//cout << "Ns: " << N << endl;
			for(int j = 0; j < N; j++){
				int ind = chList[i].ind[j];
				if(ind < 0 || ind > NRPART){
					cout << "ind error" << endl;
					exit(15);
				}
				chCM[i] = chCM[i] + pos[ind];
			}
			chCM[i] = chCM[i] * (1./(double)N);
			//test << i << "\t"<< chCM[i] << endl;
		}
	}
	//test.close();
}
void CalChainRG(Vec *pos, chain *chList, Vec *chCM, Vec* chRG){
	//ofstream test("testCH.dat");
	for(int i = 0; i < NRChains; i++){
		int N = chList[i].nr;
		if(i == 1000){
			chRG[i].reset();
		//	test << "CM: " << chCM[i] << endl;
			for(int j = 0; j < N; j++){
				int ind = chList[i].ind[j];
			//	test << ind << "\t" << pos[ind] << endl;
				if(ind < 0 || ind > NRPART){
					cout << "ind error" << endl;
					exit(15);
				}
				chRG[i].x += (pos[ind].x - chCM[i].x)*(pos[ind].x - chCM[i].x);
				chRG[i].y += (pos[ind].y - chCM[i].y)*(pos[ind].y - chCM[i].y);
				chRG[i].z += (pos[ind].z - chCM[i].z)*(pos[ind].z - chCM[i].z);
			}
			chRG[i].x = sqrt(chRG[i].x * (1./(double)N));
			chRG[i].y = sqrt(chRG[i].y * (1./(double)N));
			chRG[i].z = sqrt(chRG[i].z * (1./(double)N));
			//test << "endlich "<< chRG[i] << endl;
		}else{
			chRG[i].reset();
			//int N = chList[i].nr;
			for(int j = 0; j < N; j++){
				int ind = chList[i].ind[j];
				if(ind < 0 || ind > NRPART){
					cout << "ind error" << endl;
					exit(15);
				}
				chRG[i].x += (pos[ind].x - chCM[i].x)*(pos[ind].x - chCM[i].x);
				chRG[i].y += (pos[ind].y - chCM[i].y)*(pos[ind].y - chCM[i].y);
				chRG[i].z += (pos[ind].z - chCM[i].z)*(pos[ind].z - chCM[i].z);
			}
			chRG[i].x = sqrt(chRG[i].x * (1./(double)N));
			chRG[i].y = sqrt(chRG[i].y * (1./(double)N));
			chRG[i].z = sqrt(chRG[i].z * (1./(double)N));
		}
	}
	//test.close();
}
Vec timeAvg(Vec *val, int N){
	Vec out;
	out.reset();
	cout << "N: " << N << endl;
	for(int i = 0; i < N; i++){
		out = out + val[i];
	}
	out = out * (1./(double)N);
	cout << "out " << out << endl;
	return(out);
}
