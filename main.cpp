#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <time.h>
#include "p_vec_ten.hh"
#include "chain.h"
#include <string.h>
//#include "particle.h"
#define SQ(x)     ((x) * (x))
#define CORREL_STR
void * memset(void *dest, int c, size_t count);
typedef p_vec<> Vec;
using namespace std;

const int NLIST = 10;
const int NLIST_digit = 5;
string listInput[] = {"num_times","intrvl","NRPART","NRBONDS","NR_ATOM_TYPE", "init_data", "time_series_file","dump_path", "atom_style", "periodicity"};
//"LX","NRPART","NRBONDS","NRChains","t_first","t_step"
string bond_path, init_data, time_series_file, dump_path, atom_style;
bool dynamic_box;
int flag_num_times(0),flag_bondPath(0), flag_initData(0), flag_time_series(0), flag_dump(0), flag_NRPART(0), flag_Type(0), flag_mass(0), flag_nr_bond(0), flag_nr_chain(0),flag_intrvl(0), flag_cubic(0), flag_atom_style(0), flag_NR_ATOM_TYPE(0);
int flag_per[3];
int NRPART;
int NRBONDS;
int NRChains(1);
int nrBead_perChain;

int NR_ATOM_TYPE;
int NR_BOND_TYPE;
int NR_ANGLE;
int NR_ANGLE_TYPE;
int NR_DIHEDRAL;
int NR_DIHEDRAL_TYPE;
int NR_IMPROPER;
int NR_IMPROPER_TYPE;
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


double intrvl(0.005);//interval
double gdot;
int nr_hist_bins = 100;// default value for the number of bins

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


void compare_assign(int caseNr, string val);
void compare_assign(int caseNr, double val);
void readInput(char *fname);
void finalControlInput();
void alloc_Vec(Vec *&v, int n);
void alloc_Vec(Vec **&v, int n);
void alloc_double(double *&x, int n);
void alloc_int(int *&ix, int n);
void alloc_int(int **&ix, int n);
void alloc_chain(chain *& ch, int n);
void testInput();
void loadData0( Vec *pos, Vec *posU,  int *typen, int **bonds, int *chID, double *masses);
int analys_chID(int *chID);
void makeUpTheChains(chain *& chains, int *chID);
void loadTimeSeries(int *t_series);
void loadDump(int t, Vec* pos, int *typen, Vec* posU);
Vec calculate_msd(chain *chains,Vec *posU, Vec *cen0, double *masses, int *typen);
void store_init_cen_masses(chain *chains, Vec *cen0, Vec *posU, double *masses, int *typen);
Vec calculate_msd(chain *chains,Vec *posU, Vec *cen0, double *masses, int *typen);
Vec calculate_Rg(chain *chains, Vec *posU, double *masses, int *typen, Vec *&rG_chain);
void calEndToEnd(chain *chains, Vec *posU, Vec *rEE_chain);
void calEndToEndCorr(Vec **rEndEnd_t, int *t_series);
void store_time_arrays( string fname0, Vec *arr, int *t_series );

int main(int argc, char **argv){
	clock_t start, end;
	start = clock();

  char path[1000];

	strcpy(path, argv[1]);
	readInput(path);
	finalControlInput();
	/*if(flag_LX){
		XMIN = -LX/2.;
		XMAX = LX/2.;
		LXINV = 1./LX;
		XMIN_0 = XMIN;
		XMAX_0 = XMAX;
		LXINV_0 = LXINV;
	}
	if(flag_LY){
		YMIN = -LY/2.;
		YMAX = LY/2.;
		LYINV = 1./LY;
		YMIN_0 = YMIN;
		YMAX_0 = YMAX;
		LYINV_0 = LYINV;
	}
	if(flag_LZ){
		ZMIN = -LZ/2.;
		ZMAX = LZ/2.;
		LZINV = 1./LZ;
		ZMIN_0 = ZMIN;
		ZMAX_0 = ZMAX;
		LZINV_0 = LZINV;
	}*/

	Vec *pos_ref(NULL), *posU_ref(NULL), *pos0(NULL), *pos1(NULL), *posU0(NULL), *posU1(NULL), *cen0(NULL), *msd_avg(NULL), **radGyr_t(NULL), **rEndEnd_t(NULL), *radGyr_avg(NULL);
	int *typen(NULL), *t_series(NULL);
	chain *chains(NULL);
	int **bonds(NULL);
	double *masses(NULL);
	cout << "allocating memory" << endl;
	alloc_Vec(pos_ref, NRPART);
	alloc_Vec(posU_ref, NRPART);
	alloc_Vec(pos0, NRPART);
	alloc_Vec(pos1, NRPART);
	alloc_Vec(posU0, NRPART);
	alloc_Vec(posU1, NRPART);
	alloc_Vec(cen0, NRPART);
	alloc_Vec(msd_avg , num_times);
	alloc_Vec(radGyr_t , num_times);
	alloc_Vec(rEndEnd_t , num_times);
	alloc_Vec(radGyr_avg , num_times);
	//alloc_Vec(rEndEnd_avg, num_times);

	alloc_int(typen, NRPART);
	alloc_int(bonds, NRBONDS);
	alloc_int(t_series, num_times);

	alloc_double(masses,NR_ATOM_TYPE);

	for(int i = 0; i < NRBONDS; i++){
		alloc_int(bonds[i], 3);
	}

	int *chID;
	alloc_int(chID, NRPART);
	loadTimeSeries(t_series);
	cout << "loading the initial data file" << endl;
	loadData0( pos_ref, posU_ref, typen, bonds, chID, masses);
	/*
	if(flag_initData){

		if(!dynamic_box){
			cout << "the box is not dynamic: the right setting is being applied!" << endl;
			LX = LX_0;
			XMIN = XMIN_0;
			XMAX = XMAX_0;
			LXINV = LXINV_0;

			LY = LY_0;
			YMIN = YMIN_0;
			YMAX = YMAX_0;
			LYINV = LYINV_0;

			LZ = LZ_0;
			ZMIN = ZMIN_0;
			ZMAX = ZMAX_0;
			LZINV = LZINV_0;
		}
		cout << "calculating the number of available chains" << endl;
		NRChains = analys_chID(chID);
		//cout << "number of chains: " << NRChains << endl;
		if(NRPART % NRChains != 0){
			cout << "Fatal error: the number of particles not devidable to the number of chains" << endl;
			cout << "Here it is assumed that the chains all have the same number of particles inside" << endl;
			exit(5);
		}
		nrBead_perChain = NRPART / NRChains;
		cout << "number of chains: " << NRChains << endl;
		cout << "number of particle per chain: " << nrBead_perChain << endl;
		cout << "allocating the chains!" << endl;
		alloc_chain(chains, NRChains);
		for(int i = 0; i < NRChains; i++){
			chains[i].allocInd(nrBead_perChain);
		}

		for(int t = 0; t < num_times; t++){
			alloc_Vec(radGyr_t[t], NRChains);
			alloc_Vec(rEndEnd_t[t], NRChains);
		}

		cout << "chains are allocated" << endl;
		cout << "assigning particles to the chains!" << endl;
		makeUpTheChains(chains, chID);
		for(int t = 0; t < num_times; t++){
			cout << "time: " << t_series[t] << endl;
			cout << "dump file is being loaded!" << endl;
			loadDump(t_series[t], pos1, typen, posU1);
			if(t == 0){
				pos0 = pos1;
				posU0 = posU1;
				cout << "calculating the initial center of mass!" << endl;
				store_init_cen_masses(chains, cen0, posU0, masses,typen);
				cout << "calculating the mean square displacement!" << endl;
				msd_avg[t] = calculate_msd(chains,posU1, cen0, masses, typen);
				cout << "calculating the radius of gyration!" << endl;
				radGyr_avg[t] = calculate_Rg(chains, posU1, masses, typen, radGyr_t[t]);
				cout << "calculating the end-to-end vector!" << endl;
				calEndToEnd(chains, posU1,rEndEnd_t[t]);
			}else{
				cout << "calculating the initial center of mass!" << endl;
				msd_avg[t] = calculate_msd(chains,posU1, cen0, masses, typen);
				cout << "calculating the radius of gyration!" << endl;
				radGyr_avg[t] = calculate_Rg(chains, posU1, masses, typen,radGyr_t[t]);
				cout << "calculating the end-to-end vector!" << endl;
				calEndToEnd(chains, posU1, rEndEnd_t[t]);
			}
		}
		store_time_arrays( "msd_molecular.dat", msd_avg, t_series );
		store_time_arrays( "Averaged_RadiusGyration_time.dat", radGyr_avg, t_series);
		//calEndToEndCorr(rEndEnd_t, t_series);

	}else if(flag_bondPath){
		cout << "Bond file was chosen as an initial input!" << endl;
		cout << "Unfortunately, this feature is under development."<< endl;
		cout << "Please give as an input an initial data file from lammps." << endl;
		cout << "This is essential for having the bonds between the atoms" << endl;
		cout << "Job aborted!" << endl;
		return 1;
	}
	free(t_series);
	free(pos_ref);
	free(posU_ref);
	free(pos0);
	free(pos1);
	free(posU0);
	free(posU1);
	free(cen0);
	free(msd_avg);
	for(int t = 0; t < NRPART; t++){
		free(radGyr_t[t]);
		free(rEndEnd_t[t]);
	}
	free(radGyr_t);
	free(rEndEnd_t);
	free(radGyr_avg);
	free(typen);
	for(int i = 0; i < NRBONDS; i++){
		free(bonds[i]);
	}
	free(bonds);
	free(t_series);
	free(masses);
	free(chID);
*/
	end = clock();
	cout << "Time required for execution:"<< (double)(end-start)/CLOCKS_PER_SEC<< " seconds." << endl;
	cout << "DONE." << endl;
	return 0;
}
void loadTimeSeries(int *t_series){
	const int size_string = time_series_file.length();
	char fname[size_string+1];
	strcpy(fname, time_series_file.c_str());
	ifstream inputFile(fname);
	if (!inputFile)
	{
		cout << "ERROR: no input file: " << fname << endl;
		exit(1);
	}
	ofstream oFile("testTime.dat");
	for(int t = 0; t < num_times; t++){
		inputFile >> t_series[t];
		oFile << t << " " << t_series[t] << endl;
	}
	inputFile.close();
}
void alloc_Vec(Vec *& v, int n){
	v = (Vec*) malloc(n *sizeof(Vec));
}
void alloc_Vec(Vec **& v, int n){
	v = (Vec**) malloc(n *sizeof(Vec*));
}
void alloc_double(double *& x, int n){
	x = (double*) malloc(n *sizeof(double));
}
void alloc_int(int *& ix, int n){
	ix = (int*) malloc(n *sizeof(int));
}
void alloc_int(int **&ix, int n){
	ix = (int**) malloc(n *sizeof(int*));
}
void alloc_chain(chain *& ch, int n){
	ch = (chain*) malloc(n *  sizeof(chain));
	/*for(int i = 0; i < n; i++){
		ch[i]=chain(nrBead_perChain,i+1);
	}*/
}
void finalControlInput(){
	int error(0);
	if(!flag_NRPART){
		cout << "Fatal error: number of particles is not defined" << endl;
		error = 1;
	}
	if(!flag_num_times){
		cout << "Fatal error: number of times is not given" << endl;
		error = 1;
	}
	if(!flag_time_series){
		cout << "Fatal error: the path to the time series is not given" << endl;
		error = 1;
	}
	if(!flag_dump){
		cout << "Fatal error: the path to the dump files is not given" << endl;
		error = 1;
	}
	if(!flag_nr_bond){
		cout << "Fatal error: number of bonds is not defined" << endl;
		error = 1;
	}
	if(!flag_NR_ATOM_TYPE){
		cout << "fatal error: the number of atom type is not given." << endl;
		//cout << "therfore, the defaulf value of 1 is selected" << endl;
		error = 1;
	}
	if(!flag_atom_style){
		cout << "Fatal error: no atomic style is given." << endl;
		error = 1;
	}
	if(!flag_intrvl){
		cout << "Warning: no interval is defined" << endl;
		cout << "Warning: the default value for the interval is set to 0.005" << endl;
	}
	if(!flag_per[0] || !flag_per[1] || !flag_per[2]){
		cout << "the input for the periodicity is not given correctly" << endl;
		error = 1;
	}
	if(error){
		cout << "A fatal error has occured.\nThe program stops now" << endl;
		exit(3);
	}
}
void compare_assign(int caseNr, string val){
	//int flag_bondPath, flag_initData, flag_time_series, flag_dump;
	switch (caseNr) {
		case NLIST_digit: //NRBONDS = (int) val;
		init_data = val;
		cout << "init_data= " << init_data<< endl;
		flag_initData = 1;
		break;
		case NLIST_digit+1:
		time_series_file = val;
		cout << "time_series_file= " << time_series_file<< endl;
		flag_time_series = 1;
		break;
		case NLIST_digit+2:
		dump_path = val;
		cout << "dump_path= " << dump_path<< endl;
		flag_dump = 1;
		break;
		case NLIST_digit+3:
		atom_style = val;
		cout << "atom_style= " << atom_style << endl;
		if(atom_style == "bond"){
				flag_atom_style = 1;
		}else if(atom_style == "angle"){
				flag_atom_style = 2;
		}else if(atom_style == "molecular"){
				flag_atom_style = 3;
		}else if(atom_style == "full"){
				flag_atom_style = 4;
		}else if(atom_style == "atomic"){
			cout << "Fatal error: this program analyses the molecular properties, therefore atomic style are not supported for this particular analysis." << endl;
			cout << "Suggestion: please use other package from our website" << endl;
			exit(2);
		}else{
			cout << "Fatal error: unknown or unsuported atomic style for this analysis" << endl;
			cout << "given atom_type: " << atom_style << endl;
			exit(3);
		}
		break;
		case NLIST_digit+4:
    {
      int ind_flag = 0;
      for(int i = 0; i < val.length(); i++){
        if(val[i] != ' '){
          if(val[i] == 'p' || val[i] == 's' || val[i]=='f'){
            if(val[i] == 'p'){
              flag_per[ind_flag] = 1;
            }else if(val[i]=='s' || val[i]=='f'){
              flag_per[ind_flag]=2;
            }
            ind_flag++;
          }else{
            cout << "Error: "<< endl;
            cout << "the given input for the periodicity is not correct" << endl;
            cout << "the given values for each direction must be from {\'s\', \'p\',\'f\'}" << endl;
            cout << "given input: "<< val << endl;
            exit(16);
          }
        }
      }
      if(flag_per[0] == 1){
        cout << "periodic in x direction" << endl;
      }else{
        cout << "non-periodic in x direction" << endl;
      }
      if(flag_per[1] == 1){
        cout << "periodic in y direction" << endl;
      }else{
        cout << "non-periodic in y direction" << endl;
      }
      if(flag_per[2] == 1){
        cout << "periodic in z direction" << endl;
      }else{
        cout << "non-periodic in z direction" << endl;
      }
      break;
    }
	}
}
void compare_assign(int caseNr, double val){
	switch (caseNr) {
		case 0:	//LX = val;
		num_times = (int)val;
		cout << "num_times= " << num_times << endl;
		flag_num_times = 1;
		break;
		case 1:
		intrvl = (int) val;
		cout << "intrvl= " << intrvl<< endl;
		flag_intrvl=1;
		break;
		case 2:
		NRPART = (int) val;
		cout << "NRPART= " << NRPART<< endl;
		flag_NRPART = 1;
		break;
		case 3:
		NRBONDS =  (int) val;
		cout << "NRBONDS= " << NRBONDS<< endl;
		flag_nr_bond = 1;
		break;
		case 4:
		NR_ATOM_TYPE = (int)val;
		cout << "NR_ATOM_TYPE= " << NR_ATOM_TYPE << endl;
		flag_NR_ATOM_TYPE = 1;
		break;
	}
}
void readInput(char *fname){
	ifstream inFile(fname);
	if(!inFile){
		cout << "the initial input file could not be found" << endl;
		cout << "file name: "<< fname << endl;
		exit(1);
 	}
	string str;
	double val;
	int i = 0;//line number
	int counter = 0;// number of correct inputs
	size_t spac;
	string sub1;
	string sub2;
	while(!inFile.eof()){
	 	getline(inFile, str);
		spac = str.find("=");
		sub1 = str.substr(0,spac);
		sub2 = str.substr(spac+1);
		int flag = 0;
		for(int j = 0; j < NLIST; j++){
			//if(sub1.compare("NR_ATOM_TYPE") == 0)cout <<"here0" << endl;
			if(sub1.compare(listInput[j]) == 0){
				if(j < NLIST_digit){
					val = stod(sub2);
					compare_assign(j,val);
				}else{
					compare_assign(j,sub2);
				}
				flag = 1;
				counter++;
			 }
		}
		if(flag == 0 && str.length() != 0){
			cout << "length: " << str.length() << endl;
			cout << "The input file at line " << i+1 << " is not defined well" << endl;
			cout << "Take a look!" << endl;
			cout << "str: " << str << endl;
			spac = str.find("*");
			sub1 = str.substr(0,spac);
			sub2 = str.substr(spac+1);
			cout << "res: " << sub1.compare("NR_ATOM_TYPE") << endl;
			cout << "sub1: " << sub1 << " sub2:" << sub2 << endl;
			cout << "spac: "<< spac << endl;
			exit(2);
		}
		i++;
	 }
}
void apply_pbc_0(int time, Vec& mol, Vec& out)
{
	out = mol;
	if(flag_per[0] == 1){
		out.x = mol.x - ((int) floor((mol.x - XMIN) / LX)) * LX;
		if ( out.x < XMIN ) out.x += LX ;
	}
	if(flag_per[1] == 1){
		out.y = mol.y - ((int) floor((mol.y - YMIN) / LY) ) * LY;
		if ( out.y < YMIN ) out.y += LY ;
	}
	if(flag_per[2] == 1){
		out.z = mol.z - ((int) floor((mol.z - ZMIN) / LZ))* LZ;
		if ( out.z < ZMIN ) out.z += LZ ;
	}
}
void apply_pbc(Vec& mol)
{
	//out = mol;
	if(flag_per[0] == 1){
		mol.x = mol.x - ((int) floor((mol.x - XMIN) / LX)) * LX;
		if( mol.x < XMIN ) mol.x += LX;
	}
	if(flag_per[1] == 1){
		mol.y = mol.y - ((int) floor((mol.y - YMIN) / LY)) * LY;
		if( mol.y < YMIN ) mol.y += LY;
	}
	if(flag_per[2] == 1){
		mol.z = mol.z - ((int) floor((mol.z - ZMIN) / LZ)) * LZ;
		if( mol.z < ZMIN ) mol.z += LZ;
	}
}
int compare_everyPhrase(string satz, string phrase){
	//int nr = 0;
	size_t spac = 0;
	string subOld;
	while(satz.length() != 0){
		spac = satz.find(" ");//, str.length()
		string sub1 = satz.substr(0,spac);
		string sub2 = satz.substr(spac+1);
		//cout << "position of return: " << sub2.find("\n") << endl;
		if(sub1 == phrase){
			return 1;
		}else if(sub2 == phrase){
			return 1;
		}
		satz = sub2;
		if(sub2 == subOld){
			break;
		}
		subOld = sub2;
	}
	return 0;
}
int analys_chID(int *chID){
	int max = chID[0] ;
	for(int i = 1; i < NRPART; i++){
		if(chID[i] > max){
			max = chID[i];
		}
	}
	return max;
}
void makeUpTheChains(chain *& chains, int *chID){
	for(int i = 0; i < NRChains; i++){
		for(int j = 0; j < nrBead_perChain; j++){
			chains[i].ind[j]= -1;
		}
		chains[i].nrAdded = 0;
	}
	for(int i = 0 ; i < NRPART; i++){

		if(chID[i] > NRChains || chID[i] < 0){
			cout << "Fatal error: the chain id out of the range" << endl;
			cout << "i: " << i << ", chID: " << chID[i] - 1 << endl;
			exit(2);
		}
		int nr;
		nr = chains[chID[i]-1].nrAdded;
		if(nr < 0 || nr >= nrBead_perChain){
			cout << "the number of added particles out of range" << endl;
			cout << "chain: " << chID[i] - 1 << ", nrAdded: " << nr << endl;
		}
		chains[chID[i]-1].ind[nr] = i;
		chains[chID[i]-1].nrAdded++;
	}
}
void loadData0( Vec *pos, Vec *posU, int *typen, int **bonds, int *chID, double *masses)
{
	cout << "in the load function" << endl;
	int parInd;
	double parCharge;
	int ix,iy,iz;
	//sprintf(fname, "d.%d",t);
	//ofstream test("testInp.dat");
	const int size_string = init_data.length();
	char fname[size_string+1];
	strcpy(fname, init_data.c_str());
	ifstream inputFile(fname);
	if (!inputFile)
	{
		cout << "ERROR: no input file: " << fname << endl;
		exit(1);
	}
	//ofstream testFile("test-input.dat");
	//string temp;
	//double tempD;
	//int chID[NRPART];
	cout << "loading the header" << endl;
	int flag_found = 0;
	string str, sub1, sub2;
	size_t spac;

	//if(!flag_LX){
	while(!flag_found && !inputFile.eof()){
			getline(inputFile, str);
			if(compare_everyPhrase(str,"xlo xhi")){
				cout << "str: " << str << endl;
				spac = str.find(" ");//, str.length()
				sub1 = str.substr(0,spac);
				sub2 = str.substr(spac+1);
				cout << "sub1: " << sub1 << endl;
				XMIN_0 = stod(sub1);
				str = sub2;
				spac = str.find(" ");//, str.length()
				sub1 = str.substr(0,spac);
				sub2 = str.substr(spac+1);
				cout << "sub1: " << sub1 << endl;
				XMAX_0 = stod(sub1);
				flag_found = 1;
		}
	}
	//}
	LX_0 = XMAX_0 - XMIN_0;
	flag_found = 0;
//	if(!flag_LY){
	while(!flag_found && !inputFile.eof()){
			getline(inputFile, str);
			if(compare_everyPhrase(str,"ylo yhi")){
				cout << "str: " << str << endl;
				spac = str.find(" ");//, str.length()
				sub1 = str.substr(0,spac);
				sub2 = str.substr(spac+1);
				cout << "sub1: " << sub1 << endl;
				YMIN_0 = stod(sub1);
				str = sub2;
				spac = str.find(" ");//, str.length()
				sub1 = str.substr(0,spac);
				sub2 = str.substr(spac+1);
				cout << "sub1: " << sub1 << endl;
				YMAX_0 = stod(sub1);
				flag_found = 1;
		}
	}
//	}
	LY_0 = YMAX_0 - YMIN_0;
	flag_found = 0;
	//if(!flag_LZ){
	while(!flag_found && !inputFile.eof()){
			getline(inputFile, str);
			if(compare_everyPhrase(str,"zlo zhi")){
				cout << "str: " << str << endl;
				spac = str.find(" ");//, str.length()
				sub1 = str.substr(0,spac);
				sub2 = str.substr(spac+1);
				cout << "sub1: " << sub1 << endl;
				ZMIN_0 = stod(sub1);
				str = sub2;
				spac = str.find(" ");//, str.length()
				sub1 = str.substr(0,spac);
				sub2 = str.substr(spac+1);
				cout << "sub1: " << sub1 << endl;
				ZMAX_0 = stod(sub1);
				flag_found = 1;
		}
	}
	//}
	LZ_0 = ZMAX_0 - ZMIN_0;
	cout << "x: " << XMIN_0 << " " << XMAX_0 << endl;
	cout << "y: " << YMIN_0 << " " << YMAX_0 << endl;
	cout << "z: " << ZMIN_0 << " " << ZMAX_0 << endl;
	LXINV_0 = 1./LX_0;
	LYINV_0 = 1./LY_0;
	LZINV_0 = 1./LZ_0;

	if(flag_NR_ATOM_TYPE){
		flag_found = 0;
		while(!flag_found && !inputFile.eof()){
			getline(inputFile, str);
			spac = str.find(" ");
			sub1 = str.substr(0,spac);
			sub2 = str.substr(spac+1);
			if(spac < str.length()){
				if(sub1 == "Masses"){
					cout << "flag found" << endl;
					cout << "str: " << str << endl;
					flag_found = 1;
				}
			}else{
				if(str == "Masses"){
					cout << "flag found" << endl;
					cout << "str: "<< str << endl;
					flag_found = 1;
				}
			}
		}
		getline(inputFile, str);
		int tmp;
		cout << "NR_ATOM_TYPE: " << NR_ATOM_TYPE << endl;
		for(int i = 0; i < NR_ATOM_TYPE; i++){
			inputFile >> tmp;
			if(tmp > NR_ATOM_TYPE){
				cout << "fatal error: atom_type out of range" << endl;
				exit(3);
			}
			inputFile >> masses[tmp - 1];
			cout << tmp-1 << ": " << masses[tmp - 1] << endl;
			//testFile << masses[tmp - 1] << endl;
		}
	}
	flag_found = 0;
	while(!flag_found && !inputFile.eof()){
		getline(inputFile, str);
		spac = str.find(" ");
		sub1 = str.substr(0,spac);
		sub2 = str.substr(spac+1);
		if(spac < str.length()){
			if(sub1 == "Atoms"){
				cout << "flag found" << endl;
				cout << "str: " << str << endl;
				flag_found = 1;
			}
		}else{
			if(str == "Atoms"){
				cout << "flag found" << endl;
				cout << "str: "<< str << endl;
				flag_found = 1;
			}
		}
	}
	getline(inputFile, str);
	//getline(inputFile, str);
	//cout << str << endl;
  //int nochTmp, boxInd;
	inputFile.precision(5);
	cout << "loading the particles" << endl;
	cout << "atom_style: " << atom_style << endl;
	for(int i = 0; i < NRPART; i++){

		if(flag_atom_style < 3){
			inputFile >> parInd ;
			if(parInd > NRPART || parInd < 0){
				cout << "problem in reading the atomic data" << endl;
				cout << "particle index exceeds the expected values" << endl;
				cout << "parInd: "<< parInd << endl;
				exit(16);
			}
			// if(chInd >= NRChains || chInd < 0){
			// 	cout << "problem in reading the atomic data" << endl;
			// 	cout << "chain index exceeds the expected values" << endl;
			// 	cout << "chInd: "<< chInd << endl;
			// 	exit(16);
			// }
			inputFile >> chID[parInd - 1] >> typen[parInd - 1] >> pos[parInd - 1].x >> pos[parInd - 1].y >> pos[parInd - 1].z >> ix >> iy >> iz;
			posU[parInd - 1].x = XMIN_0 + ix*LX_0;
			posU[parInd - 1].y = YMIN_0 + iy*LY_0;
			posU[parInd - 1].z = ZMIN_0 + iz*LZ_0;
			typen[parInd - 1]--;
			//chains[chInd - 1].addOne(parInd - 1);
			//testFile << parInd << " " << chID[parInd - 1] << " " << typen[parInd - 1] << " " << pos[parInd - 1] << " " << ix << " " << iy << " " << iz << endl;
		}else if(flag_atom_style >=3){
			inputFile >> parInd;// >> chInd;
			if(parInd > NRPART || parInd < 0){
				cout << "problem in reading the atomic data" << endl;
				cout << "particle index exceeds the expected values" << endl;
				cout << "parInd: "<< parInd << endl;
				exit(16);
			}
			/*if(chInd >= NRChains || chInd < 0){
				cout << "problem in reading the atomic data" << endl;
				cout << "chain index exceeds the expected values" << endl;
				cout << "chInd: "<< chInd << endl;
				exit(16);
			}*/
			inputFile >> chID[parInd - 1] >> typen[parInd - 1] >> parCharge >> pos[parInd - 1].x >> pos[parInd - 1].y >> pos[parInd - 1].z >> ix >> iy >> iz;
			posU[parInd - 1].x = XMIN_0 + ix*LX_0;
			posU[parInd - 1].y = YMIN_0 + iy*LY_0;
			posU[parInd - 1].z = ZMIN_0 + iz*LZ_0;
			//testFile << parInd << " " << chID[parInd - 1] << " " << typen[parInd - 1] << " " << parCharge << " " << pos[parInd - 1] << " " << ix << " " << iy << " " << iz << endl;
			//chains[chInd - 1].addOne(parInd - 1);
		}
		// for(int i = 0; i < NRChains; i++){
		// 	if(chains[i].nrAdded < nrBead_perChain - 1){
		// 		cout << "Fatal error: there is a problem with the data file" << endl;
		// 		cout << "at least one chain has not been filled with particles" << endl;
		// 		exit(17);
		// 	}
		// }
	}
	cout << "particles position loaded" << endl;
	// getline(inputFile, str);
	// getline(inputFile, str);
	// cout << "str: " << str << endl;

	//int nline = 0;
	flag_found = 0;
	while(!flag_found){
		getline(inputFile, str);
		if(inputFile.eof()){
			cout << "end of file is reached" << endl;
			break;
		}
		spac = str.find(" ");
		sub1 = str.substr(0,spac);
		sub2 = str.substr(spac+1);
		if(spac < str.length()){
			if(sub1 == "Bonds"){
				cout << "flag found" << endl;
				cout << "str: " << str << endl;
				flag_found = 1;
			}
		}
		if(str == "Bonds"){
			cout << "flag found" << endl;
			cout << "str: "<< str << endl;
			flag_found = 1;
		}
	}
	if(!flag_found){
		cout << "fatal error: the bonds are not found in the data file" << endl;
		exit(3);
	}
	//testFile << "bonds"<< endl;
	getline(inputFile, str);
	//getline(inputFile, str);
	//cout << str << endl;
	int ind;
	for(int i = 0; i < NRBONDS; i++){
		inputFile >> ind;
		if(ind > NRBONDS || ind < 0){
			cout << "the bond id is out of range" << endl;
			cout << "ind: " << ind << endl;
			exit(64);
		}
		inputFile >> bonds[ind - 1][0] >> bonds[ind - 1][1] >> bonds[ind - 1][2];
		//testFile << bonds[ind - 1][0] << " " << bonds[ind - 1][1] << " " << bonds[ind - 1][2]<< endl;
	}
	//testFile.close();
  inputFile.close();
}
void loadDump(int t, Vec* pos, int *typen, Vec* posU)
{
	int i, iniLine,parInd, molInd;
	const int size_string = dump_path.length();
	char fname0[size_string+1], fname[size_string+100];
	strcpy(fname0, dump_path.c_str());
	sprintf(fname,"%s/dump.%d",fname0,t);
	ifstream inputFile(fname);
	if (!inputFile)
	{
		cout << "ERROR: no input file: " << fname << endl;
		exit(1);
	}
	string str;
	//double tempD;
	if(dynamic_box){
		iniLine = 5;
		for(i = 0; i < iniLine; i++){
			getline(inputFile,str);
		}
		inputFile >> XMIN >> XMAX;
		inputFile >> YMIN >> YMAX;
		inputFile >> ZMIN >> ZMAX;
		LX = XMAX - XMIN;
		LY = YMAX - YMIN;
		LZ = ZMAX - ZMIN;
		LXINV = 1. / LX;
		LYINV = 1. / LY;
		LZINV = 1. / LZ;
		cout << "LiS:\n" << LX << "\n" << LY << "\n" << LZ << endl;
		getline(inputFile,str);
		getline(inputFile,str);
	}else{
		iniLine = 9;
		for(i = 0; i < iniLine; i++){
			getline(inputFile,str);
		}
	}
	char fname5[size_string+100];
	sprintf(fname5,"%s/tDump.%d",fname0,t);
	ofstream oFile(fname5);
	//spring testDump("tDump.")
	cout << "Very Important: please note that the following structure of dump file is expected to be like: " << endl;
	cout << "\tid mol atom_type xu yu zu" << endl;
	inputFile.precision(5);
	for(i = 0; i < NRPART; i++){
		// getline(inputFile,str);
		// oFile << str << endl;
		inputFile >> parInd >> molInd;

		if(parInd > NRPART || parInd < 0){
				cout << "problem in reading the atomic data" << endl;
				cout << "particle index exceeds the expected values" << endl;
				cout << "parInd: "<< parInd << endl;
				exit(16);
		}
		inputFile >> typen[parInd - 1] >> pos[parInd - 1].x >> pos[parInd - 1].y >> pos[parInd - 1].z;

		posU[parInd - 1] = pos[parInd - 1];
		apply_pbc(pos[parInd - 1]);
		oFile << parInd - 1 << " "<< pos[parInd -  1] << endl;
	}
	oFile.close();
	inputFile.close();
}
void store_init_cen_masses(chain *chains, Vec *cen0, Vec *posU, double *masses, int *typen){
	for(int i = 0; i < NRChains; i++){
		cen0[i]= chains[i].calCenMass(posU, masses, typen);
		cout << i << " " << cen0[i] << endl;
	}
}
Vec calculate_msd(chain *chains,Vec *posU, Vec *cen0, double *masses, int *typen){
	Vec msd;//.reset();
	msd.reset();
	for(int i = 0; i < NRChains; i++){
		Vec cen;
		cen = chains[i].calCenMass(posU, masses, typen);
		msd.x += pow((cen.x - cen0[i].x),2.);
		msd.y += pow((cen.y - cen0[i].y),2.);
		msd.z += pow((cen.z - cen0[i].z),2.);
	}
	msd = msd *(1./(double)NRChains);
	return msd;
}
Vec calculate_Rg(chain *chains, Vec *posU, double *masses, int *typen, Vec *&rG_chain){
	Vec rg;//.reset();
	rg.reset();
	ofstream testFile("testInd.dat");
	for(int i = 0; i < NRChains; i++){
		chains[i].calRg(posU,masses,typen);
		rG_chain[i] = chains[i].rG;
	//	cout << "rg "<< i << " " << rG_chain[i] << " " << chains[i].nrBeads << endl;
		rg += chains[i].rG;
	}
	rg = rg * (1./(double)NRChains);
	return rg;
}
void calEndToEnd(chain *chains, Vec *posU, Vec *rEE_chain){
	for(int i = 0; i < NRChains; i++){
		rEE_chain[i] = chains[i].calEndToEnd(posU);
	}
}
void calHistRadGyr(Vec **radGyr){
	Vec min_rg=radGyr[0][0], max_rg=radGyr[0][0];
	for(int t = 0; t < num_times; t++){
		for(int i = 0; i < NRChains; i++){
			if(radGyr[t][i].x < min_rg.x){
				min_rg.x = radGyr[t][i].x;
			}
			if(radGyr[t][i].x > max_rg.x){
				max_rg.x = radGyr[t][i].x;
			}
			if(radGyr[t][i].y < min_rg.y){
				min_rg.y = radGyr[t][i].y;
			}
			if(radGyr[t][i].y > max_rg.y){
				max_rg.y = radGyr[t][i].y;
			}
			if(radGyr[t][i].z < min_rg.z){
				min_rg.z = radGyr[t][i].z;
			}
			if(radGyr[t][i].z > max_rg.z){
				max_rg.z = radGyr[t][i].z;
			}
		}
	}
	if(min_rg.x > max_rg.x){
		cout << "fatal error: there is an error with finding extermums in calculation of histograms: X" << endl;
		exit(89);
	}
	if(min_rg.y > max_rg.y){
		cout << "fatal error: there is an error with finding extermums in calculation of histograms: Y" << endl;
		exit(90);
	}
	if(min_rg.z > max_rg.z){
		cout << "fatal error: there is an error with finding extermums in calculation of histograms: Z" << endl;
		exit(91);
	}
	Vec del_hist = (max_rg - min_rg)*(1./(double)nr_hist_bins);
	Vec *hist;
	alloc_Vec(hist, nr_hist_bins);
	for(int i = 0; i < nr_hist_bins; i++){
		hist[i].reset();
	}
	int binX, binY, binZ;
	for(int t = 0; t < num_times; t++){
		for(int i = 0; i < NRChains; i++){
			binX = (int)floor((radGyr[t][i].x-min_rg.x)/del_hist.x);
			if(binX < 0 || binX >= nr_hist_bins){
				cout << "fatal error: in calculation of histogram, the bin out of range: x" << endl;
				cout << "bin: " << binX << endl;
				exit(92);
			}
			hist[binX].x ++;
			binY = (int)floor((radGyr[t][i].y-min_rg.y)/del_hist.y);
			if(binY < 0 || binY >= nr_hist_bins){
				cout << "fatal error: in calculation of histogram, the bin out of range: y" << endl;
				cout << "bin: " << binY << endl;
				exit(93);
			}
			hist[binY].y ++;
			binZ = (int)floor((radGyr[t][i].z-min_rg.z)/del_hist.z);
			if(binZ < 0 || binZ >= nr_hist_bins){
				cout << "fatal error: in calculation of histogram, the bin out of range: z" << endl;
				cout << "bin: " << binZ << endl;
				exit(92);
			}
			hist[binZ].z ++;
		}
	}
	ofstream oFileX("hist_radGyr_x.dat"),oFileY("hist_radGyr_y.dat"),oFileZ("hist_radGyr_z.dat");
	for(int i = 0; i < nr_hist_bins; i++){
		hist[i].x = hist[i].x / (double)(NRChains*num_times);
		oFileX << min_rg.x+(i+0.5)*del_hist.x << " " << hist[i].x << endl;
		hist[i].y = hist[i].y / (double)(NRChains*num_times);
		oFileY << min_rg.y+(i+0.5)*del_hist.y << " " << hist[i].y << endl;
		hist[i].z = hist[i].z / (double)(NRChains*num_times);
		oFileZ << min_rg.z+(i+0.5)*del_hist.z << " " << hist[i].z << endl;
	}
	oFileX.close();
	oFileY.close();
	oFileZ.close();
	free(hist);
}
void store_time_arrays( string fname0, Vec *arr, int *t_series ){
	fname0 = dump_path+fname0;
	ofstream oFile(fname0);
	for(int t = 0; t < num_times; t++){
		oFile << t_series[t] << " " << arr[t] << endl;
	}
	oFile.close();
}
void calEndToEndCorr(Vec **rEndEnd_t, int *t_series){
	double *corr(NULL), denom=0.;
	alloc_double(corr, num_times);
	const int size_string = dump_path.length();
	char fname0[size_string+1], fname[size_string+100];
	strcpy(fname0, dump_path.c_str());
	sprintf(fname,"%s/endToEnd_corr.dat",fname0);
	ofstream oFile(fname);
	for(int t = 0; t < num_times; t++){
		corr[t]= 0.;
		for(int i = 0; i < NRChains; i++){
			corr[t] += rEndEnd_t[t][i].innerProd(rEndEnd_t[0][i]);
			if(t == 0){
				denom += rEndEnd_t[t][i].innerProd(rEndEnd_t[0][i]);
			}
		}
		denom /= (double) NRChains;
		corr[t] /= (double) NRChains;
		corr[t] = corr[t] / denom;
		oFile << t_series[t] << " " << corr[t] << endl;
	}
	oFile.close();
	free(corr);
}
