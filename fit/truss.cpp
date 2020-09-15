#include"truss.hpp"

//Finite element analysis and structure modelling were developed based on:
//https://github.com/labicom-udesc/ECO-framework/blob/master/ProblemFunc.c

int num_elem;                                           
int num_node; 
int dim;

double *coordx; 
double *coordy;  
double *coordz;
double **fX;      
double **fY;     
double **fZ;
double *ymod;

int **member;

double *length;

int **kode;              

double E;
double ro;

double maxDisp;
double **maxStress;
double **minStress;

double minArea;
double maxArea;

double maxW;
double minW;

void (*distributeAreas)(double*, double*, int);
double (*fitChoice)(double*, double*, double*, bool, bool*, int);

void setupTruss(int nTruss, int nThreads, int *avalPerInd){
	switch (nTruss){
		case 10:
			setupTruss10(nThreads, avalPerInd);
			break;
		case 17:
			setupTruss17(nThreads, avalPerInd);
			break;
		case 25:
			setupTruss25(nThreads, avalPerInd);
			break;
		case 721:
			setupTruss72case1(nThreads, avalPerInd);
			break;
		case 722:
			setupTruss72case2(nThreads, avalPerInd);
			break;
		case 120:
			setupTruss120(nThreads, avalPerInd);
			break;
		case 200:
			setupTruss200(nThreads, avalPerInd);
			break;
	}
}

void setupTruss10(int nThreads, int *avalPerInd){
    num_elem = 10;                                           //GLOBAL
    num_node = 6;                                            //GLOBAL
    dim = 10;
	
    E = 1.0e7;
    double load = 1.0e5;
    ro = 0.1;

	*avalPerInd = 1;

	maxDisp = 2.0;
    maxStress = (double**)malloc(num_elem*sizeof(double*));
    minStress = (double**)malloc(num_elem*sizeof(double*));

	double *maxStressData = (double*)malloc(num_elem*nThreads*sizeof(double));
	double *minStressData = (double*)malloc(num_elem*nThreads*sizeof(double));

	for(int i=0;i<num_elem;i++){
		maxStress[i] = maxStressData + i*nThreads;
		minStress[i] = minStressData + i*nThreads;
	}

	for(int j=0;j<nThreads;j++){
		for(int i=0;i<num_elem;i++){
			maxStress[i][j] = 25.0e3;
			minStress[i][j] = -25.0e3;
		}
	}

	distributeAreas = distributeAreasEqual;
    fitChoice = trussFit2D;
    
    kode = (int**)malloc(num_node*sizeof(int*));              //GLOBAL
    int *kodeData = (int*)calloc(num_node*2, sizeof(int));
    for(int i=0;i<num_node;i++){
        kode[i] = kodeData + i*2;
    }
    
    coordx = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    coordy = (double*)malloc(num_node*sizeof(double));   //GLOBAL

    fX = (double**)malloc(num_node*sizeof(double*));      //GLOBAL
    fY = (double**)malloc(num_node*sizeof(double*));       //GLOBAL

	double *fXdata = (double*)calloc(num_node*nThreads, sizeof(double));
	double *fYdata = (double*)calloc(num_node*nThreads, sizeof(double));
	for(int i=0;i<num_node;i++){
		fX[i] = fXdata + i*nThreads;
		fY[i] = fYdata + i*nThreads;
	}
    
    coordz = NULL;
    fZ = NULL;
    
    int *memberData;
    member = (int**)malloc(num_elem*sizeof(int*));               //GLOBAL
    memberData = (int*)malloc(num_elem*2*sizeof(int));
    for(int i=0;i<num_elem;i++){
        member[i] = memberData + i*2;
    }
    
    length = (double*)malloc(num_elem*sizeof(double));   //GLOBAL
    
    ymod = (double*)malloc(num_elem*sizeof(double));
    
    /* fixo X e Y*/
    coordx[0] = 720.0;
	coordy[0] = 360.0;
    coordx[1] = 720.0;
	coordy[1] = 0.0;
    coordx[2] = 360.0;
	coordy[2] = 360.0;
    coordx[3] = 360.0;
	coordy[3] = 0.0;
    coordx[4] = 0.0;
	coordy[4] = 360.0;
    coordx[5] = 0.0;
	coordy[5] = 0.0;
	
	for(int i=0;i<nThreads;i++){
		fY[1][i] = -load;
		fY[3][i] = -load;
	}
	
    kode[4][0] = 1;
    kode[4][1] = 1;
    kode[5][0] = 1;
    kode[5][1] = 1;

	member[0][0] = 5; member[0][1] = 3;
	member[1][0] = 3; member[1][1] = 1;
	member[2][0] = 6; member[2][1] = 4;
	member[3][0] = 4; member[3][1] = 2;
	member[4][0] = 3; member[4][1] = 4;
	member[5][0] = 1; member[5][1] = 2;
	member[6][0] = 5; member[6][1] = 4;
	member[7][0] = 6; member[7][1] = 3;
	member[8][0] = 3; member[8][1] = 2;
	member[9][0] = 1; member[9][1] = 4;

	ymod[0] = E;
	ymod[1] = E;
	ymod[2] = E;
	ymod[3] = E;
	ymod[4] = E;
	ymod[5] = E;
	ymod[6] = E;
	ymod[7] = E;
	ymod[8] = E;
	ymod[9] = E;
    
    //compute member lengths
    for(int i=0;i<num_elem;i++){
        int node1, node2;
        
        node1 = member[i][0]-1;
        node2 = member[i][1]-1;
        
        length[i] = sqrt( pow(coordx[node1]-coordx[node2], 2) + pow(coordy[node1]-coordy[node2], 2) );
    }

	//min and max areas
	minArea = 0.1;
	maxArea = 35.0; 

	//compute max and min weights
	maxW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		maxW += ro*maxArea*length[i];
	}

	minW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		minW += ro*minArea*length[i];
	}
}

void setupTruss17(int nThreads, int *avalPerInd){
    num_elem = 17;                                           //GLOBAL
    num_node = 9;                                            //GLOBAL
	dim = 17;
    
    E = 3.0e7;
    double load = 1.0e5;
    ro = 0.268;

	*avalPerInd = 1;

	maxDisp = 2.0;

	maxStress = (double**)malloc(num_elem*sizeof(double*));
    minStress = (double**)malloc(num_elem*sizeof(double*));

	double *maxStressData = (double*)malloc(num_elem*nThreads*sizeof(double));
	double *minStressData = (double*)malloc(num_elem*nThreads*sizeof(double));

	for(int i=0;i<num_elem;i++){
		maxStress[i] = maxStressData + i*nThreads;
		minStress[i] = minStressData + i*nThreads;
	}

	for(int j=0;j<nThreads;j++){
		for(int i=0;i<num_elem;i++){
			maxStress[i][j] = 50.0e3;
			minStress[i][j] = -50.0e3;
		}
	}

	distributeAreas = distributeAreasEqual;
    fitChoice = trussFit2D;
    
    kode = (int**)malloc(num_node*sizeof(int*));              //GLOBAL
    int *kodeData = (int*)calloc(num_node*2, sizeof(int));
    for(int i=0;i<num_node;i++){
        kode[i] = kodeData + i*2;
    }
    
    coordx = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    coordy = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    
	fX = (double**)malloc(num_node*sizeof(double*));      //GLOBAL
    fY = (double**)malloc(num_node*sizeof(double*));       //GLOBAL

	double *fXdata = (double*)calloc(num_node*nThreads, sizeof(double));
	double *fYdata = (double*)calloc(num_node*nThreads, sizeof(double));
	for(int i=0;i<num_node;i++){
		fX[i] = fXdata + i*nThreads;
		fY[i] = fYdata + i*nThreads;
	}
    
    coordz = NULL;
    fZ = NULL;
    
    int *memberData;
    member = (int**)malloc(num_elem*sizeof(int*));               //GLOBAL
    memberData = (int*)malloc(num_elem*2*sizeof(int));
    for(int i=0;i<num_elem;i++){
        member[i] = memberData + i*2;
    }
    
    length = (double*)malloc(num_elem*sizeof(double));   //GLOBAL
    
    ymod = (double*)malloc(num_elem*sizeof(double));
    
    /* fixo X e Y*/
    coordx[0] = 400.0;
	coordy[0] = 0.0;
    coordx[1] = 300.0;
	coordy[1] = 100.0;
    coordx[2] = 300.0;
	coordy[2] = 0.0;
    coordx[3] = 200.0;
	coordy[3] = 100.0;
    coordx[4] = 200.0;
	coordy[4] = 0.0;
    coordx[5] = 100.0;
	coordy[5] = 100.0;
    coordx[6] = 100.0;
	coordy[6] = 0.0;
    coordx[7] = 0.0;
	coordy[7] = 100.0;
    coordx[8] = 0.0;
	coordy[8] = 0.0;
    
	for(int i=0;i<nThreads;i++){
		fY[0][i] = -load;
	}
    
    kode[7][0] = 1;
    kode[7][1] = 1;
    kode[8][0] = 1;
    kode[8][1] = 1;

	member[0][0]  = 6; member[0][1] = 8;
	member[1][0]  = 7; member[1][1] = 8;
	member[2][0]  = 7; member[2][1] = 9;
	member[3][0]  = 6; member[3][1] = 7;
	member[4][0]  = 4; member[4][1] = 6;
	member[5][0]  = 5; member[5][1] = 6;
	member[6][0]  = 5; member[6][1] = 7;
	member[7][0]  = 4; member[7][1] = 5;
	member[8][0]  = 2; member[8][1] = 4;
	member[9][0]  = 3; member[9][1] = 4;
    member[10][0] = 3; member[10][1] = 5;
    member[11][0] = 2; member[11][1] = 3;
    member[12][0] = 1; member[12][1] = 2;
    member[13][0] = 1; member[13][1] = 3;
    member[14][0] = 6; member[14][1] = 9;
    member[15][0] = 4; member[15][1] = 7;
    member[16][0] = 2; member[16][1] = 5;

	for(int i=0;i<num_elem;i++){
        ymod[i] = E;
    }
    
    //compute member lengths
    for(int i=0;i<num_elem;i++){
        int node1, node2;
        
        node1 = member[i][0]-1;
        node2 = member[i][1]-1;
        
        length[i] = sqrt( pow(coordx[node1]-coordx[node2], 2) + pow(coordy[node1]-coordy[node2], 2) );
    }

	//min and max areas
	minArea = 0.1;
	maxArea = 50.0;

	//compute max and min weights
	maxW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		maxW += ro*maxArea*length[i];
	}

	minW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		minW += ro*minArea*length[i];
	}
}

void setupTruss200(int nThreads, int *avalPerInd){
	num_elem = 200;                                           //GLOBAL
    num_node = 77;                                            //GLOBAL
	dim = 29;
    
	E = 3.0e4; //ksi

	*avalPerInd = 3; //number of structural analysis per individual

    ro = 0.283;

	maxDisp = 50.0;
    
	maxStress = (double**)malloc(num_elem*sizeof(double*));
    minStress = (double**)malloc(num_elem*sizeof(double*));

	double *maxStressData = (double*)malloc(num_elem*nThreads*sizeof(double));
	double *minStressData = (double*)malloc(num_elem*nThreads*sizeof(double));

	for(int i=0;i<num_elem;i++){
		maxStress[i] = maxStressData + i*nThreads;
		minStress[i] = minStressData + i*nThreads;
	}

	for(int j=0;j<nThreads;j++){
		for(int i=0;i<num_elem;i++){
			maxStress[i][j] = 10.0;  //ksi
			minStress[i][j] = -10.0; //ksi
		}
	}

	distributeAreas = distributeAreas200;
    fitChoice = fit200;
    
    //kode = (int*)calloc(num_node, sizeof(int));              //GLOBAL
    kode = (int**)malloc(num_node*sizeof(int*));              //GLOBAL
    int *kodeData = (int*)calloc(num_node*2, sizeof(int));
    for(int i=0;i<num_node;i++){
        kode[i] = kodeData + i*2;
    }
    
    coordx = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    coordy = (double*)malloc(num_node*sizeof(double));   //GLOBAL
	
	fX = (double**)malloc(num_node*sizeof(double*));      //GLOBAL
    fY = (double**)malloc(num_node*sizeof(double*));       //GLOBAL

	double *fXdata = (double*)calloc(num_node*nThreads, sizeof(double));
	double *fYdata = (double*)calloc(num_node*nThreads, sizeof(double));
	for(int i=0;i<num_node;i++){
		fX[i] = fXdata + i*nThreads;
		fY[i] = fYdata + i*nThreads;
	}
    
    coordz = NULL;
    fZ = NULL;
    
    int *memberData;
    member = (int**)malloc(num_elem*sizeof(int*));               //GLOBAL
    memberData = (int*)malloc(num_elem*2*sizeof(int));
    for(int i=0;i<num_elem;i++){
        member[i] = memberData + i*2;
    }
    
    length = (double*)malloc(num_elem*sizeof(double));   //GLOBAL
    
    ymod = (double*)malloc(num_elem*sizeof(double));
    
    //Y node coordinates
    int store = 10;
	for(int i=0;i<5;i++){
		for(int j=14*i;j < 14*i + 5;j++){
			coordy[j] = 360.0 + store*144.0;
		}
		store-=1;
		for(int j=14*i + 5;j < 14*i + 14;j++){
			coordy[j] = 360.0 + store*144.0;
		}
		store-=1;
	}
	for(int j=70;j<75;j++){
		coordy[j] = 360.0;
	}
	coordy[75] = 0.0;
	coordy[76] = 0.0;

	//X node coordinates
	for(int i=0;i<5;i++){
		for(int j=0;j<5;j++){
			coordx[14*i + j] = 240.0*j;
		}
		for(int j=0;j<9;j++){
			coordx[14*i + 5 + j] = 120.0*j;
		}
	}
	for(int j=0;j<5;j++){
		coordx[70 + j] = 240.0*j;
	}
	coordx[75] = 240.0;
	coordx[76] = 720.0;
	
    kode[75][0] = 1;
    kode[75][1] = 1;
    kode[76][0] = 1;
    kode[76][1] = 1;

	//node conections (elements)
	
	//horizontal elements
	for(int i=0;i<5;i++){
		for(int j=0;j<4;j++){
			member[38*i + j][0] = 14*i + j + 1;
			member[38*i + j][1] = 14*i + j + 2;
		}
		for(int j=0;j<8;j++){
			member[38*i + 17 + j][0] = 14*i + 5 + j + 1;
			member[38*i + 17 + j][1] = 14*i + 5 + j + 2;
		}
	}
	for(int j=0;j<4;j++){
		member[190 + j][0] = 71+j;
		member[190 + j][1] = 71+j+1;
	}
	//horizontal elements end
	
	//vertical elements
	for(int i=0;i<5;i++){
		for(int j=0;j<5;j++){
			member[38*i + 4 + 3*j][0] = 14*i + j + 1;
			member[38*i + 4 + 3*j][1] = 14*i + 5 + 2*j + 1;

			member[38*i + 25 + +3*j][0] = 14*i + 5 + 2*j + 1;
			member[38*i + 25 + +3*j][1] = 14*i + 14 + j + 1;
		}
	}
	//vertical elements end

	//diagonal elements
	for(int i=0;i<5;i++){
		for(int j=0;j<4;j++){
			member[38*i + 5 + 3*j][0] = 14*i + j + 1;
			member[38*i + 5 + 3*j][1] = 14*i + 7 + 2*j;

			member[38*i + 6 + 3*j][0] = 14*i + j + 2;
			member[38*i + 6 + 3*j][1] = 14*i + 7 + 2*j;

			member[38*i + 26 + 3*j][0] = 14*i + 7 + 2*j;
			member[38*i + 26 + 3*j][1] = 14*i + 15 + j;

			member[38*i + 27 + 3*j][0] = 14*i + 7 + 2*j;
			member[38*i + 27 + 3*j][1] = 14*i + 16 + j;
		}
	}

	for(int j=0;j<3;j++){
		member[194 + j][0] = 71 + j;
		member[194 + j][1] = 76;

		member[197 + j][0] = 73 + j;
		member[197 + j][1] = 77;
	}
	//diagonal elements end

	for(int i=0;i<num_elem;i++){
        ymod[i] = E;
    }
    
    //compute member lengths
    for(int i=0;i<num_elem;i++){
        int node1, node2;
        
        node1 = member[i][0]-1;
        node2 = member[i][1]-1;
        
        length[i] = sqrt( pow(coordx[node1]-coordx[node2], 2) + pow(coordy[node1]-coordy[node2], 2) );
    }

	//min and max areas
	minArea = 0.1;
	maxArea = 20.0;

	//compute max and min weights
	maxW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		maxW += ro*maxArea*length[i];
	}

	minW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		minW += ro*minArea*length[i];
	}
}

void setupTruss25(int nThreads, int *avalPerInd){
    num_elem = 25;                                           //GLOBAL
    num_node = 10;                                            //GLOBAL
	dim = 8;
    
    E = 1.0e4; //ksi

	*avalPerInd = 2;
	
    ro = 0.1;  //lb/in^3

	maxDisp = 0.35;
	maxStress = (double**)malloc(num_elem*sizeof(double*));
    minStress = (double**)malloc(num_elem*sizeof(double*));

	double *maxStressData = (double*)malloc(num_elem*nThreads*sizeof(double));
	double *minStressData = (double*)malloc(num_elem*nThreads*sizeof(double));

	for(int i=0;i<num_elem;i++){
		maxStress[i] = maxStressData + i*nThreads;
		minStress[i] = minStressData + i*nThreads;
	}

	distributeAreas = distributeAreas25;
    fitChoice = fit25;
    
    kode = (int**)malloc(num_node*sizeof(int*));              //GLOBAL
    int *kodeData = (int*)calloc(num_node*3, sizeof(int));
    for(int i=0;i<num_node;i++){
        kode[i] = kodeData + i*3;
    }
    
    coordx = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    coordy = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    coordz = (double*)malloc(num_node*sizeof(double));

	fX = (double**)malloc(num_node*sizeof(double*));
    fY = (double**)malloc(num_node*sizeof(double*));
	fZ = (double**)malloc(num_node*sizeof(double*));

	double *fXdata = (double*)calloc(num_node*nThreads, sizeof(double));
	double *fYdata = (double*)calloc(num_node*nThreads, sizeof(double));
	double *fZdata = (double*)calloc(num_node*nThreads, sizeof(double));
	for(int i=0;i<num_node;i++){
		fX[i] = fXdata + i*nThreads;
		fY[i] = fYdata + i*nThreads;
		fZ[i] = fZdata + i*nThreads;
	}
    
    int *memberData;
    member = (int**)malloc(num_elem*sizeof(int*));               //GLOBAL
    memberData = (int*)malloc(num_elem*2*sizeof(int));
    for(int i=0;i<num_elem;i++){
        member[i] = memberData + i*2;
    }
    
    length = (double*)malloc(num_elem*sizeof(double));   //GLOBAL
    
    ymod = (double*)malloc(num_elem*sizeof(double));
    
    /* fixo X e Y*/
    coordx[0] = -37.5;
	coordy[0] = 0.0;
    coordz[0] = 200.0;
    coordx[1] = 37.5;
	coordy[1] = 0.0;
    coordz[1] = 200.0;
    coordx[2] = -37.5;
	coordy[2] = 37.5;
    coordz[2] = 100.0;
    coordx[3] = 37.5;
	coordy[3] = 37.5;
    coordz[3] = 100.0;
    coordx[4] = 37.5;
	coordy[4] = -37.5;
    coordz[4] = 100.0;
    coordx[5] = -37.5;
	coordy[5] = -37.5;
    coordz[5] = 100.0;
    coordx[6] = -100.0;
	coordy[6] = 100.0;
    coordz[6] = 0.0;
    coordx[7] = 100.0;
	coordy[7] = 100.0;
    coordz[7] = 0.0;
    coordx[8] = 100.0;
	coordy[8] = -100.0;
    coordz[8] = 0.0;
    coordx[9] = -100.0;
	coordy[9] = -100.0;
    coordz[9] = 0.0;
    
    kode[6][0] = 1;
    kode[6][1] = 1;
    kode[6][2] = 1;
    kode[7][0] = 1;
    kode[7][1] = 1;
    kode[7][2] = 1;
    kode[8][0] = 1;
    kode[8][1] = 1;
    kode[8][2] = 1;
    kode[9][0] = 1;
    kode[9][1] = 1;
    kode[9][2] = 1;

	member[0][0] = 1; member[0][1] = 2;
	member[1][0] = 1; member[1][1] = 4;
	member[2][0] = 2; member[2][1] = 3;
    member[3][0] = 1; member[3][1] = 5;
    member[4][0] = 2; member[4][1] = 6;
    member[5][0] = 2; member[5][1] = 4;
    member[6][0] = 2; member[6][1] = 5;
    member[7][0] = 1; member[7][1] = 3;
    member[8][0] = 1; member[8][1] = 6;
    member[9][0] = 3; member[9][1] = 6;
    member[10][0] = 4; member[10][1] = 5;
    member[11][0] = 3; member[11][1] = 4;
    member[12][0] = 5; member[12][1] = 6;
    member[13][0] = 3; member[13][1] = 10;
    member[14][0] = 6; member[14][1] = 7;
    member[15][0] = 4; member[15][1] = 9;
    member[16][0] = 5; member[16][1] = 8;
    member[17][0] = 4; member[17][1] = 7;
    member[18][0] = 3; member[18][1] = 8;
    member[19][0] = 5; member[19][1] = 10;
    member[20][0] = 6; member[20][1] = 9;
    member[21][0] = 6; member[21][1] = 10;
    member[22][0] = 3; member[22][1] = 7;
    member[23][0] = 4; member[23][1] = 8;
    member[24][0] = 5; member[24][1] = 9;

	for(int i=0;i<num_elem;i++){
		ymod[i] = E;
	}
    
    //compute member lengths
    for(int i=0;i<num_elem;i++){
        int node1, node2;
        double dx, dy, dz;
        
        node1 = member[i][0]-1;
        node2 = member[i][1]-1;
        
        dx = coordx[node1]-coordx[node2];
        dy = coordy[node1]-coordy[node2];
        dz = coordz[node1]-coordz[node2];
        
        length[i] = sqrt( pow(dx, 2) + pow(dy, 2) + pow(dz, 2) );
    }

	//min and max areas
	minArea = 0.01;
	maxArea = 3.4; 

	//compute max and min weights
	maxW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		maxW += ro*maxArea*length[i];
	}

	minW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		minW += ro*minArea*length[i];
	}
}

void setupTruss72(int nThreads, int *avalPerInd){
	num_elem = 72;                                           //GLOBAL
    num_node = 20;                                            //GLOBAL
	dim = 16;
    
    E = 1.0e4; //ksi

	*avalPerInd = 2;

    ro = 0.1;  //lb/in^3

	maxDisp = 0.25;

	maxStress = (double**)malloc(num_elem*sizeof(double*));
    minStress = (double**)malloc(num_elem*sizeof(double*));

	double *maxStressData = (double*)malloc(num_elem*nThreads*sizeof(double));
	double *minStressData = (double*)malloc(num_elem*nThreads*sizeof(double));

	for(int i=0;i<num_elem;i++){
		maxStress[i] = maxStressData + i*nThreads;
		minStress[i] = minStressData + i*nThreads;
	}

	for(int j=0;j<nThreads;j++){
		for(int i=0;i<num_elem;i++){
			maxStress[i][j] = 25.0;  //ksi
			minStress[i][j] = -25.0; //ksi
		}
	}

	distributeAreas = distributeAreas72;
    fitChoice = fit72;
    
    kode = (int**)malloc(num_node*sizeof(int*));              //GLOBAL
    int *kodeData = (int*)calloc(num_node*3, sizeof(int));
    for(int i=0;i<num_node;i++){
        kode[i] = kodeData + i*3;
    }
    
    coordx = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    coordy = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    coordz = (double*)malloc(num_node*sizeof(double));

	fX = (double**)malloc(num_node*sizeof(double*));
    fY = (double**)malloc(num_node*sizeof(double*));
	fZ = (double**)malloc(num_node*sizeof(double*));

	double *fXdata = (double*)calloc(num_node*nThreads, sizeof(double));
	double *fYdata = (double*)calloc(num_node*nThreads, sizeof(double));
	double *fZdata = (double*)calloc(num_node*nThreads, sizeof(double));
	for(int i=0;i<num_node;i++){
		fX[i] = fXdata + i*nThreads;
		fY[i] = fYdata + i*nThreads;
		fZ[i] = fZdata + i*nThreads;
	}
    
    int *memberData;
    member = (int**)malloc(num_elem*sizeof(int*));               //GLOBAL
    memberData = (int*)malloc(num_elem*2*sizeof(int));
    for(int i=0;i<num_elem;i++){
        member[i] = memberData + i*2;
    }
    
    length = (double*)malloc(num_elem*sizeof(double));   //GLOBAL
    
    ymod = (double*)malloc(num_elem*sizeof(double));
    
    /* fixo X e Y*/
	for(int i=0;i<5;i++){
		coordx[i*4] = 0.0;
		coordy[i*4] = 0.0;
		coordz[i*4] = 60.0*i;

		coordx[i*4+1] = 120.0;
		coordy[i*4+1] = 0.0;
		coordz[i*4+1] = 60.0*i;

		coordx[i*4+2] = 120.0;
		coordy[i*4+2] = 120.0;
		coordz[i*4+2] = 60.0*i;

		coordx[i*4+3] = 0.0;
		coordy[i*4+3] = 120.0;
		coordz[i*4+3] = 60.0*i;
	}
    
    kode[0][0] = 1;
    kode[0][1] = 1;
    kode[0][2] = 1;
    kode[1][0] = 1;
    kode[1][1] = 1;
    kode[1][2] = 1;
    kode[2][0] = 1;
    kode[2][1] = 1;
    kode[2][2] = 1;
    kode[3][0] = 1;
    kode[3][1] = 1;
    kode[3][2] = 1;

	for(int i=0;i<4;i++){
		member[i*18][0] = i*4 + 1; member[i*18][1] = i*4 + 5;
		member[i*18+1][0] = i*4 + 2; member[i*18+1][1] = i*4 + 6;
		member[i*18+2][0] = i*4 + 3; member[i*18+2][1] = i*4 + 7;
		member[i*18+3][0] = i*4 + 4; member[i*18+3][1] = i*4 + 8;

		member[i*18+4][0] = i*4 + 1; member[i*18+4][1] = i*4 + 6;
		member[i*18+5][0] = i*4 + 2; member[i*18+5][1] = i*4 + 5;
		member[i*18+6][0] = i*4 + 2; member[i*18+6][1] = i*4 + 7;
		member[i*18+7][0] = i*4 + 3; member[i*18+7][1] = i*4 + 6;
		member[i*18+8][0] = i*4 + 3; member[i*18+8][1] = i*4 + 8;
		member[i*18+9][0] = i*4 + 4; member[i*18+9][1] = i*4 + 7;
		member[i*18+10][0] = i*4 + 1; member[i*18+10][1] = i*4 + 8;
		member[i*18+11][0] = i*4 + 4; member[i*18+11][1] = i*4 + 5;

		member[i*18+12][0] = i*4 + 5; member[i*18+12][1] = i*4 + 6;
		member[i*18+13][0] = i*4 + 6; member[i*18+13][1] = i*4 + 7;
		member[i*18+14][0] = i*4 + 7; member[i*18+14][1] = i*4 + 8;
		member[i*18+15][0] = i*4 + 8; member[i*18+15][1] = i*4 + 5;

		member[i*18+16][0] = i*4 + 5; member[i*18+16][1] = i*4 + 7;
		member[i*18+17][0] = i*4 + 6; member[i*18+17][1] = i*4 + 8;
	}
	
	for(int i=0;i<num_elem;i++){
		ymod[i] = E;
	}
    
    //compute member lengths
    for(int i=0;i<num_elem;i++){
        int node1, node2;
        double dx, dy, dz;
        
        node1 = member[i][0]-1;
        node2 = member[i][1]-1;
        
        dx = coordx[node1]-coordx[node2];
        dy = coordy[node1]-coordy[node2];
        dz = coordz[node1]-coordz[node2];
        
        length[i] = sqrt( pow(dx, 2) + pow(dy, 2) + pow(dz, 2) );
    }

	maxArea = 20.0; 

	//compute max and min weights
	maxW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		maxW += ro*maxArea*length[i];
	}

	minW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		minW += ro*minArea*length[i];
	}
}

void setupTruss72case1(int nThreads, int *avalPerInd){

	//load case 1
	int tempAval;

	setupTruss72(nThreads, &tempAval);
	*avalPerInd = tempAval;

	minArea = 0.1;
}

void setupTruss72case2(int nThreads, int *avalPerInd){
	
	//load case 2
	int tempAval;

	setupTruss72(nThreads, &tempAval);
	*avalPerInd = tempAval;

	minArea = 0.01;
}

void setupTruss120(int nThreads, int *avalPerInd){
    num_elem = 120;                                           //GLOBAL
    num_node = 49;                                            //GLOBAL
	dim = 7;
    
    E = 30450.0; //ksi

	*avalPerInd = 1;

    ro = 0.288;  //lb/in^3

	maxDisp = 0.1969;  //5mm, 0.1969in
	double yield = 58.0;
	
    maxStress = (double**)malloc(num_elem*sizeof(double*));
    minStress = (double**)malloc(num_elem*sizeof(double*));

	double *maxStressData = (double*)malloc(num_elem*nThreads*sizeof(double));
	double *minStressData = (double*)malloc(num_elem*nThreads*sizeof(double));

	for(int i=0;i<num_elem;i++){
		maxStress[i] = maxStressData + i*nThreads;
		minStress[i] = minStressData + i*nThreads;
	}

	for(int j=0;j<nThreads;j++){
		for(int i=0;i<num_elem;i++){
			maxStress[i][j] = yield*0.6;  //ksi
		}
	}

	distributeAreas = distributeAreas120;
    fitChoice = fit120;
    
    kode = (int**)malloc(num_node*sizeof(int*));              //GLOBAL
    int *kodeData = (int*)calloc(num_node*3, sizeof(int));
    for(int i=0;i<num_node;i++){
        kode[i] = kodeData + i*3;
    }
    
    coordx = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    coordy = (double*)malloc(num_node*sizeof(double));   //GLOBAL
    coordz = (double*)malloc(num_node*sizeof(double));

	fX = (double**)malloc(num_node*sizeof(double*));
    fY = (double**)malloc(num_node*sizeof(double*));
	fZ = (double**)malloc(num_node*sizeof(double*));

	double *fXdata = (double*)calloc(num_node*nThreads, sizeof(double));
	double *fYdata = (double*)calloc(num_node*nThreads, sizeof(double));
	double *fZdata = (double*)calloc(num_node*nThreads, sizeof(double));
	for(int i=0;i<num_node;i++){
		fX[i] = fXdata + i*nThreads;
		fY[i] = fYdata + i*nThreads;
		fZ[i] = fZdata + i*nThreads;
	}

	double load1 = 13.49; //kips
    double load2 = 6.744; //kips
	double load3 = 2.248; //kips

	for(int i=0;i<nThreads;i++){
		fZ[0][i] = -load1;
		for(int j=1;j<14;j++){
			fZ[j][i] = -load2;
		}
		for(int j=14;j<37;j++){
			fZ[j][i] = -load3;
		}
	}
    
    int *memberData;
    member = (int**)malloc(num_elem*sizeof(int*));               //GLOBAL
    memberData = (int*)malloc(num_elem*2*sizeof(int));
    for(int i=0;i<num_elem;i++){
        member[i] = memberData + i*2;
    }
    
    length = (double*)malloc(num_elem*sizeof(double));   //GLOBAL
    
    ymod = (double*)malloc(num_elem*sizeof(double));
    
    /* fixo X e Y*/
    coordx[0] = 0.0;
	coordy[0] = 0.0;
    coordz[0] = 275.59;

	for(int i=0;i<12;i++){
		coordx[i+1] = 273.26*cos((i*30.0e0)*PI/180.0e0);
		coordy[i+1] = 273.26*sin((i*30.0e0)*PI/180.0e0);
		coordz[i+1] = 196.85;

		coordx[i+37] = 625.59*cos((i*30.0e0)*PI/180.0e0);
		coordy[i+37] = 625.59*sin((i*30.0e0)*PI/180.0e0);
		coordz[i+37] = 0.0;
	}
	for(int i=0;i<24;i++){
		coordx[i+13] = 492.12*cos((i*15.0e0)*PI/180.0e0);
		coordy[i+13] = 492.12*sin((i*15.0e0)*PI/180.0e0);
		coordz[i+13] = 118.11;
	}
    
	for(int i=0;i<12;i++){
		kode[i+37][0] = 1;
		kode[i+37][1] = 1;
		kode[i+37][2] = 1;
	}

	for(int i=0;i<12;i++){
		member[i][0] = 1; member[i][1] = i+2;
		
		member[i+12][0] = i+2; member[i+12][1] = i+3;

		member[i+24][0] = i+2; member[i+24][1] = 14+2*i;

		member[i+36][0] = i+2; member[i+36][1] = 15+2*i;
		member[i+48][0] = i+3; member[i+48][1] = 15+2*i;

		member[i+60][0] = i+14; member[i+60][1] = i+15;
		member[i+72][0] = i+26; member[i+72][1] = i+27;

		member[i+84][0] = 14+2*i; member[i+84][1] = 38+i;

		member[i+96][0] = 15+2*i; member[i+96][1] = 38+i;
		member[i+108][0] = 15+2*i; member[i+108][1] = 39+i;
	}
	member[23][0] = 2; member[23][1] = 13;
	member[59][0] = 2; member[59][1] = 37;
	member[83][0] = 14; member[83][1] = 37;
	member[119][0] = 37; member[119][1] = 38;

	for(int i=0;i<num_elem;i++){
		ymod[i] = E;
	}
    
    //compute member lengths
    for(int i=0;i<num_elem;i++){
        int node1, node2;
        double dx, dy, dz;
        
        node1 = member[i][0]-1;
        node2 = member[i][1]-1;
        
        dx = coordx[node1]-coordx[node2];
        dy = coordy[node1]-coordy[node2];
        dz = coordz[node1]-coordz[node2];
        
        length[i] = sqrt( pow(dx, 2) + pow(dy, 2) + pow(dz, 2) );
    }

	//min and max areas
	minArea = 0.775;
	maxArea = 20.0; 

	//compute max and min weights
	maxW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		maxW += ro*maxArea*length[i];
	}

	minW = 0.0;
	for (int i = 0; i < num_elem; i++) {
		minW += ro*minArea*length[i];
	}
}

void distributeAreasEqual(double *sol, double *area, int id){
	for(int i=0;i<num_elem;i++){
		area[i] = sol[i];
	}
}

void distributeAreas200(double *sol, double *area, int id){
	int i, j;
	
	for(j=0;j<4;j++){
		area[j] = sol[0];
	}
	for(j=0;j<5;j++){
		area[4 + 3*j] = sol[1];
	}
	for(j=0;j<6;j++){
		area[18 + j] = sol[2];
	}
	for(j=0;j<5;j++){
		area[17 + 38*j] = sol[3];
		area[24 + 38*j] = sol[3];
	}
	//0 1 2 3

	for(i=0;i<5;i++){
		for(j=0;j<5;j++){
			area[i*38 + 3*j + 25] = sol[4 + 5*i];
		}
	}
	//4 9 14 19 24

	for(i=0;i<5;i++){
		for(j=0;j<4;j++){
			area[38*i + 3*j + 5] = sol[5 + 5*i];
			area[38*i + 3*j + 6] = sol[5 + 5*i];

			area[38*i + 3*j + 26] = sol[5 + 5*i];
			area[38*i + 3*j + 27] = sol[5 + 5*i];
		}
	}
	//5 10 15 20 25

	for(i=0;i<5;i++){
		for(j=0;j<4;j++){
			area[i*38 + j + 38] = sol[6 + 5*i];
		}
	}
	//6 11 16 21 26

	for(i=0;i<4;i++){
		for(j=0;j<5;j++){
			area[i*38 + 3*j + 42] = sol[7 + 5*i];
		}
	}
	//7 12 17 22

	for(i=0;i<4;i++){
		for(j=0;j<6;j++){
			area[i*38 + j + 56] = sol[8 + 5*i];
		}
	}
	//8 13 18 23

	area[194] = sol[27];
	area[196] = sol[27];
	area[197] = sol[27];
	area[199] = sol[27];

	area[195] = sol[28];
	area[198] = sol[28];

	/*
	for(i=0;i<num_elem;i++){
		//printf("%d: %lf\n", i, area[i]);
		if(area[i] == 0.0e0){
			printf("%d: %lf\n", i, area[i]);
		}
	}
	*/
}

void distributeAreas25(double *sol, double *area, int id){
	int i;

	area[0] = sol[0];
	maxStress[0][id] = 40.0;
	minStress[0][id] = -35.092;
	for(i=0;i<4;i++){
		area[i+1] = sol[1];
		maxStress[i+1][id] = 40.0;
		minStress[i+1][id] = -11.590;
	}
	for(i=0;i<4;i++){
		area[i+5] = sol[2];
		maxStress[i+5][id] = 40.0;
		minStress[i+5][id] = -17.305;
	}
	for(i=0;i<2;i++){
		area[i+9] = sol[3];
		maxStress[i+9][id] = 40.0;
		minStress[i+9][id] = -35.092;
	}
	for(i=0;i<2;i++){
		area[i+11] = sol[4];
		maxStress[i+11][id] = 40.0;
		minStress[i+11][id] = -35.092;
	}
	for(i=0;i<4;i++){
		area[i+13] = sol[5];
		maxStress[i+13][id] = 40.0;
		minStress[i+13][id] = -6.759;
	}
	for(i=0;i<4;i++){
		area[i+17] = sol[6];
		maxStress[i+17][id] = 40.0;
		minStress[i+17][id] = -6.959;
	}
	for(i=0;i<4;i++){
		area[i+21] = sol[7];
		maxStress[i+21][id] = 40.0;
		minStress[i+21][id] = -11.802;
	}

}

void distributeAreas72(double *sol, double *area, int id){
	int i, j;

	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			area[18*i + j] = sol[4*i];
		}
		for(j=4;j<12;j++){
			area[18*i + j] = sol[4*i + 1];
		}
		for(j=12;j<16;j++){
			area[18*i + j] = sol[4*i + 2];
		}
		for(j=16;j<18;j++){
			area[18*i + j] = sol[4*i + 3];
		}
	}
	/*
	for(i=0;i<num_elem;i++){
		printf("%d: %lf\n", i+1, area[i]);
	}
	*/
}

void distributeAreas120(double *sol, double *area, int id){
	int i;
	
	for(i=0;i<12;i++){
		area[i] = sol[0];
	}
	for(i=12;i<24;i++){
		area[i] = sol[1];
	}
	for(i=24;i<36;i++){
		area[i] = sol[2];
	}
	for(i=36;i<60;i++){
		area[i] = sol[3];
	}
	for(i=60;i<84;i++){
		area[i] = sol[4];
	}
	for(i=84;i<96;i++){
		area[i] = sol[5];
	}
	for(i=96;i<num_elem;i++){
		area[i] = sol[6];
	}
	/*
	groups:
	1: 0 - 11
	2: 12 - 23
	3: 24 - 35
	4: 36 - 59
	5: 60 - 83
	6: 84 - 95
	7: 96 - 119
	*/
}

void destroyTruss(){
    free(kode[0]);
    free(kode);
	free(maxStress[0]);
	free(maxStress);
	free(minStress[0]);
	free(minStress);
    free(coordx);
    free(coordy);
    if(coordz != NULL){
        free(coordz);
    }
	
	free(fX[0]);
    free(fX);
	free(fY[0]);
    free(fY);
    if(fZ != NULL){
        free(fZ[0]);
		free(fZ);
    }
	
    free(ymod);
    free(member[0]);
    free(member);
    free(length);
}

double trussFit(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id){
    return fitChoice(sol, weightRet, pen, verb, feasible, id);
}

double trussFit2D(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id){
   
    int i, j;

	double *area = (double*)calloc(num_elem, sizeof(double));
	(*distributeAreas)(sol, area, id);
    
    double *stress = (double*)malloc(num_elem*sizeof(double));
	double *displX = (double*)malloc(num_node*sizeof(double));
	double *displY = (double*)malloc(num_node*sizeof(double));
	
	double weight = 0.0;
	double normW = 0.0;
	double normSumPDisp = 0.0;
	double normSumPStressTens = 0.0;
	double normSumPStressComp = 0.0;
	double normTotalPenalty = 0.0;

	double penCoef = 0.5;
	double penConst = 5.0e-4;

    int lm[4] = {0,0,0,0}, 
        neq = 0;

    double *b = (double*)malloc(num_node*2*sizeof(double));
    
    double **a;
    double *aData;
    a = (double**)malloc(num_node*2*sizeof(double*));
    aData = (double*)malloc(num_node*2*num_node*2*sizeof(double));
    for(i=0;i<num_node*2;i++){
        a[i] = aData + i*num_node*2;
    }
    
    double stiffmatrix[4][4];

    int n = 0, nn = 0;

    double dx = 0.,dy = 0.;
    double du = 0.,dv = 0.,p = 0.,dl = 0.;

    int nl = 0;
    int n1 = 0,m1 = 0;
    int neq1 = 0;

    int m = 0, k = 0, g = 0;
    double cosa = 0.,sina = 0.,comm = 0.;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		       stiffmatrix[i][j] = 0.;

	//stiff
	/* Initialize Stiffness Matrix and Load Vector */
	neq = 2 * num_node;
	for (m=0; m<neq;m++)
	{
		b[m] = 0.;
		for(k=0;k<neq;k++) a[m][k]=0.;
	}
	/*  >>>>>>> loop over all elements  <<<<<<< */
	for(m=0;m<num_elem;m++)
	{
		/* Determination of Length, Direction */
        
		i=member[m][0] - 1;
		j=member[m][1] - 1;
		dx = coordx[j]-coordx[i];
		dy = coordy[j]-coordy[i];
        
		cosa=dx/length[m];
		sina=dy/length[m];
		comm=area[m]*ymod[m]/length[m];
		/* Construction of the 4X4 Element Stiffness Matrix */
		stiffmatrix[0][0] = cosa * cosa * comm;
		stiffmatrix[0][1] = cosa * sina * comm;
		stiffmatrix[0][2] = -stiffmatrix[0][0];
		stiffmatrix[0][3] = -stiffmatrix[0][1];
		stiffmatrix[1][0] =  stiffmatrix[0][1];
		stiffmatrix[1][1] = sina * sina * comm;
		stiffmatrix[1][2] = -stiffmatrix[0][1];
		stiffmatrix[1][3] = -stiffmatrix[1][1];
		stiffmatrix[2][0] =  stiffmatrix[0][2];
		stiffmatrix[2][1] =  stiffmatrix[1][2];
		stiffmatrix[2][2] =  stiffmatrix[0][0];
		stiffmatrix[2][3] =  stiffmatrix[0][1];
		stiffmatrix[3][0] =  stiffmatrix[0][3];
		stiffmatrix[3][1] =  stiffmatrix[1][3];
		stiffmatrix[3][2] =  stiffmatrix[2][3];
		stiffmatrix[3][3] =  stiffmatrix[1][1];

		/* Assembly of Element Stiffness to Overall
		   Stiffness   */
		lm[1] = 2 * member[m][0] - 1;
		lm[0] = lm[1] - 1;
		lm[3] = 2 * member[m][1] - 1;
		lm[2] = lm[3] - 1;
		for (k = 0; k < 4; k++)
		{
			for(g=0; g<4; g++)
			{
				i = lm [k];
				j = lm [g];
				a[i][j] = a[i][j] + stiffmatrix[k][g];
			}
		}

	}
	/*end stiff*/

	/*displ*/
	for(n = 0; n < num_node; n++)
	{
		nn = (2 * n) + 1;
		b[nn] = fY[n][id];
		b[nn-1] = fX[n][id];
	}

	for(n = 0; n < num_node; n++)
	{
		nn = (2 * n) + 1;
        
        if(kode[n][0] == 1){
            for(i=0;i<neq;i++){
                b[i]=b[i]-(a[i][nn-1]*fX[n][id]);
                a[i][nn-1]=0.;
                a[nn-1][i]=0.;
            }
            a[nn-1][nn-1]=1.00;
            b[nn-1]=fX[n][id];
        }
        
        if(kode[n][1] == 1){
            for(i=0;i<neq;i++){
                b[i]=b[i]-(a[i][nn]*fY[n][id]);
                a[i][nn]=0.;
                a[nn][i]=0.;
            }
            a[nn][nn]=1.00;
            b[nn]=fY[n][id];
        }
        
	};


	/*symsol();*/
	neq1 = nl = 0;
	neq1=neq-1;
	nl = neq - 1;

	for (n=0;n<nl;n++)
	{
		if (a[n][n]<=0.)
		{
			printf("zero or negative main-diagonal %d\n",n);
			exit(2);
		}

		n1 = n +1 ;
		for (j=n1;j<neq;j++) {
			a[n][j]=a[n][j]/a[n][n];
		}
		for(i=n1;i<neq;i++)
		{
			if(a[n][i]==0.0)
				b[i]=b[i]-a[n][i]*b[n];
			else
			{
				for(j=i;j<neq;j++)
				{
					a[i][j]=a[i][j]-a[i][n]*a[n][j];
					a[j][i]=a[i][j];
				}
			}
			b[i]=b[i]-a[n][i]*b[n];
		}
		b[n]=b[n]/a[n][n];
	};

	m=neq1;
	b[m]=b[m]/a[m][m];
	for(n=0;n<nl;n++)
	{
		m1=m;
		m=m-1;
		for(j=m1;j<neq;j++) {
			b[m]=b[m]-b[j]*a[m][j];
		}
	}
	/*end symbol*/

	for(i=0;i<num_node;i++)
	{
		displX[i] = b[2*i];
		displY[i] = b[2*i+1];
	}
	/*end displ*/


	/*stress*/
	/*printf("\nELEM I-NOD J-NOD    MEM-FOR     STRESS\n");*/
	for(m=0;m<num_elem;m++)
	{
		i   =  member[m][0] - 1;
		j   =  member[m][1] - 1;
		dx  =  coordx[j]-coordx[i];
		dy  =  coordy[j]-coordy[i];
		du  =  displX[j] -displX[i];
		dv  =  displY[j] -displY[i];
		dl  =  dv*dy/length[m]+du*dx/length[m];
		stress[m] =  dl*ymod[m]/length[m];
		p   =  area[m]*stress[m];
	}
	/*end stress*/

	//truss weight
	for (i=0; i<num_elem;i++) {
		weight += ro*area[i]*length[i];
	}
	//normalize weight
	normW = (weight - minW)/(maxW-minW);

	//displacement penalty
	for(i=0;i<num_node;i++){
		double diff = fabs(displX[i]) - maxDisp;
		//printf("%lf\n", fabs(displX[i]));
		if(diff > 0.0){
			normSumPDisp += diff;
		}

		diff = fabs(displY[i]) - maxDisp;
		//printf("%lf\n", fabs(displY[i]));
		if(diff > 0.0){
			normSumPDisp += diff;
		}
	}
	//normalization
	normSumPDisp = normSumPDisp/maxDisp;

	//stress penalty
	for(i=0;i<num_elem;i++){
		double diff;
		//tension stress
		if(stress[i]>=0){
			diff = stress[i] - maxStress[i][id];
			if(diff>0){
				normSumPStressTens += diff/maxStress[i][id];
			}
		}
		//compression stress
		else{
			diff = stress[i] - minStress[i][id];
			if(diff<0){
				normSumPStressComp += diff/minStress[i][id];
			}
		}
	}
	
	//normalization
	normTotalPenalty = normSumPDisp + normSumPStressTens + normSumPStressComp;
	if(normTotalPenalty > 0.0e0){
		normTotalPenalty += penConst;
	}
    
	if(verb){
		printf("Weight:\n");
		printf("%lf\n", weight);
		
		printf("Displacements:\n");
		for(i=0;i<num_node;i++){
			printf("%d: %lf, %lf\n", i, displX[i], displY[i]);
		}
		
		printf("Stresses:\n");
		for(i=0;i<num_elem;i++){
			printf("%d: %lf\n", i, stress[i]);
		}

		if(normTotalPenalty > 0.0){
			printf("Infeasible solution\n");
			*feasible = false;
			printf("Normalized penalty: %lf\n", normTotalPenalty);
		}
		else{
			printf("Feasible solution\n");
			*feasible = true;
		}

		//printf("\nDIsp pen: %lf\n", normSumPDisp);
	}
	else{
		if(normTotalPenalty > 0.0){
			*feasible = false;
		}
		else{
			*feasible = true;
		}
	}
    
	free(area);
    free(stress);
    free(displX);
    free(displY);
    free(b);
    free(a);
    free(aData);

	//printf("maxW: %20.15e\n", maxW);
	//printf("minW: %20.15e\n", minW);
    
    *weightRet = weight;
	*pen = normTotalPenalty;
	return (normW + penCoef*normTotalPenalty);
}

double trussFit3D(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id){
    
    int i, j;

	double *area = (double*)calloc(num_elem, sizeof(double));
	(*distributeAreas)(sol, area, id);
    
    double *stress = (double*)malloc(num_elem*sizeof(double));
	double *displX = (double*)malloc(num_node*sizeof(double));
	double *displY = (double*)malloc(num_node*sizeof(double));
    double *displZ = (double*)malloc(num_node*sizeof(double));
    
	double weight = 0.0;
	double normW = 0.0;
	double normSumPDisp = 0.0;
	double normSumPStressTens = 0.0;
	double normSumPStressComp = 0.0;
	double normTotalPenalty = 0.0;
    
	double penCoef = 0.5;
	double penConst = 5.0e-4;

    int lm[6] = {0,0,0,0,0,0}, 
        neq = 0;

    double *b = (double*)malloc(num_node*3*sizeof(double));
    
    
    double **a;
    double *aData;
    a = (double**)malloc(num_node*3*sizeof(double*));
    aData = (double*)malloc(num_node*3*num_node*3*sizeof(double));
    for(i=0;i<num_node*3;i++){
        a[i] = aData + i*num_node*3;
    }
    
    double stiffmatrix[6][6];

    int n = 0, nn = 0;

    double dx = 0.,dy = 0., dz = 0.;
    double du = 0.,dv = 0., dw = 0., p = 0.,dl = 0.;

    int nl = 0;
    int n1 = 0,m1 = 0;
    int neq1 = 0;

    int m = 0, k = 0, g = 0;
    double cosx = 0.,cosy = 0., cosz = 0.,comm = 0.;

	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
		       stiffmatrix[i][j] = 0.;

    //////////////////
    //stiff */
	/* Initialize Stiffness Matrix and Load Vector */
	neq = 3 * num_node;
	for (m=0; m<neq;m++)
	{
		b[m] = 0.;
		for(k=0;k<neq;k++) a[m][k]=0.;
	}
	/*  >>>>>>> loop over all elements  <<<<<<< */
	for(m=0;m<num_elem;m++)
	{
		/* Determination of Length, Direction */
        
		i=member[m][0] - 1;
		j=member[m][1] - 1;
		dx = coordx[j]-coordx[i];
		dy = coordy[j]-coordy[i];
        dz = coordz[j]-coordz[i];
        
		cosx=dx/length[m];
        cosy=dy/length[m];
        cosz=dz/length[m];
		comm=area[m]*ymod[m]/length[m];
        
        /* Construction of the 6X6 Element Stiffness Matrix */
		stiffmatrix[0][0] = cosx * cosx * comm;
		stiffmatrix[0][1] = cosx * cosy * comm;
		stiffmatrix[0][2] = cosx * cosz * comm;
        
		stiffmatrix[1][0] = stiffmatrix[0][1];
		stiffmatrix[1][1] = cosy * cosy * comm;
		stiffmatrix[1][2] = cosy * cosz * comm;
        
		stiffmatrix[2][0] =  cosx * cosz * comm;
		stiffmatrix[2][1] =  stiffmatrix[1][2];
		stiffmatrix[2][2] =  cosz * cosz * comm;
        
        for(k=0;k<3;k++){
            for(g=0;g<3;g++){
                stiffmatrix[k][g+3] = -stiffmatrix[k][g];
                stiffmatrix[k+3][g] = -stiffmatrix[k][g];
                stiffmatrix[k+3][g+3] = stiffmatrix[k][g];
            }
        }

		/* Assembly of Element Stiffness to Overall
		   Stiffness   */
		lm[2] = 3 * member[m][0] - 1;
		lm[1] = lm[2] - 1;
        lm[0] = lm[2] - 2;
		lm[5] = 3 * member[m][1] - 1;
		lm[4] = lm[5] - 1;
        lm[3] = lm[5] - 2;
        
		for (k = 0; k < 6; k++){
			for(g=0; g<6; g++){
				i = lm [k];
				j = lm [g];
				a[i][j] = a[i][j] + stiffmatrix[k][g];
			}
		}

	}
    /*end stiff*/
    
    /*displ*/
	for(n = 0; n < num_node; n++)
	{
		nn = (3 * n) + 2;
		b[nn] = fZ[n][id];
		b[nn-1] = fY[n][id];
        b[nn-2] = fX[n][id];
	}

	for(n = 0; n < num_node; n++)
	{
		nn = (3 * n) + 2;
        
        if(kode[n][0] == 1){
            for(i=0;i<neq;i++){
                b[i]=b[i]-(a[i][nn-2]*fX[n][id]);
                a[i][nn-2]=0.;
                a[nn-2][i]=0.;
            }
            a[nn-2][nn-2]=1.00;
            b[nn-2]=fX[n][id];
        }
        if(kode[n][1] == 1){
            for(i=0;i<neq;i++){
                b[i]=b[i]-(a[i][nn-1]*fY[n][id]);
                a[i][nn-1]=0.;
                a[nn-1][i]=0.;
            }
            a[nn-1][nn-1]=1.00;
            b[nn-1]=fY[n][id];
        }
        if(kode[n][2] == 1){
            for(i=0;i<neq;i++){
                b[i]=b[i]-(a[i][nn]*fZ[n][id]);
                a[i][nn]=0.;
                a[nn][i]=0.;
            }
            a[nn][nn]=1.00;
            b[nn]=fZ[n][id];
        }
        
	};

    /*symsol();*/
	neq1 = nl = 0;
	neq1=neq-1;
	nl = neq - 1;

	for (n=0;n<nl;n++)
	{
		if (a[n][n]<=0.)
		{
			printf("zero or negative main-diagonal %d\n",n);
			exit(2);
		}

		n1 = n +1 ;
		for (j=n1;j<neq;j++) {
			a[n][j]=a[n][j]/a[n][n];
		}
		for(i=n1;i<neq;i++)
		{
			if(a[n][i]==0.0)
				b[i]=b[i]-a[n][i]*b[n];
			else
			{
				for(j=i;j<neq;j++)
				{
					a[i][j]=a[i][j]-a[i][n]*a[n][j];
					a[j][i]=a[i][j];
				}
			}
			b[i]=b[i]-a[n][i]*b[n];
		}
		b[n]=b[n]/a[n][n];
	};

	m=neq1;
	b[m]=b[m]/a[m][m];
	for(n=0;n<nl;n++)
	{
		m1=m;
		m=m-1;
		for(j=m1;j<neq;j++) {
			b[m]=b[m]-b[j]*a[m][j];
		}
	}
    /*end symsol*/

	for(i=0;i<num_node;i++)
	{
		displX[i] = b[3*i];
		displY[i] = b[3*i+1];
        displZ[i] = b[3*i+2];
	}

    
    /*stress*/
	/*printf("\nELEM I-NOD J-NOD    MEM-FOR     STRESS\n");*/
	for(m=0;m<num_elem;m++)
	{
		i   =  member[m][0] - 1;
		j   =  member[m][1] - 1;
		dx  =  coordx[j]-coordx[i];
		dy  =  coordy[j]-coordy[i];
        dz  =  coordz[j]-coordz[i];
		du  =  displX[j] - displX[i];
		dv  =  displY[j] - displY[i];
        dw  =  displZ[j] - displZ[i];
        dl = (du*dx + dv*dy + dw*dz)/length[m];
		stress[m] =  dl*ymod[m]/length[m];
		p   =  area[m]*stress[m];

		/*printf("%d    %d    %d  %10.4f  %10.4f\n", (m+1), (i+1), (j+1), p,strss[m]);*/
	}
    /*end stress*/

	//truss weight
	for (i=0; i<num_elem;i++) {
		weight += ro*area[i]*length[i];
	}
	//normalize weight
	normW = (weight - minW)/(maxW-minW);

	//displacement penalty
	for(i=0;i<num_node;i++){
		double diff = fabs(displX[i]) - maxDisp;
		//printf("%lf\n", fabs(displX[i]));
		if(diff > 0.0){
			normSumPDisp += diff;
		}

		diff = fabs(displY[i]) - maxDisp;
		//printf("%lf\n", fabs(displY[i]));
		if(diff > 0.0){
			normSumPDisp += diff;
		}
        
        diff = fabs(displZ[i]) - maxDisp;
		//printf("%lf\n", fabs(displY[i]));
		if(diff > 0.0){
			normSumPDisp += diff;
		}
	}
	//normalization
	normSumPDisp = normSumPDisp/maxDisp;

	//stress penalty
	for(i=0;i<num_elem;i++){
		double diff;
		//tension stress
		if(stress[i]>=0){
			diff = stress[i] - maxStress[i][id];
			if(diff>0){
				normSumPStressTens += diff/maxStress[i][id];
			}
		}
		//compression stress
		else{
			diff = stress[i] - minStress[i][id];
			if(diff<0){
				normSumPStressComp += diff/minStress[i][id];
			}
		}
	}
	//normalization	
	normTotalPenalty = normSumPDisp + normSumPStressTens + normSumPStressComp;
	if(normTotalPenalty > 0.0e0){
		normTotalPenalty += penConst;
	}
    
	if(verb){
		printf("Weight:\n");
		printf("%lf\n", weight);
		
		printf("Displacements:\n");
		for(i=0;i<num_node;i++){
			printf("%d: %lf, %lf, %lf\n", i, displX[i], displY[i], displZ[i]);
		}
		
		printf("Stresses:\n");
		for(i=0;i<num_elem;i++){
			printf("%d: %lf\n", i, stress[i]);
		}

		if(normTotalPenalty > 0.0){
			printf("Infeasible solution\n");
			*feasible = false;
			printf("Normalized penalty: %lf\n", normTotalPenalty);
		}
		else{
			printf("Feasible solution\n");
			*feasible = true;
		}

		//printf("\nDIsp pen: %lf\n", normSumPDisp);
	}
	else{
		if(normTotalPenalty > 0.0){
			*feasible = false;
		}
		else{
			*feasible = true;
		}
	}
    
    free(area);
    free(stress);
    free(displX);
    free(displY);
    free(displZ);
    free(b);
    free(a);
    free(aData);
    
    *weightRet = weight;
	*pen = normTotalPenalty;
	return (normW + penCoef*normTotalPenalty);
}

double fit25(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id){
	double f1, f2;
	double w1, w2;
	double p1, p2;
	bool feasible1, feasible2;

	//change to case 2 after first evaluation
	double load1;
	double load2;
	double load3;
	double load4;

	//reset load case 2
	fX[0][id] = 0.0;
	fY[0][id] = 0.0;
	fZ[0][id] = 0.0;
	fY[1][id] = 0.0;
	fZ[1][id] = 0.0;
	fX[2][id] = 0.0;
	fX[5][id] = 0.0;

	//load case 1
	load1 = 20.0; //kips
    load2 = 5.0;  //kips

	fY[0][id] = load1;
    fZ[0][id] = -load2;
    fY[1][id] = -load1;
    fZ[1][id] = -load2;

	f1 = trussFit3D(sol, &w1, &p1, verb, &feasible1, id);

	//reseting load case 1
	fY[0][id] = 0.0;
    fZ[0][id] = 0.0;
    fY[1][id] = 0.0;
    fZ[1][id] = 0.0;

	load1 = 1.0;
	load2 = 10.0;
	load3 = 5.0;
	load4 = 0.5;

	//load case 2
	fX[0][id] = load1;
	fY[0][id] = load2;
	fZ[0][id] = -load3;
	fY[1][id] = load2;
	fZ[1][id] = -load3;
	fX[2][id] = load4;
	fX[5][id] = load4;

	f2 = trussFit3D(sol, &w2, &p2, verb, &feasible2, id);

	if(w1 != w2){
		printf("Error: different weights\n");
	}

	*weightRet = w1;
	*pen = (p1+p2)/2.0;
	*feasible = feasible1 & feasible2;

	return ( (f1+f2)/2.0 );
}

double fit72(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id){
	double f1, f2;
	double w1, w2;
	double p1, p2;
	bool feasible1, feasible2;

	//change to case 2 after first evaluation
	double load = 5.0;

	//reset load case 2
	fZ[16][id] = 0.0;
	fZ[17][id] = 0.0;
	fZ[18][id] = 0.0;
	fZ[19][id] = 0.0;

	//load case 1
	fX[16][id] = load;
	fY[16][id] = load;
	fZ[16][id] = -load;

	f1 = trussFit3D(sol, &w1, &p1, verb, &feasible1, id);

	//reseting load case 1
	fX[16][id] = 0.0;
	fY[16][id] = 0.0;
	fZ[16][id] = 0.0;

	//load case 2
	fZ[16][id] = -load;
	fZ[17][id] = -load;
	fZ[18][id] = -load;
	fZ[19][id] = -load;

	f2 = trussFit3D(sol, &w2, &p2, verb, &feasible2, id);

	if(w1 != w2){
		printf("Error: different weights\n");
	}

	*weightRet = w1;
	*pen = (p1+p2)/2.0;
	*feasible = feasible1 & feasible2;

	return ( (f1+f2)/2.0 );
}

double fit120(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id){
	double f1;
	double w1;
	double p1;
	bool feasible1;
	int groups[8] = {0, 12, 24, 36, 60, 84, 96, 120};
	double yield = 58.0; //ksi
	double cc, lambda, ri;
	double ar = 0.4993, br = 0.6777, k = 1.0;
	double temp1, temp2;
	double strs;

	/*
	groups:
	1: 0 - 11
	2: 12 - 23
	3: 24 - 35
	4: 36 - 59
	5: 60 - 83
	6: 84 - 95
	7: 96 - 119
	*/

	cc = sqrt(2*PI*PI*E/yield);

	for(int i=0;i<7;i++){
		ri = ar*pow(sol[i], br);
		lambda = (k*length[groups[i]+1])/ri;
		if(lambda < cc){
			temp1 = (1 - (lambda*lambda/(2*cc*cc)) )*yield;
			temp2 = (5/3) + (3*lambda/(8*cc)) - (lambda*lambda*lambda/(8*cc*cc*cc));
			strs = temp1/temp2;
		}
		else{
			strs = 12*PI*PI*E/(23*lambda*lambda);
		}
		if(strs>0.0){
			strs = -strs;
		}
		for(int j=groups[i]; j<groups[i+1];j++){
			minStress[j][id] = strs;
		}
	}

	f1 = trussFit3D(sol, &w1, &p1, verb, &feasible1, id);
	*weightRet = w1;
	*pen = p1;
	*feasible = feasible1;

	return f1;
}

double fit200(double *sol, double *weightRet, double *pen, bool verb, bool *feasible, int id){
	double f1, f2, f3;
	double w1, w2, w3;
	double p1, p2, p3;
	bool feasible1, feasible2, feasible3;

	double loadx = 1.0;
	double loady = -10.0;

	//case 1
	fX[0][id] = loadx;
	fX[5][id] = loadx;
	fX[14][id] = loadx;
	fX[19][id] = loadx;
	fX[28][id] = loadx;
	fX[33][id] = loadx;
	fX[42][id] = loadx;
	fX[47][id] = loadx;
	fX[56][id] = loadx;
	fX[61][id] = loadx;
	fX[70][id] = loadx;

	for(int i=0;i<6;i++){
		fY[i][id] = 0;
	}
	for(int i=7;i<12;i+=2){
		fY[i][id] = 0;
	}
	for(int i=13;i<20;i++){
		fY[i][id] = 0;
	}
	for(int i=21;i<26;i+=2){
		fY[i][id] = 0;
	}
	for(int i=27;i<34;i++){
		fY[i][id] = 0;
	}
	for(int i=35;i<40;i+=2){
		fY[i][id] = 0;
	}
	for(int i=41;i<48;i++){
		fY[i][id] = 0;
	}
	for(int i=49;i<54;i+=2){
		fY[i][id] = 0;
	}
	for(int i=55;i<62;i++){
		fY[i][id] = 0;
	}
	for(int i=63;i<68;i+=2){
		fY[i][id] = 0;
	}
	for(int i=69;i<75;i++){
		fY[i][id] = 0;
	}

	f1 = trussFit2D(sol, &w1, &p1, verb, &feasible1, id);

	//case 2
	fX[0][id] = 0;
	fX[5][id] = 0;
	fX[14][id] = 0;
	fX[19][id] = 0;
	fX[28][id] = 0;
	fX[33][id] = 0;
	fX[42][id] = 0;
	fX[47][id] = 0;
	fX[56][id] = 0;
	fX[61][id] = 0;
	fX[70][id] = 0;

	for(int i=0;i<6;i++){
		fY[i][id] = loady;
	}
	for(int i=7;i<12;i+=2){
		fY[i][id] = loady;
	}
	for(int i=13;i<20;i++){
		fY[i][id] = loady;
	}
	for(int i=21;i<26;i+=2){
		fY[i][id] = loady;
	}
	for(int i=27;i<34;i++){
		fY[i][id] = loady;
	}
	for(int i=35;i<40;i+=2){
		fY[i][id] = loady;
	}
	for(int i=41;i<48;i++){
		fY[i][id] = loady;
	}
	for(int i=49;i<54;i+=2){
		fY[i][id] = loady;
	}
	for(int i=55;i<62;i++){
		fY[i][id] = loady;
	}
	for(int i=63;i<68;i+=2){
		fY[i][id] = loady;
	}
	for(int i=69;i<75;i++){
		fY[i][id] = loady;
	}

	f2 = trussFit2D(sol, &w2, &p2, verb, &feasible2, id);

	//case 3
	fX[0][id] = loadx;
	fX[5][id] = loadx;
	fX[14][id] = loadx;
	fX[19][id] = loadx;
	fX[28][id] = loadx;
	fX[33][id] = loadx;
	fX[42][id] = loadx;
	fX[47][id] = loadx;
	fX[56][id] = loadx;
	fX[61][id] = loadx;
	fX[70][id] = loadx;

	for(int i=0;i<6;i++){
		fY[i][id] = loady;
	}
	for(int i=7;i<12;i+=2){
		fY[i][id] = loady;
	}
	for(int i=13;i<20;i++){
		fY[i][id] = loady;
	}
	for(int i=21;i<26;i+=2){
		fY[i][id] = loady;
	}
	for(int i=27;i<34;i++){
		fY[i][id] = loady;
	}
	for(int i=35;i<40;i+=2){
		fY[i][id] = loady;
	}
	for(int i=41;i<48;i++){
		fY[i][id] = loady;
	}
	for(int i=49;i<54;i+=2){
		fY[i][id] = loady;
	}
	for(int i=55;i<62;i++){
		fY[i][id] = loady;
	}
	for(int i=63;i<68;i+=2){
		fY[i][id] = loady;
	}
	for(int i=69;i<75;i++){
		fY[i][id] = loady;
	}

	f3 = trussFit2D(sol, &w3, &p3, verb, &feasible3, id);

	if(w1 != w2){
		printf("Error: different weights\n");
	}
	else if(w2 != w3){
		printf("Error: different weights\n");
	}

	*weightRet = w1;
	*pen = (p1+p2+p3)/3.0;
	*feasible = feasible1 & feasible2 & feasible3;

	return ( (f1+f2+f3)/3.0 );
}

double getMinArea(){
	return minArea;
}

double getMaxArea(){
	return maxArea;
}
int getDim(){
	return dim;
}