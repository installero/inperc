#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void PrintHelp (char * pName); // Prints help output

#include "include/constants.h"
#include "include/init.h"

FILE * res;
FILE * dbg_obj;
float * O2Od;

#define Stick_i &Stick[i*ParamsNum]
#define Stick_j &Stick[j*ParamsNum]
#define Site_i &Site[i*ParamsNum]
#define Site_j &Site[j*ParamsNum]

#define O2Od_i_j O2Od[i*ObjectNum+j]
#define O2Od_j_i O2Od[j*ObjectNum+i]

#define x Stick[X]
#define y Stick[Y]
#define z Stick[Z]
#define L Stick[Length]
#define W Stick[Width]
#define Th Stick[Theta]
#define F Stick[Fi]
#define x0 Point[X]
#define y0 Point[Y]
#define z0 Point[Z]
#define W0 Point[Width]

#define x1 Stick1[X]
#define y1 Stick1[Y]
#define z1 Stick1[Z]
#define L1 Stick1[Length]
#define W1 Stick1[Width]
#define Th1 Stick1[Theta]
#define F1 Stick1[Fi]
#define E1 Stick1[Energy]
#define x2 Stick2[X]
#define y2 Stick2[Y]
#define z2 Stick2[Z]
#define L2 Stick2[Length]
#define W2 Stick2[Width]
#define Th2 Stick2[Theta]
#define F2 Stick2[Fi]
#define E2 Stick2[Energy]

float F_rand(); // Returns random float number (0->1)

inline void StickRandomInit(float * Stick); // Inits stick's parameters with random values
inline void StickPrint(float * Stick); // Prints stick's parameters

inline float StickToStickDistance (float * Stick1, float * Stick2); // Finds distance between two sticks
inline float WildStickToStickDistance (float * Stick1, float * Stick2); // Finds distance between two sticks roughly
inline float RoughStickToStickDistance (float * Stick1, float * Stick2); // Finds distance between two sticks roughly

inline float StickToStickEnergy (float * Stick1, float * Stick2); // Finds distance between two sticks energies

inline float StickToPointDistance (float * Stick, float * Point); // Finds distance between a stick and a point
inline float StickToBoundaryDistance (float * Stick, int Dir); // Finds distance from a stick to a boundary at a certain direction

inline void FindDistances(float * Stick); // Finds distances between all the sticks in array and stores them to O2Od

inline float BondCriteria(float dR, float dE); // Simply returns equation for Eta
inline int CheckPercolation(float * Site, float BoundDist, int Start, int End); // Checks whether percotation is achieved
inline int CheckBondPresence(float * Stick1, float * Stick2, float BoundDist); // Checks whether the first stick touches the second

inline float CountAverageBondsAmount(float * Site, float BoundDist); // Returns the total amount of bonds in the system
inline int BelongsToBoundaryRegion(float * Stick, float BoundDist); // Returns 1 if stick belongs to boundary region within bond distance

int main(int argc, char * argv[]) {
  // Experiments zone

  /*static float Point_test[ParamsNum];*/
  /*static float Stick_test[ParamsNum];*/

  /*Point_test[X] = 0;*/
  /*Point_test[Y] = 0;*/
  /*Point_test[Z] = 0;*/
  /*Point_test[Length] = 0;*/

  /*Stick_test[X] = 2;*/
  /*Stick_test[Y] = 2;*/
  /*Stick_test[Z] = 2;*/
  /*Stick_test[Theta] = (3.1415)/2;*/
  /*Stick_test[Fi] = -(3.1415)/4;*/
  /*Stick_test[Length] = 1;*/

  /*printf("%1.5f\n", StickToPointDistance(Stick_test,Point_test));*/

  /*return 0;*/

  // end


  SetDefaultValues();
  SetParams(argc,argv);

  res = fopen ( ResFileName, "a" ) ;

  float * Stick;
  Stick = malloc(ObjectNum*ParamsNum*sizeof(float));

  int i,e, percolation_x = 0, percolation_y = 0, percolation_z = 0;

  float rs = 0;
  float BoundStep = 0, StartBoundStep = 0;
  float Eta_c_av = 0, Eta_c_disp = 0, nu_c_av = 0, nu_c_disp = 0;
  float BoundDistNew_x = 0, BoundDist_x = 0, Eta_c_x_av = 0, Eta_c_x_disp = 0, nu_c_x_av = 0, nu_c_x_disp = 0;
  float BoundDistNew_y = 0, BoundDist_y = 0, Eta_c_y_av = 0, Eta_c_y_disp = 0, nu_c_y_av = 0, nu_c_y_disp = 0;
  float BoundDistNew_z = 0, BoundDist_z = 0, Eta_c_z_av = 0, Eta_c_z_disp = 0, nu_c_z_av = 0, nu_c_z_disp = 0;

  float * Eta_c_x, * Eta_c_y, * Eta_c_z;
  Eta_c_x = malloc(ExperimentNum*sizeof(float));
  Eta_c_y = malloc(ExperimentNum*sizeof(float));
  Eta_c_z = malloc(ExperimentNum*sizeof(float));

  float * nu_c_x, * nu_c_y, * nu_c_z;
  nu_c_x = malloc(ExperimentNum*sizeof(float));
  nu_c_y = malloc(ExperimentNum*sizeof(float));
  nu_c_z = malloc(ExperimentNum*sizeof(float));

  if (ThreeDMode) {
    rs = pow((3/(4*ObjectNum*pi)),0.33333);
    StartBoundStep = 1.5f*2.0f*rs/LocalizationLength + StickEnergyDistortion/Temperature;
  } else {
    rs = sqrt(1/(ObjectNum*pi));
    StartBoundStep = 2.2f*2.0f*rs/LocalizationLength + StickEnergyDistortion/Temperature;
  }

  for (e=0; e<ExperimentNum; e++) {
    for (i=0; i<ObjectNum; i++) StickRandomInit(Stick_i);
    if(VerboseMode) {
      for (i=0; i<ObjectNum; i++) {
        printf("%d: ",i);
        StickPrint(Stick_i);
      }
    }

    BoundStep = StartBoundStep;
    BoundDistNew_x = BoundStep;
    BoundDistNew_y = BoundStep;
    BoundDistNew_z = BoundStep;

    while (BoundStep > BoundAccuracy) {
      BoundStep /= 2;
      BoundDist_x = BoundDistNew_x;
      BoundDist_y = BoundDistNew_y;
      BoundDist_z = BoundDistNew_z;

      if (ThreeDMode) {
        fprintf(stderr, "\rCritical radius acuracy: %1.5f, Eta_c_x: %1.5f, Eta_c_y: %1.5f, Eta_c_z: %1.5f", BoundStep, BoundDist_x, BoundDist_y, BoundDist_z);
      } else {
        fprintf(stderr, "\rCritical radius acuracy: %1.5f, Eta_c_x: %1.5f, Eta_c_y: %1.5f", BoundStep, BoundDist_x, BoundDist_y);
      }

      percolation_x = CheckPercolation(Stick, BoundDist_x, MIN_X, MAX_X);
      percolation_y = CheckPercolation(Stick, BoundDist_y, MIN_Y, MAX_Y);
      if (ThreeDMode) percolation_z = CheckPercolation(Stick, BoundDist_z, MIN_Z, MAX_Z);

      if (!percolation_x) BoundDistNew_x += BoundStep;
      else BoundDistNew_x -= BoundStep;
      if (!percolation_y) BoundDistNew_y += BoundStep;
      else BoundDistNew_y -= BoundStep;
      if (!percolation_z) BoundDistNew_z += BoundStep;
      else BoundDistNew_z -= BoundStep;
    }

    Eta_c_x[e] = BoundDist_x;
    Eta_c_y[e] = BoundDist_y;
    Eta_c_z[e] = BoundDist_z;

    nu_c_x[e] = CountAverageBondsAmount(Stick, BoundDist_x);
    nu_c_y[e] = CountAverageBondsAmount(Stick, BoundDist_y);
    if (ThreeDMode) nu_c_z[e] = CountAverageBondsAmount(Stick, BoundDist_z);

    fprintf(stderr,"\n");
  }

  for (e=0; e<ExperimentNum; e++) {
    Eta_c_x_av += Eta_c_x[e];
    Eta_c_y_av += Eta_c_y[e];
    Eta_c_z_av += Eta_c_z[e];

    nu_c_x_av += nu_c_x[e];
    nu_c_y_av += nu_c_y[e];
    nu_c_z_av += nu_c_z[e];
  }

  Eta_c_x_av=Eta_c_x_av/(float)ExperimentNum;
  Eta_c_y_av=Eta_c_y_av/(float)ExperimentNum;
  Eta_c_z_av=Eta_c_z_av/(float)ExperimentNum;

  nu_c_x_av=nu_c_x_av/(float)ExperimentNum;
  nu_c_y_av=nu_c_y_av/(float)ExperimentNum;
  nu_c_z_av=nu_c_z_av/(float)ExperimentNum;

  for (e=0; e<ExperimentNum; e++) {
    Eta_c_x_disp += pow(Eta_c_x[e]-Eta_c_x_av,2);
    Eta_c_y_disp += pow(Eta_c_y[e]-Eta_c_y_av,2);
    Eta_c_z_disp += pow(Eta_c_z[e]-Eta_c_z_av,2);

    nu_c_x_disp += pow(nu_c_x[e]-nu_c_x_av,2);
    nu_c_y_disp += pow(nu_c_y[e]-nu_c_y_av,2);
    nu_c_z_disp += pow(nu_c_z[e]-nu_c_z_av,2);
  }

  Eta_c_x_disp=sqrt(Eta_c_x_disp/(float)ExperimentNum);
  Eta_c_y_disp=sqrt(Eta_c_y_disp/(float)ExperimentNum);
  Eta_c_z_disp=sqrt(Eta_c_z_disp/(float)ExperimentNum);

  nu_c_x_disp=sqrt(nu_c_x_disp/(float)ExperimentNum);
  nu_c_y_disp=sqrt(nu_c_y_disp/(float)ExperimentNum);
  nu_c_z_disp=sqrt(nu_c_z_disp/(float)ExperimentNum);

  if (ThreeDMode) {
    Eta_c_av = (Eta_c_x_av + Eta_c_y_av + Eta_c_z_av)/3;
    Eta_c_disp = (Eta_c_x_disp + Eta_c_y_disp + Eta_c_z_disp)/3;

    nu_c_av = (nu_c_x_av + nu_c_y_av + nu_c_z_av)/3;
    nu_c_disp = (nu_c_x_disp + nu_c_y_disp + nu_c_z_disp)/3;
  } else {
    Eta_c_av = (Eta_c_x_av + Eta_c_y_av)/2;
    Eta_c_disp = (Eta_c_x_disp + Eta_c_y_disp)/2;

    nu_c_av = (nu_c_x_av + nu_c_y_av)/2;
    nu_c_disp = (nu_c_x_disp + nu_c_y_disp)/2;
  }

  fprintf (stderr,"Eta_c: %1.5f±%1.5f, r_s: %1.5f, nu_c: %1.5f±%1.5f\n", Eta_c_av, Eta_c_disp, rs, nu_c_av, nu_c_disp);

	fprintf (res, "%%e\t%%t\t%%n\t%%l\t%%dl\t%%w\t%%dw\t%%df\t%%dt\t%%Eta_c\t%%Eta_c_e\t%%Rc\t%%Rc_e\t%%rs\t%%Eta_c/rs\t%%Eta_c/rs_e\t%%nu_c\t%%nu_c_e\n" );
	fprintf (res, "%d\t%1.5f\t%d\t", ExperimentNum, Temperature, ObjectNum);
	fprintf (res, "%1.5f\t%1.5f\t", StickLength, StickLengthDistortion);
	fprintf (res, "%1.5f\t%1.5f\t", StickWidth, StickWidthDistortion);
	fprintf (res, "%1.5f\t%1.5f\t", StickFiDistortion, StickThetaDistortion);
	fprintf (res, "%1.7f\t%1.7f\t", Eta_c_av, Eta_c_disp);
	fprintf (res, "%1.7f\t%1.7f\t", LocalizationLength*(Eta_c_av)/2, LocalizationLength*(Eta_c_disp)/2);
	fprintf (res, "%1.7f\t", rs);
	fprintf (res, "%1.7f\t%1.7f\t", Eta_c_av/rs, Eta_c_disp/rs);
	fprintf (res, "%1.7f\t%1.7f\t", nu_c_av, nu_c_disp);
	fprintf (res, "\n");

  return 0;
}

void StickRandomInit(float * Stick) {
  Stick[X]=MIN_XSide+(MAX_XSide-MIN_XSide)*F_rand();
  Stick[Y]=MIN_YSide+(MAX_YSide-MIN_YSide)*F_rand();
  Stick[Length]=(StickLength-StickLengthDistortion)+2*StickLengthDistortion*F_rand();
  Stick[Width]=(StickWidth-StickWidthDistortion)+2*StickWidthDistortion*F_rand();
  Stick[Fi]=(2*F_rand()-1)*StickFiDistortion;
  Stick[Energy]=StickEnergyDistortion*F_rand();
  if (ThreeDMode) {
    Stick[Z]=MIN_ZSide+(MAX_ZSide-MIN_ZSide)*F_rand();
    Stick[Theta]=pi/2-StickThetaDistortion*(2*F_rand()-1);
  } else {
    Stick[Z]=0;
    Stick[Theta]=pi/2; 
  }
};

void StickPrint(float * Stick) {
  int i;
  printf("|");
  for (i=0; i<ParamsNum; i++) printf("%1.5f|",Stick[i]);
  printf("\n");
};

float F_rand() {return (float)rand()/RAND_MAX;};

float StickToStickDistance (float * Stick1, float * Stick2) {
  static float dx, dy, dz, ex1, ex2, ey1, ey2, ez1, ez2, K, Kp, Kq, p, q, Dx, Dy, Dz;

  dx = x1-x2;
  dy = y1-y2;
  dz = z1-z2;

  ex1 = sin(Th1)*cos(F1);
  ey1 = sin(Th1)*sin(F1);
  ez1 = cos(Th1);
  ex2 = sin(Th2)*cos(F2);
  ey2 = sin(Th2)*sin(F2);
  ez2 = cos(Th2);

  K = ex1*ex2+ey1*ey2+ez1*ez2;
  Kp = dx*ex1+dy*ey1+dz*ez1;
  Kq = dx*ex2+dy*ey2+dz*ez2;

  if (K == 1) { // if the stick's are parallel
    static float Distance_tmp, Distance_new;
    static float Point_tmp[ParamsNum];

    StickRandomInit(Point_tmp);
    Point_tmp[Length] = 0;
    Distance_tmp = Distance_MAX;

    Point_tmp[X] = x1+L1*ex1/2;
    Point_tmp[Y] = y1+L1*ey1/2;
    Point_tmp[Z] = z1+L1*ez1/2;
    Point_tmp[Width] = W1;
    Distance_new = StickToPointDistance(Stick2,Point_tmp);

    if (Distance_new < Distance_tmp) Distance_tmp=Distance_new;

    Point_tmp[X] = x1-L1*ex1/2;
    Point_tmp[Y] = y1-L1*ey1/2;
    Point_tmp[Z] = z1-L1*ez1/2;
    Distance_new = StickToPointDistance(Stick2,Point_tmp);

    if (Distance_new < Distance_tmp) Distance_tmp=Distance_new;

    Point_tmp[X] = x2+L2*ex2/2;
    Point_tmp[Y] = y2+L2*ey2/2;
    Point_tmp[Z] = z2+L2*ez2/2;
    Point_tmp[Width] = W2;
    Distance_new = StickToPointDistance(Stick1,Point_tmp);

    if (Distance_new < Distance_tmp) Distance_tmp=Distance_new;

    return Distance_tmp;
  }
  else { // if the stick's are not parallel
    if (L1) p=2*(K*Kq-Kp)/(L1*(1-K*K));
    else p=0;
    if (L2) q=2*(Kq-K*Kp)/(L2*(1-K*K));
    else q=0;
    
    if (p > 1) p=1;
    else if (p < -1) p=-1;
    if (q > 1) q=1;
    else if (q < -1) q=-1;

    Dx=dx+L1*p*ex1/2-L2*q*ex2/2;
    Dy=dy+L1*p*ey1/2-L2*q*ey2/2;
    Dz=dz+L1*p*ez1/2-L2*q*ez2/2;

    return sqrt(Dx*Dx+Dy*Dy+Dz*Dz)-(W1+W2);
  };
};

float WildStickToStickDistance (float * Stick1, float * Stick2) {
  return fmax(fabs(x1-x2),fmax(fabs(y1-y2),fabs(z1-z2)))-(W1+W2)-(L1+L2);
};

float RoughStickToStickDistance (float * Stick1, float * Stick2) {
  static float dx, dy, dz;

  dx = x1-x2;
  dy = y1-y2;
  dz = z1-z2;
  
  return sqrt(dx*dx+dy*dy+dz*dz)-(W1+W2)-(L1+L2);
};

float StickToStickEnergy (float * Stick1, float * Stick2) {
  return (fabs(E1-E2)+fabs(E1-StickEnergyDistortion/2)+(E2-StickEnergyDistortion/2))/2;
};

float StickToPointDistance (float * Stick, float * Point) {

  static float dx, dy, dz, ex, ey, ez, K, p, Dx, Dy, Dz;

  dx=x-x0;
  dy=y-y0;
  dz=z-z0;
  
  ex=sin(Th)*cos(F);
  ey=sin(Th)*sin(F);
  ez=cos(Th);

  K=dx*ex+dy*ey+dz*ez;
  if (L) p=-2*K/L;
  else p=0;
  
  if (p > 1) p=1;
  else if (p < -1) p=-1;

  Dx=dx+L*p*ex/2;
  Dy=dy+L*p*ey/2;
  Dz=dz+L*p*ez/2;

  return sqrt(Dx*Dx+Dy*Dy+Dz*Dz)-(W0+W);
};

inline float StickToBoundaryDistance (float * Stick, int Dir) {
  if (Dir == MIN_X) return (x-W-fabs(L*sin(Th)*cos(F)/2))-MIN_XSide;
  if (Dir == MAX_X) return MAX_XSide-(x+W+fabs(L*sin(Th)*cos(F)/2));
  if (Dir == MIN_Y) return (y-W-fabs(L*sin(Th)*sin(F)/2))-MIN_YSide;
  if (Dir == MAX_Y) return MAX_YSide-(y+W+fabs(L*sin(Th)*sin(F)/2));
  if (Dir == MIN_Z) return (z-W-fabs(L*cos(Th)/2))-MIN_ZSide;
  if (Dir == MAX_Z) return MAX_ZSide-(z+W+fabs(L*cos(Th)/2));
  fprintf (stderr,"Error: StickToBoundaryDistance, Unknown direction!");
  exit(1);
  return 0;
}

void FindDistances(float * Stick) {
  int i, j;
  for (i=0; i < ObjectNum; i++) {
    fprintf ( stderr, "\rCounting distances: %1.2f%%", (float)(i*100)/(float)(ObjectNum));
		for (j=0; j<ObjectNum; j++) {
			if (i == j) O2Od_i_j=0;
			else if (i > j) O2Od_i_j=O2Od_j_i;
			else O2Od_i_j=StickToStickDistance(Stick_i,Stick_j);
		}
	}
}

int CheckBondPresence(float * Stick1, float * Stick2, float BoundDist) {
  if (BondCriteria(WildStickToStickDistance(Stick1,Stick2),StickToStickEnergy(Stick1,Stick2)) < BoundDist) {
    if (BondCriteria(RoughStickToStickDistance(Stick1,Stick2),StickToStickEnergy(Stick1,Stick2)) < BoundDist) {
      if (BondCriteria(StickToStickDistance(Stick1,Stick2),StickToStickEnergy(Stick1,Stick2)) < BoundDist) {
        if (VerboseMode) {
          printf("Touch (Distanse = %1.5f, Energy = %1.5f, Criteria = %1.5f)!\n",
            StickToStickDistance(Stick1,Stick2),StickToStickEnergy(Stick1,Stick2),
            BondCriteria(StickToStickDistance(Stick1,Stick2),StickToStickEnergy(Stick1,Stick2))
          );
          StickPrint(Stick1);
          StickPrint(Stick2);
        }
        return 1;
      }
    }
  }
  return 0;
}

float BondCriteria(float dR, float dE) {
  return 2*dR/LocalizationLength + dE/Temperature;
};

int CheckPercolation(float * Site, float BoundDist, int Start, int End)
{
  int i,j,k;

  int SiteState[ObjectNum];
	for (i=0; i < ObjectNum; i++) SiteState[i]=StateNew;

	int InfCluster[ObjectNum];
	int InfClusterStickNum = 0;

	if (VerboseMode) printf("Start:%d, BoundDist:%1.5f\n", Start, BoundDist);

	for (i=0;i < ObjectNum;i++) {
		if (BondCriteria(StickToBoundaryDistance(Site_i, Start),0) < BoundDist) {
			if (VerboseMode) StickPrint(Site_i);
			SiteState[i] = StateInfClNew;
			InfCluster[InfClusterStickNum] = i;
			InfClusterStickNum ++;
		}
	}
 
	for (k=0; k < InfClusterStickNum; k++) {
		i = InfCluster[k];

		if (BondCriteria(StickToBoundaryDistance(Site_i, End),0) < BoundDist) {
			if (VerboseMode) {
				printf ("End (%d):\n", End);
        StickPrint(Site_i);
			}
			return 1;
		}

		if (SiteState[i] == StateInfClNew) {
			for (j=0; j < ObjectNum; j++) {
				if (SiteState[j] == StateNew) {
          if (CheckBondPresence(Site_i,Site_j,BoundDist)) {
            SiteState[j] = StateInfClNew;
            InfCluster[InfClusterStickNum] = j;
            InfClusterStickNum ++;
          }
				}
			}
			SiteState[i] = StateInfClOld;
		}
	}

	return 0;
};

int BelongsToBoundaryRegion(float * Stick, float BoundDist) {
  float DoubleBoundDist = 2*BoundDist;
  if (BondCriteria(StickToBoundaryDistance(Stick, MIN_X),0) < DoubleBoundDist) return 1;
  if (BondCriteria(StickToBoundaryDistance(Stick, MAX_X),0) < DoubleBoundDist) return 1;
  if (BondCriteria(StickToBoundaryDistance(Stick, MIN_Y),0) < DoubleBoundDist) return 1;
  if (BondCriteria(StickToBoundaryDistance(Stick, MAX_Y),0) < DoubleBoundDist) return 1;
  if (ThreeDMode) {
    if (BondCriteria(StickToBoundaryDistance(Stick, MIN_Z),0) < DoubleBoundDist) return 1;
    if (BondCriteria(StickToBoundaryDistance(Stick, MAX_Z),0) < DoubleBoundDist) return 1;
  }
  return 0;
};

float CountAverageBondsAmount(float * Site, float BoundDist) {
  int res = 0;
  int i,j;
  int invSticksAmount = 0;

	for (i = 0; i < ObjectNum; i++) {
    if(!BelongsToBoundaryRegion(Site_i, BoundDist)) {
      for (j = 0; j < ObjectNum; j++) {
        if (j != i) {
          if (CheckBondPresence(Site_i, Site_j, BoundDist)) res++;
        }
      }
    } else {
      invSticksAmount++;
    }
	}

  return (float)res/(float)(ObjectNum-invSticksAmount);
};

void PrintHelp (char * pName)
{
	printf("Usage: %s [options]\n", pName);
	printf("\t-v\t\t\tdisplay verbose output [off]\n");
	printf("\t-3\t\t\t3 dimension system mode [off]\n");
	printf("\t-e int\t\t\tnumber of realisations [1]\n");
	printf("\t-n int\t\t\tnumber of objects [10]\n");
	printf("\t-a int\t\t\tlocalization length, 10^4A [0.01]\n");
	printf("\t-t int\t\t\ttemperature,meV [10]\n");
	printf("\t-c float\t\tcritical bound dist accuracy [0.1]\n");
	printf("\t-l float\t\tstick length, 10^4A [0.0]\n");
	printf("\t-dl float\t\tstick leng dispersion, 10^4A [0.0]\n");
	printf("\t-w float\t\tstick width, 10^4A [0.0]\n");
	printf("\t-dw float\t\tstick width dispersion, 10^4A [0.0]\n");
	printf("\t-df float\t\tstick azimuth angle dispersion [pi]\n");
	printf("\t-dt float\t\tstick polar angle dispersion [pi/2]\n");
	printf("\t-de float\t\tstick energy dispersion, meV [0]\n");
	printf("\t-r FILE\t\t\tresults output file [results.txt]\n");
  printf("\t-h, --help\t\tshow this usage statement\n");
	exit (1);						
};
