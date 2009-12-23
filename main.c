#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void PrintHelp (char * pName); // Prints help output

#include "constants.h"
#include "init.h"

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
#define x2 Stick2[X]
#define y2 Stick2[Y]
#define z2 Stick2[Z]
#define L2 Stick2[Length]
#define W2 Stick2[Width]
#define Th2 Stick2[Theta]
#define F2 Stick2[Fi]

float F_rand(); // Returns random float number (0->1)

inline void StickRandomInit(float * Stick); // Inits stick's parameters with random values
inline void StickPrint(float * Stick); // Prints stick's parameters
inline float StickToStickDistance (float * Stick1, float * Stick2); // Finds distance between two sticks
inline float StickToPointDistance (float * Stick, float * Point); // Finds distance between a stick and a point
inline float StickToBoundaryDistance (float * Stick, int Dir); // Finds distance from a stick to a boundary at a certain direction
inline void FindDistances(float * Stick); // Finds distances between all the sticks in array and stores them to O2Od
inline int CheckPercolation(float * Site, float BoundDist, int Start, int End); // Checks whether percotation is achieved

int main(int argc, char * argv[]) {
  SetDefaultValues();
  SetParams(argc,argv);

  float * Stick;
  Stick = malloc(ObjectNum*ParamsNum*sizeof(float));

  O2Od = malloc(ObjectNum*ObjectNum*sizeof(float));

  int i,e, percolation_x, percolation_y, percolation_z;

  float BoundStep;
  float BoundDistNew_x, BoundDist_x, Rc_x_av = 0;
  float BoundDistNew_y, BoundDist_y, Rc_y_av = 0;
  float BoundDistNew_z, BoundDist_z, Rc_z_av = 0;

  float * Rc_x, * Rc_y, * Rc_z;
  Rc_x = malloc(ExperimentNum*sizeof(float));
  Rc_y = malloc(ExperimentNum*sizeof(float));
  Rc_z = malloc(ExperimentNum*sizeof(float));

  for (e=0; e<ExperimentNum; e++) {
    for (i=0; i<ObjectNum; i++) StickRandomInit(Stick_i);
    if(VerboseMode) {
      for (i=0; i<ObjectNum; i++) {
        printf("%d: ",i);
        StickPrint(Stick_i);
      }
    }

    FindDistances(Stick);

    percolation_x = 0;
    percolation_y = 0;
    percolation_z = 0;
    BoundStep = (float)(Distance_MAX)/2;
    BoundDistNew_x = BoundStep;
    BoundDistNew_y = BoundStep;
    BoundDistNew_z = BoundStep;
    BoundDist_x = 0;
    BoundDist_y = 0;
    BoundDist_z = 0;

    while (BoundStep > BoundAccuracy) {
      BoundStep /= 2;
      BoundDist_x = BoundDistNew_x;
      BoundDist_y = BoundDistNew_y;
      BoundDist_z = BoundDistNew_z;

      fprintf(stderr, "\rCounting distances: 100%%, Critical radius acuracy: %1.5f, B_x: %1.5f, B_y: %1.5f, B_z: %1.5f", BoundStep, BoundDist_x, BoundDist_y, BoundDist_z);

      percolation_x = CheckPercolation(Stick, BoundDist_x, MIN_X, MAX_X);
      percolation_y = CheckPercolation(Stick, BoundDist_y, MIN_Y, MAX_Y);
      percolation_z = CheckPercolation(Stick, BoundDist_z, MIN_Z, MAX_Z);

      if (!percolation_x) BoundDistNew_x += BoundStep;
      else BoundDistNew_x -= BoundStep;
      if (!percolation_y) BoundDistNew_y += BoundStep;
      else BoundDistNew_y -= BoundStep;
      if (!percolation_z) BoundDistNew_z += BoundStep;
      else BoundDistNew_z -= BoundStep;
    }

    Rc_x[e] = BoundDist_x;
    Rc_y[e] = BoundDist_y;
    Rc_z[e] = BoundDist_z;

    fprintf(stderr,"\n");
  }

  for (e=0; e<ExperimentNum; e++) {
    Rc_x_av += Rc_x[e]/(float)ExperimentNum;
    Rc_y_av += Rc_y[e]/(float)ExperimentNum;
    Rc_z_av += Rc_z[e]/(float)ExperimentNum;
  }

  fprintf (stderr,"Rc_x_av: %1.5f, Rc_y_av: %1.5f, Rc_z_av: %1.5f\n", Rc_x_av, Rc_y_av, Rc_z_av); 

  return 0;
}

void StickRandomInit(float * Stick) {
  Stick[X] = MIN_XSide + (MAX_XSide - MIN_XSide)*F_rand();
  Stick[Y] = MIN_YSide + (MAX_YSide - MIN_YSide)*F_rand();
  Stick[Z] = MIN_ZSide + (MAX_ZSide - MIN_ZSide)*F_rand();
  Stick[Length] = (StickLength - StickLengthDistortion) + 2*StickLengthDistortion*F_rand();
  Stick[Width] = (StickWidth - StickWidthDistortion) + 2*StickWidthDistortion*F_rand();
  Stick[Theta] = (pi/2 - StickThetaDistortion) + 2*StickThetaDistortion*F_rand();
  Stick[Fi] = (pi - StickFiDistortion) + 2*StickFiDistortion*F_rand();
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
  ex2 = sin(Th2)*cos(F2);
  ey1 = sin(Th1)*sin(F1);
  ey2 = sin(Th2)*sin(F2);
  ez1 = cos(Th1);
  ez2 = cos(Th2);

  K = ex1*ex2 + ey1*ey2 + ez1*ez2;
  Kp = dx*ex1 + dy*ey1 + dz*ez1;
  Kq = dx*ex2 + dy*ey2 + dz*ez2;

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
    if (L2) q=2*(K*Kp-Kq)/(L2*(1-K*K));
    
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
  fprintf (stderr,"Error: StickToDoundaryDistance, Unknown direction!");
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

int CheckPercolation(float * Site, float BoundDist, int Start, int End)
{
  int i,j,k;

  int SiteState[ObjectNum];
	for (i=0; i < ObjectNum; i++) SiteState[i]=StateNew;

	int InfCluster[ObjectNum];
	int InfClusterStickNum = 0;

	if (VerboseMode) printf("Start:%d, BoundDist:%1.5f\n", Start, BoundDist);

	for (i=0;i < ObjectNum;i++) {
		if (StickToBoundaryDistance(Site_i, Start) < BoundDist) {
			if (VerboseMode) StickPrint(Site_i);
			SiteState[i] = StateInfClNew;
			InfCluster[InfClusterStickNum] = i;
			InfClusterStickNum ++;
		}
	}
 
	for (k=0; k < InfClusterStickNum; k++) {
		i = InfCluster[k];

		if (StickToBoundaryDistance(Site_i, End) < BoundDist) {
			if (VerboseMode) {
				printf ("End (%d):\n", End);
        StickPrint(Site_i);
			}
			return 1;
		}

		if (SiteState[i] == StateInfClNew) {
			for (j=0; j < ObjectNum; j++) {
					if (SiteState[j] == StateNew && (O2Od_i_j < BoundDist)) {
            if (VerboseMode) {
              printf("Touch!\n");
              StickPrint(Site_i);
              StickPrint(Site_j);
            }
            SiteState[j] = StateInfClNew;
						InfCluster[InfClusterStickNum] = j;
						InfClusterStickNum ++;
					}
			}
			SiteState[i] = StateInfClOld;
		}
	}

	return 0;
};

void PrintHelp (char * pName)
{
	printf("Usage: %s [options]\n", pName);
	printf("\t-v\t\t\tdisplay verbose output [off]\n");
	printf("\t-e int\t\t\tnumber of realisations [1]\n");
	printf("\t-n int\t\t\tnumber of objects [10]\n");
	printf("\t-c float\t\tcritical bound dist accuracy [0.1]\n");
	printf("\t-l float\t\tstick length [0.0]\n");
	printf("\t-dl float\t\tstick leng dispersion [0.0]\n");
	printf("\t-w float\t\tstick width [0.0]\n");
	printf("\t-dw float\t\tstick width dispersion [0.0]\n");
	printf("\t-df float\t\tstick azimuth angle dispersion [pi]\n");
	printf("\t-dt float\t\tstick polar angle dispersion [pi/2]\n");
	printf("\t-r FILE\t\t\tresults output file [results.txt]\n");
  printf("\t-h, --help\t\tshow this usage statement\n");
	exit (1);						
};
