#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "init.h"

FILE * res;
FILE * dbg_obj;
char * resFileName;

float * O2Od;

float F_rand(); // Returns random float number (0->1)
void PrintHelp (char * pName); // Prints help output

void StickRandomInit(float * Stick); // Inits stick's parameters with random values
void StickPrint(float * Stick); // Prints stick's parameters
inline float StickToStickDistance (float * Stick1, float * Stick2); // Finds distance between two sticks
inline float StickToPointDistance (float * Stick, float * Point); // Finds distance between a stick and a point
inline float StickToBoundaryDistance (float * Stick, int Dir); // Finds distance from a stick to a boundary at a certain direction

int main(int argc, char * argv[]) {
  SetDefaultValues();

  int i;

  for (i=1;i < argc;i++) {
		if (argv[i][0] == '-') {
			switch ( argv[i][1] ) {
				case 'v':	{
					VerboseMode=1;
					break;
				};
				case 'c': {
					if ((i+1) < argc) BoundAccuracy=atof(argv[i+1]);
					else PrintHelp (argv[0]);
				  break;
				};
				case 'e': {
					if ((i+1) < argc) ExperimentNum=atoi(argv[i+1]);
					else PrintHelp (argv[0]);
					break;
				};
				case 'r': {
					if ((i+1) < argc) resFileName=argv[i+1];
					else PrintHelp (argv[0]);
					break;
				};
				case 'n': {
					if ((i+1) < argc) ObjectNum=atoi(argv[i+1]);
					else PrintHelp (argv[0]);
					break;
				};
				case 'w': {
					if ((i+1) < argc ) StickWidth = atof (argv[i+1]);
					else PrintHelp (argv[0]);
					break;
				};
				case 'l':	{
					if ((i+1) < argc) StickLength = atof(argv[i+1]);
					else PrintHelp (argv[0]);
					break;
				};
				case 'd': {
					switch (argv[i][2]) {
            case 'l':	{
              if ((i+1) < argc) StickLengthDistortion=atof(argv[i+1]);
					    else PrintHelp (argv[0]);
              break;
						};
            case 'w':	{
              if ((i+1) < argc) StickWidthDistortion=atof(argv[i+1]);
					    else PrintHelp (argv[0]);
              break;
						};
            case 't':	{
              if ((i+1) < argc)	StickThetaDistortion=pi*atof(argv[i+1]);
              else PrintHelp (argv[0]);
              break;
            };
            case 'f':	{
              if ((i+1) < argc)	StickFiDistortion=pi*atof(argv[i+1]);
              else PrintHelp (argv[0]);
              break;
            };
						default: PrintHelp (argv[0]);
					}
					break;
				};
				default: PrintHelp (argv[0]);
			};
		};
	};

  float Stick1[ParamsNum];
  float Stick2[ParamsNum];

  StickRandomInit(Stick1);
  StickRandomInit(Stick2);

  StickPrint(Stick1);
  StickPrint(Stick2);

  float * Stick;

  Stick = (float *) malloc(ObjectNum*ParamsNum);
  if (Stick==NULL) exit(1);

  for (i=0; i<ObjectNum*ParamsNum; i+=ParamsNum) {
    #define Stick_i Stick+i*sizeof(float)
    StickRandomInit(Stick_i);
    //StickPrint(Stick_i);
  }

  free (Stick);

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

/*
int CheckPercolation(float * Stick[ParamsNum], float BoundDist, int Start, int End)
{
	int InfCluster[ObjectNum];
	int InfClusterStickNum = 0;

  int StickState[ObjectNum];
  int i;

	for (i=0; i < SiteNumber; i++) StickState[i]=StateNew;

	if (Verbose) printf("\nThese objects touched the start (%d), BoundDist = %1.5f:\n", Start, BoundDist);

	for (i=0;i < SiteNumber;i++) {
		if (StickToBoundaryDistance(Stick[i], Start)<BoundDist)
		{
			if (Verbose) Site[i].Print();
			Site[i].History = InfClAlive;
			InfCluster[InfClusterStickNum] = i;
			InfClusterStickNum ++;
		}
	}

	for ( int k = 0; k < InfClusterStickNum; k++ )
	{
		int i = InfCluster[k];

		if ( Site[i].History == InfClAlive )
		{
			for ( int j = 0; j < ObjectNumber; j++ )
			{
				if ( !SeekLength )
				{
					if ( Site[j].History == New && ( O2Od[j*ObjectNumber + i] <= (BoundDist + 2*StickWidth) ) )
					{
						Site[j].History = InfClAlive;
						InfCluster[InfClusterStickNum] = j;
						InfClusterStickNum ++;
					}
				}
				else
				{
					if ( Site[j].History == New && ( Site[j].O2Odist(Site[i]) <= (BoundDist + 2*StickWidth) ) )
					{
						Site[j].History = InfClAlive;
						InfCluster[InfClusterStickNum] = j;
						InfClusterStickNum ++;
					}

				}

			}
			Site[i].History = InfClDead;
		}
		
		if ( Site[i].BoundaryTouch ( DirEnd, BoundDist + StickWidth) )
		{
			if ( Verbose )
			{
				printf ( "\nThis object touched the end (%d):\n", DirEnd );
				Site[i].Print();
			}
			return true;
		}
	}

	return false;
};
*/

void PrintHelp (char * pName)
{
	printf("Usage: %s [options]\n", pName);
	printf("\t-v\t\t\tdisplays verbose output [no]\n");
	printf("\t-e int\t\t\tnumber of realisations for each system [2]\n");
	printf("\t-n int\t\t\tnumber of objects [10]\n");
	printf("\t-w float\t\tstick width [0.0]\n");
	printf("\t-a float\t\tangle dispersion value (in pi) [0.5]\n");
	printf("\t-l float\t\tstick length value [0.0]\n");
	printf("\t-c float\t\tcritical bound dist accuracy [0.01]\n");
	printf("\t-sl\t\tseek for critical length on each situation [no]\n");
	printf("\t-su\t\tshow the critical av. bounds per node on each situation [no]\n");
	printf("\t-si float\t\tStick length distortion [0]\n");
	printf("\t-r FILE\t\t\tresults output file [results.txt]\n");
  printf("\t-h, --help\t\tshow this usage statement\n");
	exit (1);						
};
