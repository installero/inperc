float * O2Od;

float F_rand(); // Returns random float number (0->1)

inline void StickRandomInit(float * Stick); // Inits stick's parameters with random values
inline void StickNullLength(float * Stick); // Sets stick's length to zero
inline void StickPrint(float * Stick); // Prints stick's parameters

inline float BondCriteria(float dR, float dE); // Simply returns equation for Eta

inline float StickToPointDistance (float * Stick, float * Point); // Finds distance between a stick and a point
inline float StickToBoundaryDistance (float * Stick, int Dir); // Finds distance from a stick to a boundary at a certain direction
inline int BelongsToBoundaryRegion(float * Stick, float BoundDist); // Returns 1 if a stick belongs to boundary region within bond distance
inline float StickToStickEnergy (float * Stick1, float * Stick2); // Finds distance between two sticks energies
inline float WildStickToStickDistance (float * Stick1, float * Stick2); // Finds distance between two sticks wildly
inline float RoughStickToStickDistance (float * Stick1, float * Stick2); // Finds distance between two sticks roughly
inline float StickToStickDistance (float * Stick1, float * Stick2); // Finds distance between two sticks
inline int CheckBondPresence(float * Stick1, float * Stick2, float BoundDist); // Checks whether the first stick touches the second

inline void FindDistances(float * Stick); // Finds distances between all the sticks in array and stores them to O2Od
inline int CheckPercolation(float * Site, float BoundDist, int Start, int End); // Checks whether percotation is achieved
inline float CountAverageBondsAmount(float * Site, float BoundDist); // Returns the total amount of bonds in the system

/* -- */

void StickNullLength(float * Stick) {
  Stick[Length]=0;
};

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

  if(ThreeDMode){
    K = ex1*ex2+ey1*ey2+ez1*ez2;
    Kp = dx*ex1+dy*ey1+dz*ez1;
    Kq = dx*ex2+dy*ey2+dz*ez2;
  }
  else {K=1;}

  if (K == 1) { // if the sticks are parallel
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
  return fmax(fabs(x1-x2),fmax(fabs(y1-y2),fabs(z1-z2)))-(W1+W2)-(L1+L2)/2;
};

float RoughStickToStickDistance (float * Stick1, float * Stick2) {
  static float dx, dy, dz;

  dx = x1-x2;
  dy = y1-y2;
  dz = z1-z2;
  
  return sqrt(dx*dx+dy*dy+dz*dz)-(W1+W2)-(L1+L2)/2;
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
  float dE = StickToStickEnergy(Stick1,Stick2);
  if (BondCriteria(WildStickToStickDistance(Stick1,Stick2),dE) < BoundDist) {
    if (BondCriteria(RoughStickToStickDistance(Stick1,Stick2),dE) < BoundDist) {
      if (BondCriteria(StickToStickDistance(Stick1,Stick2),dE) < BoundDist) return 1;
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
    }
    else invSticksAmount++;
	}

  return (float)res/(float)(ObjectNum-invSticksAmount);
};
