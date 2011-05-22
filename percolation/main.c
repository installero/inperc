#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void PrintHelp (char * pName); // Prints help output

#include "include/constants.h"
#include "include/init.h"
#include "include/defines.h"
#include "include/functions.c"

FILE * res;
FILE * dbg_obj;

int main(int argc, char * argv[]) {
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
