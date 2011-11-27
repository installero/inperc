#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void PrintHelp (char * pName); // Prints help output

#include "include/constants.h"
#include "include/init.h"
#include "include/defines.h"
#include "include/functions.c"
#include "include/functions_test.c"

FILE * res;
FILE * dbg_obj;

int main(int argc, char * argv[]) {
  SetDefaultValues();
  SetParams(argc,argv);

  if(TestMode) {
    TestAllFunctions();
    return 0;
  };

  res = fopen (ResFileName, "a") ;

  float * Stick;
  Stick = malloc(ObjectNum*ParamsNum*sizeof(float));

  int e, i;

  float nu_av = 0, nu_disp = 0;

  float * nu;
  nu = malloc(ExperimentNum*sizeof(float));

  for (e=0; e<ExperimentNum; e++) {
    for (i=0; i<ObjectNum; i++) StickRandomInit(Stick_i);
    nu[e] = CountAverageBondsAmount(Stick, BondDistance);
    fprintf(stderr,"%1.5f, %1.5f\n", BondDistance, nu[e]);
  }

  nu_av = 0, nu_disp = 0;
  for (e=0; e<ExperimentNum; e++) nu_av += nu[e];
  nu_av = nu_av/(float)ExperimentNum;
  for (e=0; e<ExperimentNum; e++) nu_disp += pow(nu[e]-nu_av,2);
  nu_disp=sqrt(nu_disp/(float)ExperimentNum);

  fprintf(stderr,"%1.5f, %1.5f±%1.5f\r", BondDistance, nu_av, nu_disp);

  float rs = 0;
  if (ThreeDMode) rs = pow((3/(4*ObjectNum*pi)),0.33333);
  else rs = sqrt(1/(ObjectNum*pi));

  float nu_t = 0;
  if (ThreeDMode) {
    nu_t = ObjectNum*(
            (4*pi*BondDistance*BondDistance*BondDistance/3) +
            (2*pi*StickLength*BondDistance*BondDistance)
          );
  } else {
    nu_t = ObjectNum*(
            (pi*BondDistance*BondDistance) +
            (4*StickLength*BondDistance)+
            (StickFiDistortion ? StickLength*StickLength*(1-(sin(2*StickFiDistortion)/(2*StickFiDistortion)))/StickFiDistortion : 0)
          );
  }
  int color = (nu_t < nu_av + nu_disp && nu_t > nu_av - nu_disp ) ? 32 : 31;
  
  fprintf(res,"%%l\t%%n\t%%rs\t%%a\t%%df\t%%bd\t%%nu\t%%nu_err\t%%nu_t\n");
  fprintf(res,"%1.7f\t%d\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\n", StickLength, ObjectNum, rs, LocalizationLength, StickFiDistortion, BondDistance, nu_av, nu_disp, nu_t);

  fprintf(stderr,"\n%1.7f\t%d\t%1.7f\t%1.7f\t%1.7f\t%1.7f %1.7f±%1.7f", StickLength, ObjectNum, rs, LocalizationLength, StickFiDistortion, BondDistance, nu_av, nu_disp);
  fprintf(stderr, "%c[%d;%d;%dm (ν_t = %1.7f)", 0x1B, 1, color, 40, nu_t);
  fprintf(stderr, "%c[%dm\n", 0x1B, 0);

  
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
	printf("\t-df float\t\tstick azimuth angle dispersion [pi/2]\n");
	printf("\t-dt float\t\tstick polar angle dispersion [pi/2]\n");
	printf("\t-de float\t\tstick energy dispersion, meV [0]\n");
	printf("\t-r FILE\t\t\tresults output file [results.txt]\n");
  printf("\t-z\t\trun self-testing\n");
  printf("\t-h, --help\t\tshow this usage statement\n");
	exit (1);						
};
