// System
#define MIN_XSide 0.0f
#define MAX_XSide 1.0f
#define MIN_YSide 0.0f
#define MAX_YSide 1.0f
#define MIN_ZSide 0.0f
#define MAX_ZSide 1.0f

#define Distance_MAX (MAX_XSide - MIN_XSide + MAX_YSide - MIN_YSide + MAX_ZSide - MIN_ZSide)

// Stick states
#define StateNew		      0
#define StateInfClNew	    1
#define StateInfClOld	    2

// Experiment
#define VerboseModeDef		          0
#define ThreDModeDef		                0

#define StickLengthDef		          0.0f
#define StickLengthDistortionDef		0.0f
#define StickWidthDef		            0.0f
#define StickWidthDistortionDef		  0.0f
#define StickThetaDistortionDef	    pi*0.5f
#define StickFiDistortionDef	      pi*1.0f

#define TemperatureDef              10
#define StickEnergyDistortionDef    0
#define LocalizationLengthDef       0.01

#define ExperimentNumDef	          1
#define ObjectNumDef		            10
#define BoundAccuracyDef            0.1f
#define ResFileNameDef		          "results.txt"

int VerboseMode;			          // Verbose mode on/off
int ThreeDMode;			            // 3 dimensions mode on/off

float StickLength;			        // Stick length
float StickLengthDistortion;		// Stick length is distorted with that dispercy
float StickWidth;			          // Stick width
float StickWidthDistortion;		  // Stick width is distorted with that dispercy
float StickThetaDistortion;     // Stick polar angle distortion
float StickFiDistortion;        // Stick azimuth angle distortion

float Temperature;              // Systems temperature
float StickEnergyDistortion;    // Density of states width
float LocalizationLength;       // Wave function decay radio

int ExperimentNum;		          // Number of realisations of each system
int ObjectNum;		              // Number of objects
float BoundAccuracy;            // Minimum value that two bound distance steps can vary one from each other
char * ResFileName;             // Results file name

void SetDefaultValues();
void SetParams();

#include "init.c"
