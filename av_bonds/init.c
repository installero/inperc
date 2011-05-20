void SetDefaultValues() {
  // Setting modes
	VerboseMode = VerboseModeDef;
  ThreeDMode = ThreDModeDef;

  // Setting default stick parameters 
	StickLength = StickLengthDef;
  StickLengthDistortion = StickLengthDistortionDef;
	StickWidth = StickWidthDef;
  StickWidthDistortion = StickWidthDistortionDef;
  StickThetaDistortion = StickThetaDistortionDef;
  StickFiDistortion = StickFiDistortionDef;

  // Setting experiment parameters
  ExperimentNum = ExperimentNumDef;
	ObjectNum = ObjectNumDef;
	BoundAccuracy = BoundAccuracyDef;
  ResFileName = ResFileNameDef;

  Temperature = TemperatureDef;
  StickEnergyDistortion = StickEnergyDistortionDef;
  LocalizationLength = LocalizationLengthDef;

  BondDistance = BondDistanceDef;
}

void SetParams(int argc, char * argv[]) {
  int i;
  for (i=1;i < argc;i++) {
		if (argv[i][0] == '-') {
			switch ( argv[i][1] ) {
				case 'v':	{
					VerboseMode=1;
					break;
				};
				case '3':	{
					ThreeDMode=1;
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
					if ((i+1) < argc) ResFileName=argv[i+1];
					else PrintHelp (argv[0]);
					break;
				};
				case 'n': {
					if ((i+1) < argc) ObjectNum=atoi(argv[i+1]);
					else PrintHelp (argv[0]);
					break;
				};
				case 't': {
					if ((i+1) < argc) Temperature=atof(argv[i+1]);
					else PrintHelp (argv[0]);
					break;
				};
				case 'a': {
					if ((i+1) < argc) LocalizationLength=atof(argv[i+1]);
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
				case 'b':	{
					if ((i+1) < argc) BondDistance = atof(argv[i+1]);
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
            case 'e':	{
              if ((i+1) < argc) StickEnergyDistortion=atof(argv[i+1]);
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
}
