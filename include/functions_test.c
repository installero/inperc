void TestAllFunctions (); // Runs all function tests
int TestExpression(int expression, char * message); // Returns expression value & prints message if value is 0
int ContextBegin(char * message); // Opens new context level
int ContextEnd(); // Closes last context level

int ContextLevel;

void TestAllFunctions () {

    ContextLevel = 1;
    
    ContextBegin("F_rand() returns random float number (0->1)..."); {
      float f = F_rand();
      TestExpression(f<1, "F_rand should return smth less than 1");
      TestExpression(f>0, "F_rand should return smth greater than 0");
    };
    ContextEnd();
    
    ContextBegin("StickRandomInit(float * Stick) inits stick's parameters with random values..."); {
      float * TestStick;
      TestStick = malloc(ParamsNum*sizeof(float));

      ContextBegin("Checking coordinates..."); {
        ContextBegin("ThreeDMode is off..."); {
          ThreeDMode = 0;
          StickRandomInit(TestStick);
          TestExpression(TestStick[X]>0, "TestStick's x coordinate should be greater than 0");
          TestExpression(TestStick[X]<1, "TestStick's x coordinate should be less than 1");
          TestExpression(TestStick[Y]>0, "TestStick's y coordinate should be greater than 0");
          TestExpression(TestStick[Y]<1, "TestStick's y coordinate should be less than 1");
          TestExpression(TestStick[Z] == 0, "TestStick's z coordinate should be 0");
        }
        ContextEnd();

        ContextBegin("ThreeDMode is on..."); {
          ThreeDMode = 1;
          StickRandomInit(TestStick);
          TestExpression(TestStick[Z]>0, "TestStick's z coordinate should be greater than 0");
          TestExpression(TestStick[Z]>0, "TestStick's z coordinate should be less than 1");
        }
        ContextEnd();
      }
      ContextEnd();

      ContextBegin("Testing stick length"); {
        StickLength = 0;
        StickLengthDistortion = 0;
        StickRandomInit(TestStick);
        TestExpression(TestStick[Length] == 0.0f, "TestStick's length should be equal to 0.0");

        StickLength = 0.01;
        StickLengthDistortion = 0;
        StickRandomInit(TestStick);
        TestExpression(TestStick[Length] == 0.01f, "TestStick's length should be equal to 0.01");

        StickLength = 0.2;
        StickLengthDistortion = 0.05;
        StickRandomInit(TestStick);
        TestExpression(TestStick[Length] > 0.15, "TestStick's length should be greater than 0.15");
        TestExpression(TestStick[Length] < 0.25, "TestStick's length should be less than 0.25");
      }
      ContextEnd();

      ContextBegin("Testing stick width"); {
        StickWidth = 0;
        StickWidthDistortion = 0;
        StickRandomInit(TestStick);
        TestExpression(TestStick[Width] == 0.0f, "TestStick's width should be equal to 0.0");

        StickWidth = 0.01;
        StickWidthDistortion = 0;
        StickRandomInit(TestStick);
        TestExpression(TestStick[Width] == 0.01f, "TestStick's width should be equal to 0.01");

        StickWidth = 0.2;
        StickWidthDistortion = 0.05;
        StickRandomInit(TestStick);
        TestExpression(TestStick[Width] > 0.15, "TestStick's width should be greater than 0.15");
        TestExpression(TestStick[Width] < 0.25, "TestStick's width should be less than 0.25");
      }
      ContextEnd();

      ContextBegin("Testing stick angles"); {
        ContextBegin("ThreeDMode is off..."); {
          StickFiDistortion = 0;
          StickRandomInit(TestStick);
          TestExpression(TestStick[Fi] == 0.0f, "TestStick's fi should be equal to 0.0");

          StickFiDistortion = pi/2;
          StickRandomInit(TestStick);
          TestExpression(TestStick[Fi] > -pi/2, "TestStick's fi should be greater than -pi/2");
          TestExpression(TestStick[Fi] < pi/2, "TestStick's fi should be less than pi/2");
        };
        ContextEnd();

        ContextBegin("ThreeDMode is on..."); {
          ThreeDMode = 1;
          StickThetaDistortion = 0.0f;
          StickRandomInit(TestStick);
          TestExpression(TestStick[Theta] == pi/2, "TestStick's theta should be equal to pi/2");

          StickThetaDistortion = pi/2;
          StickRandomInit(TestStick);
          TestExpression(TestStick[Theta] > 0.0f, "TestStick's theta should be greater than 0.0");
          TestExpression(TestStick[Theta] < pi, "TestStick's theta should be less than pi");
        };
        ContextEnd();
      }
      ContextEnd();

      ContextBegin("Testing stick energy"); {
        StickEnergyDistortion = 0;
        StickRandomInit(TestStick);
        TestExpression(TestStick[Energy] == 0.0f, "TestStick's energy should be equal to 0.0");

        StickEnergyDistortion = 10;
        StickRandomInit(TestStick);
        TestExpression(TestStick[Energy] > 0, "TestStick's energy should be greater than 0");
        TestExpression(TestStick[Energy] < 10, "TestStick's energy should be less than 10");
      }
      ContextEnd();

      free(TestStick);
    };
    ContextEnd();

    ContextBegin("BondCriteria simply returns equation for Eta"); {
      LocalizationLength = 0.5;
      Temperature = 10;
      TestExpression(BondCriteria(0,0) == 0, "should return 0 if both dE & dR are 0s");
      TestExpression(BondCriteria(1,0) == 4.0f, "should return 2dR/a if dE is 0 and dR is not");
      TestExpression(BondCriteria(0,2) == 0.2f, "should return dE/T if dR is 0 and dE is not");
      TestExpression(BondCriteria(5,7) == 20.7f, "should return 2Dr/a + dE/T if both dR and dE are not 0s");
    };
    ContextEnd();

    ContextBegin("StickToPointDistance finds distance between a stick and a point"); {
      float * TestStick;
      float * TestPoint;
      TestStick = malloc(ParamsNum*sizeof(float));
      TestPoint = malloc(ParamsNum*sizeof(float));

      TestPoint[X] = TestPoint[Y] = TestPoint[Z] = TestPoint[Length] = TestPoint[Width] = TestPoint[Fi] = TestPoint[Theta] = TestPoint[Energy] = 0;
      TestStick[X] = TestStick[Y] = TestStick[Z] = TestStick[Length] = TestStick[Width] = TestStick[Fi] = TestStick[Theta] = TestStick[Energy] = 0;

      TestExpression(StickToPointDistance(TestStick, TestPoint) == 0.0f, "should return 0 if both stick & point are zero lenght points at {0,0,0}");

      TestStick[X] = 0.3;
      TestExpression(StickToPointDistance(TestStick, TestPoint) == 0.3f, "should return 0.3 if sticks X is 0.3, and length is 0");

      TestStick[Y] = 0.4;
      TestExpression(StickToPointDistance(TestStick, TestPoint) == 0.5f, "should return 0.5 if sticks X is 0.3, Y is 0.4, and length is 0");

      TestStick[X] = 0.03;
      TestStick[Y] = 0.04;
      TestStick[Z] = 0.12;
      TestExpression(StickToPointDistance(TestStick, TestPoint) == 0.13f, "should return 0.13 if sticks X=0.03, Y=0.04, Z=0.12, L=0");

      TestStick[X] = 0.0f;
      TestStick[Y] = 0.0f;
      TestStick[Z] = 0.25f;
      TestStick[Length] = 0.10f;
      TestExpression(StickToPointDistance(TestStick, TestPoint) == 0.2f, "should return 0.1 if sticks X=0, Y=0, Z=0.25, L=0.1 (Theta,Fi=0)");

      TestStick[Theta] = pi/2;
      TestExpression(StickToPointDistance(TestStick, TestPoint) == 0.25f, "should return 0.25 if sticks X=0, Y=0, Z=0.25, L=0.1 (Theta=pi/2,Fi=0)");

      TestStick[X] = 0.3f;
      TestStick[Y] = 0.2f;
      TestStick[Z] = 0.1f;
      TestStick[Theta] = pi/4;
      TestStick[Fi] = pi/2;
      TestStick[Length] = 1.5f;
      TestExpression(
              StickToPointDistance(TestStick, TestPoint) > 0.30822f && StickToPointDistance(TestStick, TestPoint) < 0.30823f,
              "should return approx 0.308 for stick {0.3,0.2,0.1}, L=0.15, Theta=pi/4,Fi=pi/2");

      free(TestStick);
      free(TestPoint);
    };
    ContextEnd();

    ContextBegin("StickToBoundaryDistance finds distance from a stick to a boundary at a certain direction"); {
      float * TestStick;
      TestStick = malloc(ParamsNum*sizeof(float));
      TestStick[X] = 0.6;
      TestStick[Y] = 0.3;
      TestStick[Z] = 0.4;
      TestStick[Length] = 0.1;
      TestStick[Width] = TestStick[Energy] = 0;
      TestStick[Fi] = pi/3;
      TestStick[Theta] = pi/6;

      TestExpression(
              StickToBoundaryDistance(TestStick, MIN_X) > 0.5875 && StickToBoundaryDistance(TestStick, MIN_X) < 0.5876,
              "should return approx 0.5875 for stick {0.5,0.3,0.6}, L=0.1, Theta=pi/6,Fi=pi/3 and MIN_X border");
      TestExpression(
              StickToBoundaryDistance(TestStick, MAX_X) > 0.3874 && StickToBoundaryDistance(TestStick, MAX_X) < 0.3875,
              "should return approx 0.3875 for stick {0.5,0.3,0.6}, L=0.1, Theta=pi/6,Fi=pi/3 and MAX_X border");
      TestExpression(
              StickToBoundaryDistance(TestStick, MIN_Y) > 0.27834 && StickToBoundaryDistance(TestStick, MIN_Y) < 0.27835,
              "should return approx 0.27834 for stick {0.5,0.3,0.6}, L=0.1, Theta=pi/6,Fi=pi/3 and MIN_Y border");
      TestExpression(
              StickToBoundaryDistance(TestStick, MAX_Y) > 0.67834 && StickToBoundaryDistance(TestStick, MAX_Y) < 0.67835,
              "should return approx 0.67834 for stick {0.5,0.3,0.6}, L=0.1, Theta=pi/6,Fi=pi/3 and MAX_Y border");
      TestExpression(
              StickToBoundaryDistance(TestStick, MIN_Z) > 0.3566 && StickToBoundaryDistance(TestStick, MIN_Z) < 0.3567,
              "should return approx 0.3566 for stick {0.5,0.3,0.6}, L=0.1, Theta=pi/6,Fi=pi/3 and MIN_Z border");
      TestExpression(
              StickToBoundaryDistance(TestStick, MAX_Z) > 0.5566 && StickToBoundaryDistance(TestStick, MAX_Z) < 0.5567,
              "should return approx 0.5566 for stick {0.5,0.3,0.6}, L=0.1, Theta=pi/6,Fi=pi/3 and MAX_Z border");
    };
    ContextEnd();

    ContextBegin("Wild & Rough StickToStickDistance finds distance between two sticks wildly & roughly"); {
      float * TestStick1;
      float * TestStick2;
      TestStick1 = malloc(ParamsNum*sizeof(float));
      TestStick2 = malloc(ParamsNum*sizeof(float));
      TestStick1[Fi] = TestStick1[Theta] = TestStick2[Fi] = TestStick2[Theta] = 0;
      TestStick1[X] = 0.6;
      TestStick1[Y] = 0.5;
      TestStick1[Z] = 0.4;
      TestStick2[X] = 0.1;
      TestStick2[Y] = 0.2;
      TestStick2[Z] = 0.3;
      TestStick1[Length] = 0.1;
      TestStick2[Length] = 0.3;

      TestExpression(
              WildStickToStickDistance(TestStick1, TestStick2) == 0.3f,
              "wildSTSD should return 0.3 for sticks {0.6,0.5,0.3} L1 = 0.1 and {0.1,0.2,0.3} L = 0.3");

      TestExpression(
              RoughStickToStickDistance(TestStick1, TestStick2) > 0.3916f && RoughStickToStickDistance(TestStick1, TestStick2) < 0.3917f,
              "roughSTSD should return appr 0.31961 for sticks {0.6,0.5,0.3} L1 = 0.1 and {0.1,0.2,0.3} L = 0.3");
    }
    ContextEnd();

    ContextBegin("StickToStickDistance finds distance between two sticks"); {
      float * TestStick1;
      float * TestStick2;
      TestStick1 = malloc(ParamsNum*sizeof(float));
      TestStick2 = malloc(ParamsNum*sizeof(float));
      TestStick1[X] = TestStick1[Y] = TestStick1[Z] = 0;
      TestStick2[X] = TestStick2[Y] = TestStick2[Z] = 0.5;
      TestStick1[Width] = TestStick2[Width] = TestStick1[Fi] = TestStick2[Fi] = TestStick1[Theta] = TestStick2[Theta] = 0;
      TestStick1[Length] = TestStick2[Length] = 1;

      TestExpression(
              StickToStickDistance(TestStick1, TestStick2) > 0.7071 && StickToStickDistance(TestStick1, TestStick2) < 0.7072,
              "should return approx 0.7071 for {0,0,0,1,0,0} & {1/2,1/2,1/2,1,0,0}");

      TestStick2[Z] = 0;

      TestExpression(
              StickToStickDistance(TestStick1, TestStick2) > 0.7071 && StickToStickDistance(TestStick1, TestStick2) < 0.7072,
              "should return approx 0.7071 for {0,0,0,1,0,0} & {1/2,1/2,0,1,0,0}");

      TestStick1[Theta] = TestStick2[Theta] = pi/2;

      TestExpression(
              StickToStickDistance(TestStick1, TestStick2) == 0.5,
              "should return 0.5 for {0,0,0,1,0,pi/2} & {1/2,1/2,0,1,0,pi/2}");

      TestStick1[Fi] = TestStick2[Fi] = -pi/4;

      TestExpression(
              StickToStickDistance(TestStick1, TestStick2) > 0.7071 && StickToStickDistance(TestStick1, TestStick2) < 0.7072,
              "should return approx 0.7071 for {0,0,0,1,-pi/4,pi/2} & {0.5,0.5,0,1,-pi/4,pi/2}");

      TestStick2[Z] = 0.5;

      TestExpression(
              StickToStickDistance(TestStick1, TestStick2) > 0.866025 && StickToStickDistance(TestStick1, TestStick2) < 0.866026,
              "should return approx 0.866025 for {0,0,0,1,-pi/4,pi/2} & {0.5,0.5,0.5,1,-pi/4,pi/2}");

      TestStick2[Theta] = 0;

      TestExpression(
              StickToStickDistance(TestStick1, TestStick2) > 0.7071 && StickToStickDistance(TestStick1, TestStick2) < 0.7072,
              "should return approx 0.7071 for {0,0,0,1,-pi/4,pi/2} & {0.5,0.5,0.5,1,-pi/4,0}");

      TestStick2[Fi] = pi/4;

      TestExpression(
              StickToStickDistance(TestStick1, TestStick2) > 0.7071 && StickToStickDistance(TestStick1, TestStick2) < 0.7072,
              "should return approx 0.7071 for {0,0,0,1,-pi/4,pi/2} & {0.5,0.5,0.5,1,pi/4,0}");

      TestStick2[Theta] = -pi/4;

      TestExpression(
              StickToStickDistance(TestStick1, TestStick2) > 0.8535 && StickToStickDistance(TestStick1, TestStick2) < 0.8536,
              "should return approx 0.85355 for {0,0,0,1,-pi/4,pi/2} & {0.5,0.5,0.5,1,pi/4,-pi/4}");

      TestStick2[Theta] = pi/4;

      TestExpression(
              StickToStickDistance(TestStick1, TestStick2) > 0.38268 && StickToStickDistance(TestStick1, TestStick2) < 0.38269,
              "should return approx 0.38268 for {0,0,0,1,-pi/4,pi/2} & {0.5,0.5,0.5,1,pi/4,pi/4}");
    };
    ContextEnd();

    /*ContextBegin(); {*/
    /*};*/
    /*ContextEnd();*/

};

int TestExpression (int expression, char * message){
    int i;
    for(i = 0; i < ContextLevel; i++) fprintf(stderr, "\t");
    int color = expression ? 32 : 31;
    fprintf(stderr, "%c[%d;%d;%dm%s...", 0x1B, 1, color, 40, message);
    fprintf(stderr, "%c[%dm", 0x1B, 0);
    fprintf(stderr, "\n");
    return expression;
};

int ContextBegin (char * message){
    int i;
    for(i = 0; i < ContextLevel; i++) fprintf(stderr, "\t");
    fprintf(stderr, "%s\n", message);
    ContextLevel++;
    return 0;    
}

int ContextEnd (){
    ContextLevel--;
    int i;
    for(i = 0; i < ContextLevel; i++) fprintf(stderr, "\t");
    fprintf(stderr, "Done\n");
    return 0;    
};
