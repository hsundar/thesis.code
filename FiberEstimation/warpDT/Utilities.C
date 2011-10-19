#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

#include "Global.h"


void printTensor(char *explanation, DTensor dt)
{
  int i;
  cout <<explanation;
  for (i=0; i< 6; i++)
    cout <<"  "<<dt.a[i]<<"  ";
  cout << endl;
}
