#ifndef _GLOBAL_H
#define _GLOBAL_H

#define ZERO 		0.00001
#define PI   		3.1415926
#define ORIGIN          fPoint(0,0,0)
#define MAX_VAL 	65500 

#define UNDEFINED_PT    -5000

#define FIBER_CONT_THRESHOLD 0.125
#define OPACITY_ZERO	0.02f
#define OPAC_MAX	0.6f  //used in VData.C, to be the max value of opacity
#define PIXEL_BYTE	4     //define how many bytes are used for each pixel in the images
#define  noView		160
#define  sampGap	0.5f


//for fixing hole: special identifier color
#define SP_COLOR_R	255
#define SP_COLOR_G      0 
#define SP_COLOR_B      0
#define SP_COLOR_ALP    255

enum {atlasZ=57,atlasX=256,atlasY=256};

enum {FALSE=0, TRUE=1}; //BOOL type
enum {False=0, True=1};
enum {NO=0, YES=1};
enum {ABNORMAL, NORMAL}; // tResult type

struct Pt3d{
  float x,y,z;
};

struct DTensor{
  float a[6];
};

//for class RGBAcolorMap 
struct byteRGBA{
  char rgba[4];
}; 

typedef int BOOL;
typedef int tResult;
typedef unsigned char* pChar;


class Matrix;
class fPoint;
class Vector;
class Plane;
class Line;
class Image;
class Space;
class VolumeData;
class TrackFiber;
class RegularCube;
class DTVolume;
class Ellipsoid;
class VoxMap;

class tNode; 

class RBufferUnit;
class RBuffer;

class ISO;
class MC;
class MCapprox;
class TETRA;

class List;
typedef Space* pSpace;


#define min(x,y) ((x<y)?x:y)
#define max(x,y) ((x>y)?x:y)
//-----------------------


double getSeconds(void);
void exitProgram(char* message);

void out3componentsI(char* msg, int a, int b, int c);// definitions in PDcolorMap.C
void out3componentsF(char* msg, float a, float b, float c);
void out4componentsF(char* msg, float a, float b, float c,float d);

#endif //_GLOBAL_H_
