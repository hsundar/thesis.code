
//----------------------------------------------------------------------
// This program is used to distribute certain number of sample
//   points evenly on a template corpus collosum slice image,
//   then map it to an individual corpus collosum image based
//   on a given vector-field data file
//
// Related file:
//   Makefile
//
// Run Command:
//   edis
//   edis -h for help
//----------------------------------------------------------------------


#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream.h>
#include <strstream.h>
#include <fstream.h>
#include <fmclient.h>



#include "RGBAcolorMap.h"

int  iWidth=256, iHeight=256; //slice image size
int  pixVal = 250; //default value of the pixel value of the area of interest
int  wTemp; //window for Template and Target
char * bufRGBA = NULL;
int dx, dy, nZslice;
class RGBAcolorMap *map=NULL;
char infile [200]; //input color map file to be displayed


GLubyte space[] =
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
GLubyte letters[][13] = {
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c, 0x18},
    {0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
    {0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
    {0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce, 0xfc},
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0, 0xff},
    {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xff},
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
    {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e},
    {0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06},
    {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6, 0xc3},
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0},
    {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7, 0xc3},
    {0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3, 0xe3},
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7, 0x7e},
    {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
    {0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66, 0x3c},
    {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
    {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7, 0x7e},
    {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0xff},
    {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
    {0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
    {0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
    {0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
    {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
    {0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03, 0xff}
};

GLuint fontOffset;

void helpInfo();

void makeRasterFont(void)
{
   GLuint i, j;
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   fontOffset = glGenLists (128);
   for (i = 0,j = 'A'; i < 26; i++,j++) {
      glNewList(fontOffset + j, GL_COMPILE);
      glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, letters[i]);
      glEndList();
   }
   glNewList(fontOffset + ' ', GL_COMPILE);
   glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, space);
   glEndList();
}
void init(void)
{
   glShadeModel (GL_FLAT);
   makeRasterFont();
}

void outtext(char *s)
{
   glPushAttrib (GL_LIST_BIT);
   glListBase(fontOffset);
   glCallLists(strlen(s), GL_UNSIGNED_BYTE, (GLubyte *) s);
   glPopAttrib ();
}



//----------------------
void InitEnvironment()
{
  //glEnable(GL_DEPTH_TEST);
  //glDepthFunc(GL_LESS);
  glDisable(GL_CULL_FACE);
  glClearColor(0,0,0,0.0);
  glClearDepth(1.0);
  glShadeModel(GL_SMOOTH);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  //glEnable(GL_LINE_SMOOTH);
  //glEnable(GL_POLYGON_SMOOTH);

  makeRasterFont(); //for initializing font
}
void onIdle(void)
{ 
}
void onSpecialKey(int key,int x,int y)
{
  switch (key){
  case GLUT_KEY_UP:
	nZslice ++;
	if (nZslice > map->getDimZ()-1)
	   nZslice = map->getDimZ()-1;
	break;
  case GLUT_KEY_DOWN:
	if (nZslice > 0) 
	  nZslice --;
	break;
  case GLUT_KEY_LEFT:
	break;
  case GLUT_KEY_RIGHT:
	break;
  }
  glutPostRedisplay();

}
void onReshape(int w, int h)
{

   //setReshapeToView(0,w,h);

   glViewport(0,0,w,h);   

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(0.0,w, 0.0, h, -1.0, 1.0);
   //gluPerspective(45,GLdouble(w)/(float)h,0.5,500);
   //gluLookAt(0,0,300,0,0,0,0,1,0);
   //gluLookAt(0,0,1,0,0,0,0,1,0);

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();

   glClearColor(0,0,0,0.0);
   glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);  
   glColor3f(1.0,1.0,0.0);
}
void onDisplay()
{
  GLfloat ff[4],gg[4];
  glColor3f(0.0,1.0,0.0);

  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
  glDisable(GL_LIGHTING);
  glPushMatrix ();
  //-----------------
    glDrawBuffer(GL_BACK);    
    //glDrawPixels(dx,dy,GL_BLUE,GL_UNSIGNED_BYTE,app->getImage());
    map->getZslice(nZslice, bufRGBA);
    glDrawPixels(dx,dy,GL_RGBA,GL_UNSIGNED_BYTE,bufRGBA);
    //glClear(GL_DEPTH_BUFFER_BIT);

    glColor3f(1.0,1.0,1.0);
    glClear(GL_DEPTH_BUFFER_BIT);
    //glColor3f(1.0,1.0,1.0);
    //glClear(GL_COLOR_BUFFER_BIT);
    glGetFloatv(GL_CURRENT_RASTER_POSITION, ff);
      cout << " --1.(xyzw):"<<ff[0]<<","<<ff[1]<<","<<ff[2]<<","<<ff[3]<<endl;
    glRasterPos2i(0,0);
    glGetFloatv(GL_CURRENT_RASTER_COLOR,gg);
      cout << "   --raster (rgba):"<<gg[0]<<","<<gg[1]<<","<<gg[2]<<","<<gg[3]<<endl;
    outtext("xxx yyy");
    glRasterPos4f(ff[0],ff[1],0,ff[3]);
    glGetFloatv(GL_CURRENT_RASTER_POSITION, ff);
      cout << " --2.(xyzw):"<<ff[0]<<","<<ff[1]<<","<<ff[2]<<","<<ff[3]<<endl;
    //glFlush();

    glutSwapBuffers();
  //-----------------  
  glPopMatrix();
  cout <<" Slice = "<< nZslice<<endl;  
  //fmprstr("xxx yyy");
}

void onKeyStroke(unsigned char key, int x, int y)
{
  switch (key){
  case 27:
     exit(0);
     break;
  case 'h':
     helpInfo();
     break;
  }
}

void onMouseClick(int button, int status, int xx, int yy)
{
   /*
   if (button == GLUT_LEFT_BUTTON && status == GLUT_UP){    
   } else 
   if (button == GLUT_RIGHT_BUTTON && status == GLUT_UP){
   } 
   if (button == GLUT_MIDDLE_BUTTON && status == GLUT_UP){
   } 
   if (button == GLUT_LEFT_BUTTON && status == GLUT_DOWN){    
   } else 
   if (button == GLUT_LEFT_BUTTON && status == GLUT_UP){
   }else 
   if (button == GLUT_RIGHT_BUTTON && status == GLUT_UP){
   }else
   if (button == GLUT_MIDDLE_BUTTON && status == GLUT_UP){
   } */
}

void onMouseMove(int newX,int newY)
{
 
}
void helpInfo()
{
}
void usage(char* runname)
{
  cout <<"\n\nThis program is for displaying RGBA color map slice by slice\n";
  cout <<"  Usage:\n    "<<runname<< " inputFile [ -D dx,dy]\n";
  cout <<"	   inputFile :  the input RGBA colormap file"<<endl;
  cout <<"	   -D<dx,dy> :  dimensions of the volume in X and Y, \n";
  cout <<"                      default: 256x256xZ"<<endl;   
  cout <<"		        !!NOTE: dx and dy must be the power of 2 !!\n";
  cout <<"	   -H: Help information" <<endl;
  cout << endl <<endl;

  cout << " An example:\n";
  cout << "   " << runname <<" Laura.cmp  -D 256,256\n";

  cout << "------------------" <<endl;
  cout << " -by XU, Dongrong " <<endl;
  cout << "     Feb 15, 2002 " <<endl;
  cout << "------------------\n\n" <<endl;
  exit(0);

}
void setDefaultValues()
{
    dx = dy = 256;
}
void getParameters(int argc, char ** argv)
{   
   extern char *optarg;
   int c;

   setDefaultValues();
   if (argc <2){
      usage(argv[0]);
      return;
   } else {
     for (int i=0; i< argc; i++)
        cout <<i<<": "<< argv[i] <<endl;
     strcpy(infile,argv[1]);

     while ((c=getopt(argc-1,argv+1,"D:H")) != -1) {
      //cout << " -- --  char = " << char(c) << endl;
      switch (c) {
      case 'D':
	sscanf(optarg,"%d,%d",&dx,&dy);
	//sscanf(optarg,"%f,%f,%f",&resX,&resY,&resZ);
        break;
      case 'H':
      defualts:
	usage(argv[0]);
	break;
      }
     }
     cout <<"input file = "<<infile<<endl;
     cout << "(dx,dy)=("<<dx<<","<<dy<<")"<<endl;
  }
  //exit(0);
}
  
int main(int argc, char **argv)
{
  char wCaption[240];

  getParameters(argc,argv);
  map = new RGBAcolorMap(infile, dx,dy);
  
  strcpy(wCaption, "Color Map: ");
  strcat(wCaption, infile);

  bufRGBA = new char[4*dx*dy];
  if (bufRGBA == NULL){
	cout <<" memory allocation failed."<<endl;
	exit(0);
  }
  nZslice = map->getDimZ()/2;
  //fminit(); //init font
  //--------------------
  
  glutInit(&argc, argv);
  glutInitWindowSize(dx,dy);
  
  glutInitWindowPosition(200,200);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  wTemp=glutCreateWindow(wCaption); 
  InitEnvironment();
  glutReshapeFunc(onReshape);
  glutKeyboardFunc(onKeyStroke);
  glutSpecialFunc(onSpecialKey);  
  glutMouseFunc(onMouseClick);
  glutMotionFunc(onMouseMove);
  glutDisplayFunc(onDisplay);
  //glutIdleFunc(onIdleSurface);

  glutMainLoop();

  delete []bufRGBA;
  delete map;
  
  return 0;
}


