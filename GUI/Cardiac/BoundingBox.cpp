//////////////////////////////////////////////////////////////////////////////
//
//  Author      : Josh Grant
//  Date        : November 12th, 2001
//  File        : BoundingBox.cpp
//  Description : Implementation of the BoundingBox class
//
//////////////////////////////////////////////////////////////////////////////

#include <Inventor/SbBox.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/bundles/SoMaterialBundle.h>
#include <Inventor/elements/SoGLTextureCoordinateElement.h>
#include <Inventor/elements/SoGLTextureEnabledElement.h>
#include <Inventor/elements/SoLightModelElement.h>
#include <Inventor/elements/SoMaterialBindingElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/misc/SoState.h>
#include <cstring>

#include <iostream>
//#include <config.h>

#include "BoundingBox.h"

//////////////////////////////////////////////////////////////////////////////
//
// BoundingBox class
//
//////////////////////////////////////////////////////////////////////////////

// Inventor macro
SO_NODE_SOURCE(BoundingBox);

// vertices used to create a cube
static const int32_t bbCubeVerts[8][3] =
{
  {0, 0, 0},
  {1, 0, 0},
  {1, 1, 0},
  {0, 1, 0},
  {0, 0, 1},
  {1, 0, 1},
  {1, 1, 1},
  {0, 1, 1}
};

// indexes used to create the edges, -1 is used to deliminate between edges
static const int32_t bbCubeEdges[36] =
{
  0,1,-1, 1,2,-1, 2,3,-1, 3,0,-1,
  4,5,-1, 5,6,-1, 6,7,-1, 7,4,-1,
  0,4,-1, 1,5,-1, 2,6,-1, 3,7,-1
};

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Inventor method which must be called before an object of this type can
//    be created.
//
//////////////////////////////////////////////////////////////////////////////
void BoundingBox::initClass ()
{
  SO_NODE_INIT_CLASS(BoundingBox, SoShape, "Shape");
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Create the fields and then make the nodes required by this node.
//
//////////////////////////////////////////////////////////////////////////////
BoundingBox::BoundingBox ()
{
  // Inventor macro
  SO_NODE_CONSTRUCTOR(BoundingBox);

  // Inventor macro used to register a field with the database
  SO_NODE_ADD_FIELD(minBounds,  (-1.0, -1.0, -1.0));
  SO_NODE_ADD_FIELD(maxBounds,  ( 1.0,  1.0,  1.0));
  SO_NODE_ADD_FIELD(textOn,                 (TRUE));

  // the root of the scene for this node
  root = new SoSeparator();

  // create the box separator
  SoSeparator *boxSep = new SoSeparator();

  // the coordinates for the box
  boxCoords = new SoCoordinate3();
  boxCoords->point.setNum(8);
  boxSep->addChild(boxCoords);

  // the lines of the box
  boxLines  = new SoIndexedLineSet();
  boxLines->coordIndex.setNum(36);
  // setting indexes for each edge of the cube
  boxLines->coordIndex.setValues(0, 36, bbCubeEdges);
  boxSep->addChild(boxLines);
  
  root->addChild(boxSep);

  // create the text nodes, including a transform for each vertice offset
  textSep = new SoSeparator();
  for (int i = 0; i < 8; i++) {
    SoSeparator *temp = new SoSeparator();
    SoTransform *trans = new SoTransform();
    temp->addChild(trans);
    SoText2 *text = new SoText2();
    text->justification.setValue(SoText2::CENTER);
    temp->addChild(text);
    textSep->addChild(temp);
  }
  root->addChild(textSep);

  // make sure tree doesn't get deleted
  root->ref();

}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Delete the scene within the node
//
//////////////////////////////////////////////////////////////////////////////
BoundingBox::~BoundingBox ()
{
  root->unref();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Render the box and text
//
//////////////////////////////////////////////////////////////////////////////
void
BoundingBox::GLRender (SoGLRenderAction *action)
{
  SbVec3f corner[2], ctr, offset, *vptr;
  int i, j;
  bool text;
  char str[50], buf[10];

  // grab the current state
  SoState *state = action->getState();

  if (!shouldGLRender(action))
    return;

  // get the latest values from the fields
  corner[0] = minBounds.getValue();
  corner[1] = maxBounds.getValue();
  text      = textOn.getValue();

  // set the coordinates for the LineSet to point to
  vptr = boxCoords->point.startEditing();
  for (i = 0; i < 8; i++)
    for (j = 0; j < 3; j++)
      vptr[i][j] = corner[bbCubeVerts[i][j]][j];

  // a hack to prevent the SoText2 node from bombing on Linux
#ifdef NO_SOTEXT2
  text = false;
#endif

  // if text is true then set the text nodes
  if (text) {
	ctr = (corner[1] - corner[0]);
	ctr /= 2.0;
    for (i = 0; i < 8; i++) {
      // create the string for the text
      strcpy(str, "(");
      sprintf(buf, "%6.2f", vptr[i][0]);
      strcat(str, buf);
      strcat(str, ",");
      sprintf(buf, "%6.2f", vptr[i][1]);
      strcat(str, buf);
      strcat(str, ",");
      sprintf(buf, "%6.2f", vptr[i][2]);
      strcat(str, buf);
      strcat(str, ")");

      SoSeparator *sep   = (SoSeparator *)textSep->getChild(i);
      SoTransform *trans = (SoTransform *)sep->getChild(0);

      offset = vptr[i] - ctr;
      //offset.normalize();
      trans->translation.setValue(vptr[i].getValue());// + offset * 0.1);
      SoText2 *t         = (SoText2 *)sep->getChild(1);
      t->string.setValue(str);
      //printf("%s\n", str);
    }

    textSep->ref();
    if (root->findChild(textSep) < 0)
      root->addChild(textSep);
  }
  // otherwise unreference the text nodes
  else {
    textSep->unrefNoDelete();
    if (root->findChild(textSep) >= 0)
      root->removeChild(textSep);
  }
  
  boxCoords->point.finishEditing();

  // render the node
  root->GLRender(action);
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Not implemented yet
//
//////////////////////////////////////////////////////////////////////////////
void
BoundingBox::generatePrimitives (SoAction *action)
{
  
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Set the bounding box of the scene
//
//////////////////////////////////////////////////////////////////////////////
void
BoundingBox::computeBBox (SoAction *action,
		          SbBox3f &box, SbVec3f &center)
{
	center = (minBounds.getValue() + maxBounds.getValue());
	center /= 2.0;
  box.setBounds(minBounds.getValue(), maxBounds.getValue());
}


  
