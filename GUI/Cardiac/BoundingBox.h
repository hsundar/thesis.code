//////////////////////////////////////////////////////////////////////////////
//
//  Author      : Josh Grant
//  Date        : November 12th, 2001
//  File        : BoundingBox.h
//  Description : Header file defining the BoundingBox class
//
//////////////////////////////////////////////////////////////////////////////

#ifndef _BOUNDING_BOX_
#define _BOUNDING_BOX_

#include <Inventor/SbLinear.h>
#include <Inventor/nodes/SoShape.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFVec4f.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/nodes/SoSeparator.h>

//////////////////////////////////////////////////////////////////////////////
//
//  Class: BoundingBox
//
//  A subclass of SoShape used to create an axis aligned wire frame box based
//  on the minBounds and maxBounds fields.  The class also has a field which
//  can be toggled on or off for displaying text coordinate labels at the
//  vertices of the box.
//
//////////////////////////////////////////////////////////////////////////////

class BoundingBox : public SoShape {

  // required by Inventor
  SO_NODE_HEADER(BoundingBox);

public:
  // minimum box coordinates
  SoSFVec3f            minBounds;
  // maximum box coordinates
  SoSFVec3f           maxBounds;
  // true if the coordinates are to be displayed at each vertex
  SoSFBool            textOn;

                      BoundingBox        ();

  // Inventor method, must be called before an object of this type can be
  // created.  
  static void         initClass          ();

protected:
  // render the node
  virtual void        GLRender           (SoGLRenderAction *action);

  // generate primitives, not implemented yet
  virtual void        generatePrimitives (SoAction *action);

  // compute the bounding box
  virtual void        computeBBox        (SoAction *action,
					  SbBox3f &box, SbVec3f &center);

private:
                      ~BoundingBox       ();

  // internal variables
  SoSeparator        *root, *textSep;
  SoCoordinate3      *boxCoords;
  SoIndexedLineSet   *boxLines;

};

#endif /* _BOUNDING_BOX_ */

