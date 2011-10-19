//////////////////////////////////////////////////////////////////////////////
//
//  Author      : Josh Grant
//  Date        : December 16th, 2001
//  File        : SFScalarField.h
//  Description : Header file defining the SFScalarField class
//
//////////////////////////////////////////////////////////////////////////////

#ifndef _SF_SCALAR_FIELD_
#define _SF_SCALAR_FIELD_

#include <Inventor/fields/SoSubField.h>
#include <Inventor/SbLinear.h>

//////////////////////////////////////////////////////////////////////////////
//
//  Class: SFScalarField
//
//  Base class for all data fields.  Even though this class can be
//  instantiated it was meant to only serve as an abstract class.  It offers
//  the ability to add a data field to an Inventor node or engine so that the
//  Database detects when it has been updated.  It is very similar to the
//  SoSFImage field, except for one major difference.  It also has an option in
//  the setValue() function for copying the data or not.  This is very useful
//  when setting a field from an Engine.  It can greatly improve the
//  performance time of interactive tools.
//
//  The copy option works by either copying all the data or just saving the
//  pointer to the data.  Of course this will only work if the memory the
//  pointer references isn't deleted.
//
//////////////////////////////////////////////////////////////////////////////

class SFScalarField : public SoSField {

  // Inventor defined macros for subclassing SoSField
  SO_SFIELD_REQUIRED_HEADER(SFScalarField);
  SO_SFIELD_CONSTRUCTOR_HEADER(SFScalarField);

public:
  static void  initClass ();

  // set the dimensions, number of dimensions, number of data components and
  // the pointer to the data, triggers a valueChanged() to Inventor
  void         setValue (const int dm[3], float *ptr, bool copyData=true);
  // get the data value, triggers an evaluate() to Inventor
  float *      getValue (int dm[3]) const;

  // the following functions do not trigger any Inventor calls
  float *      getDataPtr       () const {return data;}
  const int *  getDims          () const {return dims;}
  int          getNumValues     () const {return numValues;}
  bool         getDataCopied    () const {return dataCopied;}
  void         getMinMax        (float &min, float &max) const;
  void         allocateSpace    (const int dm[3]);

  int operator == (const SFScalarField &d) const;
  int operator != (const SFScalarField &d) const
                    { return !((*this) == d); }

protected:
  // used to maintain the dimensions of the data.
  int     dims[3];
  // these variables are here just to keep the code general, currently they
  // are just constants
  int     numComponents;
  int     numDimensions;

  // current number of values being managed, this is the product of all the
  // dimension values times the number of components.
  int     numValues;
  // pointer to data
  float   *data;

private:
  // used to determine if the data was copied over to this field or if this
  // field just points to it somewhere else
  bool    dataCopied;
  
  // whenever the data is copied space is allocated.  this variable is used
  // to keep track of how many values have been allocated.  This way if data
  // of the same or smaller size needs to be copied, we don't need to
  // allocate more memory.
  int     valuesAllocated;
  
  // pointer to allocated array.
  float   *myData;

  virtual SbBool readValue  (SoInput *in);
  virtual void   writeValue (SoOutput *out) const;

};


#endif /* _SF_SCALAR_FIELD_ */

