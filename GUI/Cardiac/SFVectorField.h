//////////////////////////////////////////////////////////////////////////////
//
//  Author      : Hari Sundar
//  Date        : August 23rd, 2006
//  File        : SFVectorField.h
//  Description : Header file defining the SFVectorField class
//
//////////////////////////////////////////////////////////////////////////////

#ifndef _SF_VECTOR_FIELD_
#define _SF_VECTOR_FIELD_

#include <Inventor/fields/SoSubField.h>
#include <Inventor/SbLinear.h>

class SFVectorField : public SoSField {

  // Inventor defined macros for subclassing SoSField
  SO_SFIELD_REQUIRED_HEADER(SFVectorField);
  SO_SFIELD_CONSTRUCTOR_HEADER(SFVectorField);

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

  int operator == (const SFVectorField &d) const;
  int operator != (const SFVectorField &d) const
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


#endif /* _SF_VECTOR_FIELD_ */

