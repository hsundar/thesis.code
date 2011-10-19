//////////////////////////////////////////////////////////////////////////////
//
//  Author      : Hari Sundar
//  Date        : August 23rd, 2006
//  File        : SFVectorField.cpp
//  Description : Implementation of SFVectorField class
//
//////////////////////////////////////////////////////////////////////////////

#include "SFVectorField.h"
#include <cstring>

// The below block of macros is needed if the code is going to be portable
// across Coin and Inventor.  If you don't care about this portability
// problem then either comment out the conditionals.  Also if you are no
// longer using the config.h file then comment it out as well.  It is only
// referenced in order to see if HAVE_COIN is defined.
// #include <config.h>
#define HAVE_COIN 1
// Required Inventor macros for subclassing SoSField
#if HAVE_COIN
// This means we are using Coin and these macros must be referenced
PRIVATE_TYPEID_SOURCE(SFVectorField);
PRIVATE_EQUALITY_SOURCE(SFVectorField);
#else
// This means we are using Inventor and these macros must be referenced
SO__FIELD_ID_SOURCE(SFVectorField);
SO__FIELD_EQ_SAME_SOURCE(SFVectorField);
#endif

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Default constructor which sets all variables to zero.
//
//////////////////////////////////////////////////////////////////////////////
SFVectorField::SFVectorField ()
{
  dims[0] = dims[1] = dims[2] = 0;
  data = myData   = NULL;
  
  // these are set to constants, but the rest of the code is set up to
  // generalize to n dimensions and m components
  numDimensions   = 3;
  numComponents   = 3;
  
  numValues       = 0;
  valuesAllocated = 0;
  dataCopied      = false;
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    deletes all allocated arrays
//
//////////////////////////////////////////////////////////////////////////////
SFVectorField::~SFVectorField ()
{
  if (myData)
    delete [] myData;
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    needed for Inventor to recognize the class.  This must be called at
//    startup time and is best to do it immediately following the
//    initialization of Inventor.
//
//////////////////////////////////////////////////////////////////////////////
void
SFVectorField::initClass ()
{
  classTypeId = SoType::createType(SoSField::getClassTypeId(),
				   "SFVectorField",
				   &SFVectorField::createInstance);
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    set the value of the field.  This will set the dimensions and the data
//    itself.  If the copyData argument is false then only a 
//    pointer to the data is saved, otherwise all the data is copied over to
//    this field.  (Note: copying the data can slow performance when dealing
//    with large arrays.  If the data is being passed from an Engine it is
//    best to set copyData to false.)
//
//////////////////////////////////////////////////////////////////////////////
void
SFVectorField::setValue (const int dm[3], float *ptr, bool copyData)
{
  int i, nextValuesAlloc;

  dims[0] = dm[0];
  dims[1] = dm[1];
  dims[2] = dm[2];

  // compute how many bytes the data needs
  nextValuesAlloc = numComponents;
  for (i = 0; i < numDimensions; i++)
    nextValuesAlloc *= dims[i];

  // if we are not to copy the data then just save the pointer
  if (!copyData) {
    data       = ptr;
    dataCopied = false;
  }
  // otherwise allocate memory if needed, and copy the data
  else {

    // check to see if we need more memory
    if (nextValuesAlloc > valuesAllocated) {
      valuesAllocated = nextValuesAlloc;

      if (myData)
	delete [] myData;

      myData = new float[nextValuesAlloc];
    }

    // copy the data to the local array
    memcpy(myData, ptr, sizeof(float)*nextValuesAlloc);

    data = myData;

    dataCopied = true;
  }

  // set the number of values currently being used
  numValues = nextValuesAlloc;

  // notify database of changes
  valueChanged();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    return the data and notify Inventor that there is a request for this
//    field.
//
//////////////////////////////////////////////////////////////////////////////
float *
SFVectorField::getValue (int dm[3]) const
{
  // Inventor function
  evaluate();

  dm[0] = dims[0];
  dm[1] = dims[1];
  dm[2] = dims[2];

  return data;
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    return the minimum and maximum value of the data.  This is to be
//    overloaded by subclasses of this class.
//
//////////////////////////////////////////////////////////////////////////////
void
SFVectorField::getMinMax (float &min, float &max) const
{
  if (numValues > 0) {
    min = data[0];
    max = data[0];
  }
  
  for (int i = 1; i < numValues; i++) {
    if (data[i] < min)
      min = data[i];
    else if (data[i] > max)
      max = data[i];
  }
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    allocates enough memory needed to support the input dimensions.
//
//////////////////////////////////////////////////////////////////////////////
void
SFVectorField::allocateSpace (const int dm[3])
{
  int i, nextValuesAlloc;

  int nd = numDimensions;
  int nc = numComponents;

  // compute how many bytes the data needs
  nextValuesAlloc = nc;
  for (i = 0; i < nd; i++)
    nextValuesAlloc *= dm[i];

  // check to see if we need more memory
  if (nextValuesAlloc > valuesAllocated) {
    valuesAllocated = nextValuesAlloc;

    if (myData)
      delete [] myData;

    myData = new float[nextValuesAlloc];
  }
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    assignment operator.  If the data from 'd' was copied then it is copied
//    to this field
//
//////////////////////////////////////////////////////////////////////////////
const SFVectorField &
SFVectorField::operator = (const SFVectorField &d)
{
  int dm[3];
  float *ptr = d.getValue(dm);

  setValue(dm, ptr, d.getDataCopied());

  return *this;
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    equality operator.
//
//////////////////////////////////////////////////////////////////////////////
int
SFVectorField::operator == (const SFVectorField &d) const
{
  // check the easy ones first
  if (numComponents != d.numComponents ||
      numDimensions != d.numDimensions ||
      numValues     != d.numValues     ||
      dataCopied    != d.dataCopied)
    return false;

  // make sure the dimension values are the same
  for (int i = 0; i < numDimensions; i++)
    if (dims[i] != d.dims[i])
      return false;

  if (!dataCopied && data != d.data)
    return false;
  else if (dataCopied && memcmp(data, d.data, sizeof(float)*numValues))
    return false;

  return true;
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    not yet implemented.
//
//////////////////////////////////////////////////////////////////////////////
SbBool
SFVectorField::readValue (SoInput *in)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    not yet implemented.
//
//////////////////////////////////////////////////////////////////////////////
void
SFVectorField::writeValue (SoOutput *out) const
{
  ;
}
