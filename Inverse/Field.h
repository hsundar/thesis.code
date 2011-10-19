/**
*  @file	Field.h
*  @brief	Container class to hold a 3D vector image, like velocity or displacement. 
*  @author	Hari Sundar
*  @date	1/10/06
*
*/

#ifndef _FIELD_H_
#define _FIELD_H_

#include <string>
#include "Point.h"

/**
*  @brief	Container class to hold a 3D vector image,  like velocity, displacement 
*		or a tensor. This is the parent class for all derived field types.
* 		
*		The class is templated on the datatype of the vector component, and the 
*		number of components. The class	assumes that double and float will be the 
*		common types.
*  @author	Hari Sundar
*  @date	1/10/06, 8/6/06
* 
* Container class to hold a 3D vector image,  like velocity, displacement or a tensor. 
* This is the parent class for all derived field types.		
*
* The class is templated on the datatype of the vector component, and the  number of 
* components. The class	assumes that double and float will be the common type. The data 
* is stored in a single C-array of dimension n*x*y*z.
* 
**/

template<typename T, int n>
class Field {
  public:

    /** @name Constructors and destructors **/
    //@{
    Field();
    Field(const Field<T, n>& rhs);
    virtual ~Field();
    //@}

    /**
     *  @brief		initizlizes the Field, from disk 
     *  @param		fileName string containing the Field header file in the MetaIO
     * 				format.
     *  @return		true if successfull, false otherwise
    **/
    virtual bool init(std::string fileName);

    /**
     *  @brief		initizlizes the Field 
     *  @param		sizeX the dimension of the field along the X axis
     *  @param		sizeY the dimension of the field along the Y axis
     *  @param		sizeZ the dimension of the field along the Z axis 
     *  @return		true if successfull, false otherwise
     *
     *  Follow by reading in the raw field file if required.
    **/
    bool init(int sizeX, int sizeY, int sizeZ);

    /**
     *  @brief		return the size of the Field in the specified direction.
     *  @param		dim	Dimension, must be 0, 1, or 2.
     *  @return		dimension size.
     *
    **/
    int getSize( int dim ) const { return m_Dims[dim]; };

    /**
     *  @brief		return the size of a voxel in the specified direction.
     *  @param		dim	Dimension, must be 0, 1, or 2.
     *  @return		voxel size in specified direction.
     *
    **/
    double getVoxelSize( int dim ) const { return m_VoxelSize[dim]; };

    /**
     *  @brief		set the size of a voxel in the specified direction.
     *  @param		dim	Dimension, must be 0, 1, or 2.
     *  @param		voxelSize	the desired voxel size.	
    **/
    void setVoxelSize( int dim, double voxelSize ) { m_VoxelSize[dim] = voxelSize; };

    /**
     *  @brief		returns the file name of the raw volume.
     *  @return		CString for the filename.
    **/
    const std::string& getRawFileName() const { return m_RawFileName; };

    /**
     *  @brief		return pointer to data.
     *  @return		pointer to the field data.
     *
    **/
    T* getDataPtr() const { return m_Data; };

    /**
     *  @brief		return pointer to data at voxel location.
     *  @param		x,y,z the voxel location.
     *  @return		pointer to the field data.
     *
    **/
    T* getAt(int x, int y, int z) const { 
		if ( (x < m_Dims[0]) && (y < m_Dims[1]) && (z < m_Dims[2]) )  
			return m_Data + n*(m_Dims[0]*(z*m_Dims[1] + y) +x) ;
		else 
			return static_cast<T*>(0);
	}
	
	Point getPt(int x, int y, int z) const {
		if ( (x>=0) && (y>=0) &&(z>=0) && (x < m_Dims[0]) && (y < m_Dims[1]) && (z < m_Dims[2]) )  {
			unsigned int idx = n*(m_Dims[0]*(z*m_Dims[1] + y) + x);
			return Point(m_Data[idx], m_Data[idx+1], m_Data[idx+2]) ;
		} else 
			return Point();
	}
	
	Point getAt(double x, double y, double z) const { 
		if ( (x < 0) || (y < 0) || (z < 0) ) {
			return Point();
		}
		// Function takes a float index and returns velocity at that point after
		// interpolation. Currently, linear interpolation is used. Might need to
		// change that eventually.
		int i, j, k;
		double dx, dy, dz;
		i = (int)x; j = (int)y; k = (int)z;
		dx = x- i;  dy = y - j;  dz = z - k;

		// Simple tri-linear Interpolation ....
		
		Point p1,p2,p3,p4,p5,p6,p7,p8;
		double curv;
		
		/*
		if ( fabs(curv) < 0.8) p2=Point();
		if ( fabs(curv) < 0.8) p3=Point();
		if ( fabs(curv) < 0.8) p4=Point();
		if ( fabs(curv) < 0.8) p5=Point();
		if ( fabs(curv) < 0.8) p6=Point();
		if ( fabs(curv) < 0.8) p7=Point();
		if ( fabs(curv) < 0.8) p8=Point();
		*/
		
		p1 = getPt( i,j,k ); p1.normalize();
		p2 = getPt( i+1,j,k ); p2.normalize(); 			curv = p2.dot(p1); 	if (curv < 0) p2 = -p2; 		 
		p3 = getPt( i,j+1,k ); p3.normalize(); 			curv = p3.dot(p1); 	if (curv < 0) p3 = -p3; 		
		p4 = getPt( i,j,k+1 ); p4.normalize(); 			curv = p4.dot(p1); 	if (curv < 0) p4 = -p4; 		
		p5 = getPt( i+1,j,k+1 ); p5.normalize(); 		curv = p5.dot(p1); 	if (curv < 0) p5 = -p5; 
		p6 = getPt( i,j+1,k+1 ); p6.normalize(); 		curv = p6.dot(p1); 	if (curv < 0) p6 = -p6; 
		p7 = getPt( i+1,j+1,k ); p7.normalize(); 		curv = p7.dot(p1); 	if (curv < 0) p7 = -p7; 
		p8 = getPt( i+1,j+1,k+1 ); p8.normalize(); 	curv = p8.dot(p1); 	if (curv < 0) p8 = -p8; 
		
		Point val = p1*( (1-dx)*(1-dy)*(1-dz) ) +
								p2*( dx*(1-dy)*(1-dz) ) +
								p3*( (1-dx)*dy*(1-dz) ) +
								p4*( (1-dx)*(1-dy)*dz ) +
								p5*( dx*(1-dy)*dz ) +
								p6*( (1-dx)*dy*dz ) +
								p7*( dx*dy*(1-dz) ) +
								p8*( dx*dy*dz ) ;

		val.normalize();
		return val;
	}

    /**
     *  @brief		Exports the raw field into a file.
     *  @param		fileName The filename as a STL string.
     *  @return		true if successful, false otherwise.
    **/
    bool exportRawFieldToFile(std::string fileName);

    /**
     *  @brief		Imports a raw field from a file.
     *  @param		fileName The filename as a STL string.
     *  @return		true if successful, false otherwise.
    **/
    bool importRawFieldFromFile(std::string fileName);

    // member variable definitions
  protected:
    /// Filename of Field data (with full path)
    std::string m_RawFileName;
    /// x,y,z dimensions of field
    int m_Dims[3];
    /// Pointer to the field data
    T* m_Data;
    /// Voxel size
    double m_VoxelSize[3];

};

template <typename T, int n>
Field<T, n>::Field() {
  m_RawFileName = std::string("none");
  m_Dims[0] = 0; m_Dims[1] = 0; m_Dims[2] = 0;
  m_VoxelSize[0] = 1.0; m_VoxelSize[1] = 1.0; m_VoxelSize[2] = 1.0;
  m_Data = NULL;
};

template <typename T, int n>
Field<T, n>::Field(const Field<T, n>& rhs) {
  m_RawFileName = rhs.m_RawFileName;
  m_Dims[0] = rhs.m_Dims[0]; m_Dims[1] = rhs.m_Dims[1]; m_Dims[2] = rhs.m_Dims[0];
  m_VoxelSize[0] = rhs.m_VoxelSize[0];m_VoxelSize[1] = rhs.m_VoxelSize[1];m_VoxelSize[2] = rhs.m_VoxelSize[2];
  m_Data = new T[n*m_Dims[0]*m_Dims[1]*m_Dims[2]];
  memcpy(m_Data, rhs.getDataPtr(), n*m_Dims[0]*m_Dims[1]*m_Dims[2]*sizeof(T));
}

template <typename T, int n>
Field<T,n>::~Field() {
  // delete the memory ...
  if (m_Data)
    delete [] m_Data;
  m_Data = NULL;
}

template <typename T, int n>
bool Field<T, n>::init(int sizeX, int sizeY, int sizeZ) {
  m_RawFileName = std::string("none");
  m_Dims[0] = sizeX; m_Dims[1] = sizeY; m_Dims[2] = sizeZ;
  m_VoxelSize[0] = 1.0; m_VoxelSize[1] = 1.0; m_VoxelSize[2] = 1.0;
  m_Data = new T[n*sizeX*sizeY*sizeZ];
  return true;
}

template <typename T, int n>
bool Field<T, n>::importRawFieldFromFile(std::string fileName) {
  FILE * pFile;
  unsigned long lSize;

  pFile = fopen ( fileName.c_str() , "rb" );
  if (pFile==NULL) 
    return false;

  // obtain file size.
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  if ( lSize != n*m_Dims[0]*m_Dims[1]*m_Dims[2]*sizeof(T) ) {
    fclose(pFile);
    return false;
  }

  // copy the file into the buffer.
  fread (m_Data, 1, lSize, pFile);

  // terminate
  fclose (pFile);

  return true;
}

template <typename T, int n>
bool Field<T, n>::exportRawFieldToFile(std::string fileName) {
  FILE * pFile;
  long lSize;

  lSize = n*m_Dims[0]*m_Dims[1]*m_Dims[2]*sizeof(T);

  pFile = fopen ( fileName.c_str() , "wb" );
  if (pFile==NULL) 
    return false;

  // write the data to file ...
  fwrite (m_Data, 1, lSize, pFile);

  // terminate
  fclose (pFile);

  return true;
}

template <typename T, int n>
bool Field<T, n>::init(std::string filename) {
  char tmp[256];
  char tmp2[256]; /* for scanning the .mhd file and getting the volume
                     dimensions */
  char temp[256]; /* a dummy string */

  // std::cout << "Opening file: " << str << std::endl;

  FILE *fp_dat;  /* pointer to the FILE structure for the corresponding .mhd
  */
  if ((fp_dat = fopen(filename.c_str(), "r")) == NULL) 
  { 
    // std::cout << "Could not open file" << std::endl;
    return false;
  } 

  char fName[256];
  char c[8];
  while (!feof(fp_dat)) 
  {
    fgets(temp, 255, fp_dat); 
    // std::cout << ">>: " << std::string(temp) << std::endl;
    if (ferror(fp_dat))
    {
      //wxLogMessage(wxT("Error reading Volume header file"));
		// QMessageBox::information(0, "Field Init", "Error reading Volume header file");
      return false;
    } 
    if (strncmp(temp, "ElementDataFile =", strlen("ElementDataFile =")) ==
        0) 
    {
      sscanf(temp, "%s %s %s", tmp2, c, fName);
    }
    if (strncmp(temp, "DimSize =", strlen("DimSize =")) == 0) 
    {
      sscanf(temp, "%s %s %d %d %d", tmp2, c, &m_Dims[0], &m_Dims[1],
          &m_Dims[2]);
    }
    if (strncmp(temp, "ElementSpacing =", strlen("ElementSpacing =")) == 0) 
    {
      sscanf(temp, "%s %s %lf %lf %lf", tmp2, c, &m_VoxelSize[0],
          &m_VoxelSize[1], &m_VoxelSize[2]);
    }
	
    if (strncmp(temp, "ElementType =", strlen("ElementType =")) == 0) 
    {
      sscanf(temp, "%s %s %s", tmp2, c, tmp);
      if ( !strncmp(tmp, "MET_FLOAT", strlen("MET_FLOAT")) && ( sizeof(T) != sizeof(float) ) ) {
	   // QMessageBox::information(0, "Field Init", "Data Type is not float");
        return false; // crudely done to support only float or double ... can crash at runtime with wrong type ...
      } else if ( !strncmp(tmp, "MET_DOUBLE", strlen("MET_DOUBLE")) && ( sizeof(T) != sizeof(double) ) ) {
	   // QMessageBox::information(0, "Field Init", "Data Type is not double");
        return false; // crudely done to support only float or double ... can crash at runtime with wrong type ...
      }
	}
	 /* else {
		  QMessageBox::information(0, "Field Init", "Unknown Data Type"+QString(tmp));
        return false;
      }
    }
    if (strncmp(temp, "NDims =", strlen("NDims =")) == 0) 
    {
      int ndim;
      sscanf(temp, "%s %s %d", tmp2, c, &ndim);
      if (ndim != n)
        return false;
    }*/

	/*
         if (strncmp(temp, "Position =", strlen("Position =")) == 0) 
         {
         sscanf(temp, "%s %s %lf %lf %lf", tmp2, c, &m_Corner000[0],
         &m_Corner000[1], &m_Corner000[2]);
         }
         if (strncmp(temp, "Orientation =", strlen("Orientation =")) == 0) 
         {
         sscanf(temp, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf", tmp2,
         c, &m_OrientationX[0], &m_OrientationX[1],
         &m_OrientationX[2],
         &m_OrientationY[0], &m_OrientationY[1],
         &m_OrientationY[2], &m_OrientationZ[0], &m_OrientationZ[1],
         &m_OrientationZ[2]);
         } */
  }

  fclose(fp_dat);

  // QMessageBox::information(0, "Header Read", "Successfully read MHD Image header");

  m_Data = 0;
  // int numOfEntries = m_Dims[0] * m_Dims[1] * m_Dims[2];

  // allocate memory ... 
  m_Data = new T[n*m_Dims[0]*m_Dims[1]*m_Dims[2]];

  // strip the header name on the end and replace with image name.
  std::string imgname = filename;
  std::string::size_type pos = imgname.rfind ("/", imgname.size());
	if (pos != std::string::npos)
		imgname = imgname.substr(0, pos);
	else 
		imgname = ".";
  imgname += std::string("/");
  imgname += std::string(fName);

	// std::cout << "Reading raw image " << imgname << std::endl;
  // QMessageBox::information(0, "Reading image", "Image file is " + QString(imgname.c_str()) );

  return importRawFieldFromFile( imgname );

}

#endif /*_FIELD_H_*/
