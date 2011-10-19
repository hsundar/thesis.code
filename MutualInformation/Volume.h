/***************************************************************************
*   Copyright (C) 2005 by Hari sundar   *
*   hsundar@seas.upenn.edu   *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along with this program; if not, write to the                         *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/

#ifndef VOLUME_H
#define VOLUME_H

#include "Point.h"
#include <vector>
#include <string>
#include <qstring.h>

#define maxDepth 8

class Node {
public:
	Node():x(0), y(0), z(0), d(0) {};
	Node(int i, int j, int k, int l): x(i), y(j), z(k), d(l) {};
	~Node() {};

	bool addChildren(std::vector <Node> &children) { 
		// unsigned int dim = d;

		if( d == maxDepth) {
			return false;
		}

		//The check that lev < maxD is taken care of in the constructor.

		unsigned int len = (unsigned int)(1<<( maxDepth - (d+1) ) ) ;

		Node   tmpNode0( x, y, z, d+1);
		children.push_back(tmpNode0);

		Node   tmpNode1( x+len, y, z, d+1);
		children.push_back(tmpNode1);

		Node   tmpNode2( x, y+len, z, d+1);
		children.push_back(tmpNode2);

		Node   tmpNode3( x+len, y+len, z, d+1);
		children.push_back(tmpNode3);

		Node   tmpNode4( x, y, z+len, d+1);
		children.push_back(tmpNode4);

		Node   tmpNode5( x+len, y, z+len, d+1);
		children.push_back(tmpNode5);

		Node   tmpNode6( x, y+len, z+len, d+1);
		children.push_back(tmpNode6);

		Node   tmpNode7( x+len, y+len, z+len, d+1);
		children.push_back(tmpNode7);

		return true;
	}//end function

	Point center() {
		Point p1(x,y,z);
		int len = (unsigned int)(1<<( maxDepth - d ) ) ;
		Point p2(len, len, len);
		return p1 + p2/2;
	}

	bool operator == ( Node   const & other)  const{

		if( (x == other.x) && (y == other.y) && (z == other.z) && ( d == other.d ) ) {
			return true;
		}else {
			return false;
		}
	}//end fn.

	bool operator != (Node  const  & other)  const{
		return (!((*this) == other));
	}//end fn.

	bool operator  < (Node   const & other)  const{
		// first compare the x, y, and z to determine which one dominates ...
		//Ancestor is smaller in all spaces.
		if( ( x == other.x) && (y == other.y) && (z == other.z) ) {
			return ( d < other.d );	
		}//end if	

		unsigned int _x = (x^other.x);
		unsigned int _y = (y^other.y);
		unsigned int _z = (z^other.z);

		//Default pref: y > z > x in all spaces.
		unsigned int maxC = _y;
		unsigned int zOrx = _z;
		if(zOrx < _x) { if( (_x^zOrx) >= zOrx ) {zOrx = _x;} }
		if(maxC < zOrx) { if( (maxC^zOrx) >= maxC ) {maxC = zOrx;} }

		if(maxC == _y) { return (y < other.y); }
		else if(maxC == _z) {  return (z < other.z); }
		else {  return (x < other.x); }
	}//end function

	bool operator  <= ( Node  const  & other)  const{
		return ( ((*this) < other) || ((*this) == other) );
	}//end fn.

	bool operator  > ( Node  const  & other)  const{
		return ( (!((*this) < other)) && (!((*this) == other)) );
	}//end fn.

	bool operator  >= ( Node  const & other)  const{
		return ( !((*this) < other) ) ;
	}//end fn.


	int x;
	int y;
	int z;
	int d;
};

class Volume  
{
public:	
	/**
	*  @brief		Constructor
	*
	**/
	Volume();

	/**
	*  @brief		Destructor
	*
	**/
	virtual ~Volume();

	/**
	*  @brief		initizlizes the volume 
	*  @param		the filename
	containing the volume settings
	*  @return		true if
	successfull, false otherwisr
	*
	**/
	bool Init(std::string fileName);
	bool InitMhd(QString fileName);
	bool InitAnalyze(QString fileName);

	/**
	*  @brief		initizlizes
	the volume 
	*  @param	
	three sizes and type of the volume volBits either 8 or 16
	*  @return	
	true if successfull, false otherwisr
	*
	**/
	bool
		Init(int sizeX, int sizeY, int sizeZ, int allocatedVolBits, int					
		usedVolBits, bool isSigned);

	/**
	*  @brief	
	return the size of the volume in the specified direction.

	*  @param		dim	Dimension, must be 0, 1, or 2.

	*  @return		dimension size.

	*
	**/

	int GetSize( int dim ) const { return m_Dims[dim]; };

	/**

	*  @brief		return the the center of the original volume in pixels.

	*  @return		Ramp3DPoint for center.
	**/

	Point GetCenterInPixel() {return												
		Point(m_Dims[0]/2.0,m_Dims[1]/2.0,m_Dims[2]/2.0);};

	/**

	*  @brief		return the the center of the original volume in pixels.

	*  @return		Ramp3DPoint for center.
	**/

	Point GetCenterInMM() {return Pixel2MM(GetCenterInPixel());};


	int GetVolBitsAllocated() const { return m_VolBitsAllocated; };


	int GetVolBitsStored() const { return m_VolBitsStored; };

	/**

	*  @brief		return the size of a voxel in the specified direction.

	*  @param		dim	Dimension, must be 0, 1, or 2.

	*  @return		voxel size in specified direction.

	*
	**/

	double GetVoxelSize( int dim ) const { return m_VoxelSize[dim]; };

	/**

	*  @brief		set the size of a voxel in the specified direction.

	*  @param		dim	Dimension, must be 0, 1, or 2.

	*
	**/

	void SetVoxelSize( int dim, double voxelSize ) { m_VoxelSize[dim] =				
		voxelSize; };

	/**

	*  @brief		return the coordinates of the volume corner 000

	*  @param		dim	Dimension, must be 0, 1, or 2.

	*  @return		coordinate in specified direction.

	*
	**/

	double GetCorner000( int dim ) const { return m_Corner000[dim]; };

	Point GetCorner000() const { return												
		Point(m_Corner000[0],m_Corner000[1],m_Corner000[2]);
	};

	/**

	*  @brief		set the coordinates of the volume corner 000

	*  @param		dim	Dimension, must be 0, 1, or 2.

	*
	**/

	void SetCorner000( int dim, double corner ) { m_Corner000[dim] = corner; };

	void SetCorner000( double* corner ) { m_Corner000[0] =							
		corner[0];m_Corner000[1] =
		corner[1];m_Corner000[2] = corner[2]; };

	void SetCorner000( Point inPoint ) { m_Corner000[0] =							
		inPoint.x();m_Corner000[1] =
		inPoint.y();m_Corner000[2] = inPoint.z(); };

	/**

	*  @brief		set the orientation of the volume it should be normalized!
	**/

	void SetOrientationX( double* orientation ) { m_OrientationX[0] =				

		orientation[0];m_OrientationX[1] = orientation[1];m_OrientationX[2] =			
		orientation[2];
	};

	void SetOrientationY( double* orientation ) { m_OrientationY[0] =				

		orientation[0];m_OrientationY[1] = orientation[1];m_OrientationY[2] =			
		orientation[2];
	};

	void SetOrientationZ( double* orientation ) { m_OrientationZ[0] =				

		orientation[0];m_OrientationZ[1] = orientation[1];m_OrientationZ[2] =			
		orientation[2];
	};


	void SetOrientationX( Point orientation ) { m_OrientationX[0] =					

		orientation.x();m_OrientationX[1] = orientation.y();m_OrientationX[2] =			
		orientation.z();
	};

	void SetOrientationY( Point orientation ) { m_OrientationY[0] =					

		orientation.x();m_OrientationY[1] = orientation.y();m_OrientationY[2] =			
		orientation.z();
	};

	void SetOrientationZ( Point orientation ) { m_OrientationZ[0] =					

		orientation.x();m_OrientationZ[1] = orientation.y();m_OrientationZ[2] =			
		orientation.z();
	};



	double GetOrientationX( int dim ) const { return m_OrientationX[dim]; };

	double GetOrientationY( int dim ) const { return m_OrientationY[dim]; };

	double GetOrientationZ( int dim ) const { return m_OrientationZ[dim]; };


	Point GetOrientationX() const { return											

		Point(m_OrientationX[0],m_OrientationX[1],m_OrientationX[2]); };

	Point GetOrientationY() const { return											

		Point(m_OrientationY[0],m_OrientationY[1],m_OrientationY[2]); };

	Point GetOrientationZ() const { return											

		Point(m_OrientationZ[0],m_OrientationZ[1],m_OrientationZ[2]); };

	/**

	*  @brief		returns the file name of the raw volume.

	*  @return		CString for the filename.
	**/

	const std::string& GetRawFileName() const { return m_RawFileName; };

	/**

	*  @brief		return pointer to data.

	*  @return		pointer to the volume data.

	*
	**/

	void* GetDataPtr() const { return m_VolData; };

	/**

	*  @brief		set or get specific pixel point in the volume container.

	*  @param		dimension and the value

	*  @note        for this change to take effect in displaying CVRTObject

	*               has to be re-initialized.
	**/

	void SetAt(int x, int y, int z, unsigned char intensity);

	void SetAt(int x, int y, int z, unsigned short intensity);

	void* GetAt(int x, int y, int z);
	double GetAt(Point pt);

	/**

	*  @brief		converts physical coords. of a point into pixel coordinate.

	*  @param		input point (in imaging coordinate system).

	*  @return		pixel coordinate of the input point.

	*
	**/

	Point Pixel2MM(Point inPoint);

	/**

	*  @brief		converts pixel coords of the input point to the physical		
	coordinate.

	*  @param		input point in pixel within the volume bounds.

	*  @return		output point in imaging coordinate system.

	*
	**/

	Point MM2Pixel(Point inPoint, bool roundingFlag);

	/**

	*  @brief		converts physical coords. of a point normalized into pixel		
	coordinate.

	*  @param		input point in imaging coordinate system.

	*  @return		output point in normalized pixel coordinate system (bounds		
	of the
	volume is mapped to 0 and 1).

	*
	**/

	Point MM2NormalizedPixelCoords(const Point& inPoint);

	/**

	*  @brief		converts physical coords. of a point normalized into pixel		
	coordinate.

	*  @param		input point in imaging coordinate system.

	*  @return		output point in normalized pixel coordinate system (bounds		
	of the
	volume is mapped to 0 and 1).

	*  @note        useful for the texture mapping.
	**/

	Point MM2GL3DTexCoord( const Point& inPoint);

	/**

	*  @brief		Exports the raw volume into a file.

	*  @param		CString filename.

	*  @return		true if successful, false otherwise.

	*
	**/

	bool ExportRawVolumeToFile(std::string fileName);


	int GetHistSize() const {return (1 << m_VolBitsAllocated);};


	bool IsSigned(){return m_Signed;}


	void ComputeHistogram(unsigned int * histogram, int nrBins);

	// Octree ...
	void ComputeLinearOctree();
	void SetMaxDepth(unsigned int d) { m_uiMaxDepth = d; };
	std::vector<Node>& getOctree() { return m_Octree; };

	unsigned int octSize() {return m_Octree.size(); }
	Node& getNode(unsigned int i) { return m_Octree[i]; }
	Point getNodeCenter(unsigned int i) { return m_Octree[i].center() ; }


	// Landmarks related functions ...

	void AddLandmark(Point pt, int i) {

		if (i) { 

			m_vecLandmarks2.push_back(pt); 

		} else {

			m_vecLandmarks1.push_back(pt); 

			m_iNumberOfLandmarks++; 

		}

	}


	int GetNumberOfLandmarks() { return m_iNumberOfLandmarks; }

	Point GetLandmark(unsigned int i, int j) 

	{


		if (j) {

			if ( i < m_vecLandmarks2.size() )

				return m_vecLandmarks2[i];

			else

				return Point(0,0,0);

		} else {

			if ( i < m_vecLandmarks1.size() )

				return m_vecLandmarks1[i];

			else

				return Point(0,0,0);

		}

	}	

protected:
	bool isUnique(Point _off, Point _sz);

	/// Filename of volume data (with full path)
	std::string m_RawFileName;
	/// x,y,z dimensions of volume
	int m_Dims[3];
	/// Number of storage bits per voxel
	int m_VolBitsStored;
	/// Number of bits used in original volume
	int m_VolBitsAllocated;
	/// Pointer to the volume data
	void* m_VolData;
	/// Voxel size
	double m_VoxelSize[3];
	/// Coordinates of center of voxel 000 in imaging coordinate system
	double m_Corner000[3];
	/// Coordinates of orientation vector (Normalized) 
	double m_OrientationX[3];
	double m_OrientationY[3];
	double m_OrientationZ[3];

	// For Landmarks ....
	std::vector<Point> m_vecLandmarks1;
	std::vector<Point> m_vecLandmarks2;
	int m_iNumberOfLandmarks;

	// Octree stuff ...
	std::vector<Node> 			m_Octree;
	unsigned int						m_uiMaxDepth;

	bool m_Signed;
};

#endif // VOLUME_H
