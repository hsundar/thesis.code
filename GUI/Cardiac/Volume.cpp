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
#include "Volume.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

// analyze support ...
#include "dbh.h"

Volume::Volume() :
m_VolBitsStored(0),
m_VolBitsAllocated(0),
m_VolData(0),
m_ucpLowRes(0),
m_ucpHighRes(0)
{
	m_Dims[0] = 0;
	m_Dims[1] = 0;
	m_Dims[2] = 0;

	m_VoxelSize[0] = 1.0;
	m_VoxelSize[1] = 1.0;
	m_VoxelSize[2] = 1.0;

	m_Corner000[0] = 0.0;
	m_Corner000[1] = 0.0;
	m_Corner000[2] = 0.0;

	m_OrientationX[0] = 1.0;
	m_OrientationY[0] = 0.0;
	m_OrientationZ[0] = 0.0;

	m_OrientationX[1] = 0.0;
	m_OrientationY[1] = 1.0;
	m_OrientationZ[1] = 0.0;

	m_OrientationX[2] = 0.0;
	m_OrientationY[2] = 0.0;
	m_OrientationZ[2] = 1.0;

	m_iNumberOfLandmarks = 0;

	m_Signed = false;
	m_uiMaxDepth = maxDepth;
}

Volume::~Volume()
{
	if ( m_VolData != 0 ) {
		switch ( m_VolBitsAllocated ) {
		case 8: {
				unsigned char* vol = (unsigned char*) m_VolData;
				delete [] vol;
				break;}
		case 16: {
				unsigned short* vol = (unsigned short*) m_VolData;
				delete [] vol;
				break;}
		}
	}
	// clean up the low and high res images.
	if (m_ucpLowRes != 0) {
		delete [] m_ucpLowRes;
	}
	if (m_ucpHighRes != 0) {
		delete [] m_ucpHighRes;
	}
}

bool Volume::Init(int sizeX, int sizeY, int sizeZ, int allocatedVolBits, int
									usedVolBits, bool isSigned)
{
	m_Dims[0] = sizeX;
	m_Dims[1] = sizeY;
	m_Dims[2] = sizeZ;

	m_VolBitsAllocated      = allocatedVolBits;
	m_VolBitsStored  = usedVolBits;
	m_Signed         = isSigned;

	if ( m_VolData != 0 ) {
		switch ( m_VolBitsAllocated ) {
		case 8: {
				unsigned char* vol = (unsigned char*) m_VolData;
				delete [] vol;
				break;}
		case 16: {
				unsigned short* vol = (unsigned short*) m_VolData;
				delete [] vol;
				break;}
		}
	}

	m_VolData = 0;

	int numOfEntries = m_Dims[0] * m_Dims[1] * m_Dims[2];

	switch ( m_VolBitsAllocated ) {
	case 8: {
			unsigned char* vol = new unsigned char[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;}
	case 16: {
			unsigned short* vol = new unsigned short[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;}
	}

	return true;
}

Point Volume::MM2Pixel(Point inPoint, bool roundingFlag)
{
	Point outPoint;
	Point tmpPoint;

	tmpPoint.x() =  inPoint.x()-m_Corner000[0];
	tmpPoint.y() =  inPoint.y()-m_Corner000[1];
	tmpPoint.z() =  inPoint.z()-m_Corner000[2];

	outPoint.x() =

	tmpPoint.dot(Point(m_OrientationX[0],m_OrientationX[1],m_OrientationX[2]))/   
	m_VoxelSize[0];
	outPoint.y() =

	tmpPoint.dot(Point(m_OrientationY[0],m_OrientationY[1],m_OrientationY[2]))/   
	m_VoxelSize[1];
	outPoint.z() =

	tmpPoint.dot(Point(m_OrientationZ[0],m_OrientationZ[1],m_OrientationZ[2]))/   
	m_VoxelSize[2];

	if (roundingFlag) {
		outPoint.x() = (int) (outPoint.x() + 0.5);
		outPoint.y() = (int) (outPoint.y() + 0.5);
		outPoint.z() = (int) (outPoint.z() + 0.5);
	}

	return outPoint;
}

Point Volume::Pixel2MM(Point inPoint)
{
	Point outPoint;

	outPoint.x()= m_Corner000[0] + inPoint.x() * m_VoxelSize[0] *
								m_OrientationX[0]+ inPoint.y() * m_VoxelSize[1] * m_OrientationY[0]+
								inPoint.z()     * m_VoxelSize[2] * m_OrientationZ[0];
	outPoint.y()= m_Corner000[1] + inPoint.x() * m_VoxelSize[0] *
								m_OrientationX[1]+ inPoint.y() * m_VoxelSize[1] * m_OrientationY[1]+
								inPoint.z()     * m_VoxelSize[2] * m_OrientationZ[1];
	outPoint.z()= m_Corner000[2] + inPoint.x() * m_VoxelSize[0] *
								m_OrientationX[2]+ inPoint.y() * m_VoxelSize[1] * m_OrientationY[2]+
								inPoint.z()     * m_VoxelSize[2] * m_OrientationZ[2];

	return outPoint;
}

void Volume::SetAt(int x, int y, int z, unsigned char intensity)
{
	if ( x <0 || x>= m_Dims[0] || y <0 || y>= m_Dims[1] || z<0 || z>= m_Dims[2])
		return;

	int index = z*(m_Dims[0]*m_Dims[1]) + y* m_Dims[0] + x;

	unsigned char* vol;

	vol = (unsigned char*) m_VolData;

	vol[index] = intensity;

}
void Volume::SetAt(int x, int y, int z, unsigned short intensity)
{
	if ( x <0 || x>= m_Dims[0] || y <0 || y>= m_Dims[1] || z<0 || z>= m_Dims[2])
		return;

	int index = z*(m_Dims[0]*m_Dims[1]) + y* m_Dims[0] + x;

	unsigned short* vol;

	vol = (unsigned short*) m_VolData;

	vol[index] = intensity;

}

double Volume::GetAt(Point pt) {
	int x = pt.xint();
	int y = pt.yint();
	int z = pt.zint();

	if ( x <0 || x>= m_Dims[0] || y <0 || y>=m_Dims[1] || z<0 || z>= m_Dims[2])
		return 0.;

	int index = z*(m_Dims[0]*m_Dims[1]) + y* m_Dims[0] + x;

	switch ( m_VolBitsAllocated ) {
	case 8: 
		return(double)( ((unsigned char *) m_VolData)[index]);
	case 16: 
		return(double)( ((unsigned short *) m_VolData)[index]);
	}

	return 0.;
}

void* Volume::GetAt(int x, int y, int z)
{
	if ( x <0 || x>= m_Dims[0] || y <0 || y>=m_Dims[1] || z<0 || z>= m_Dims[2])
		return NULL;

	int index = z*(m_Dims[0]*m_Dims[1]) + y* m_Dims[0] + x;

	switch ( m_VolBitsAllocated ) {
	case 8: 
		return(void *)((unsigned char*)m_VolData+index);
	case 16: 
		return(void *)((unsigned short*)m_VolData+index);
	}

	return NULL;
}

Point Volume::MM2NormalizedPixelCoords(const Point& inPoint)
{
	Point outPoint;

	outPoint  = MM2Pixel(inPoint, false);

	outPoint.x() = outPoint.x() / (double) m_Dims[0];
	outPoint.y() = outPoint.y() / (double) m_Dims[1];
	outPoint.z() = outPoint.z() / (double) m_Dims[2];

	return outPoint;
}


bool Volume::ExportRawVolumeToFile(std::string fileName)
{

	std::fstream volStream;

	volStream.open( fileName.c_str() , std::ios::out | std::ios::binary );

	if ( !volStream ) {
		// wxLogMessage("Could NOT open the volume file: %s \n",
		// fileName.c_str());
		return false;
	}

	volStream.write( (char*)m_VolData, (m_VolBitsAllocated / 8) * m_Dims[0] *
									 m_Dims[1] * m_Dims[2] );

	volStream.close();

	return true;
}

Point Volume::MM2GL3DTexCoord(const Point &inPoint)
{
	return MM2NormalizedPixelCoords(inPoint);
}

bool Volume::Init(std::string fileName)
{
	std::ifstream volInf( fileName.c_str() );

	if (!volInf) {
		m_Dims[0] = 0;
		m_Dims[1] = 0;
		m_Dims[2] = 0;

		m_VoxelSize[0] = 1.0;
		m_VoxelSize[1] = 1.0;
		m_VoxelSize[2] = 1.0;

		m_Corner000[0] = 0.0;
		m_Corner000[1] = 0.0;
		m_Corner000[2] = 0.0;

		m_OrientationX[0] = 1.0;
		m_OrientationY[0] = 0.0;
		m_OrientationZ[0] = 0.0;

		m_OrientationX[1] = 0.0;
		m_OrientationY[1] = 1.0;
		m_OrientationZ[1] = 0.0;

		m_OrientationX[2] = 0.0;
		m_OrientationY[2] = 0.0;
		m_OrientationZ[2] = 1.0;

		// wxLogMessage("Could NOT open the volume file: %s \n",
		// fileName.c_str());

		return false;
	}

	char volFileName[2048];
	int ignoreBytes;

	volInf >> volFileName;
	volInf >> ignoreBytes;
	volInf >> m_Dims[0] >> m_Dims[1] >> m_Dims[2];
	volInf >> m_VolBitsAllocated >> m_VolBitsStored;
	volInf >> m_VoxelSize[0] >> m_VoxelSize[1] >> m_VoxelSize[2];
	volInf >> m_Corner000[0] >> m_Corner000[1] >> m_Corner000[2];
	volInf >> m_OrientationX[0] >> m_OrientationX[1] >> m_OrientationX[2];
	volInf >> m_OrientationY[0] >> m_OrientationY[1] >> m_OrientationY[2];
	volInf >> m_OrientationZ[0] >> m_OrientationZ[1] >> m_OrientationZ[2];
	volInf.close();


	m_VolData = 0;
	int numOfEntries = m_Dims[0] * m_Dims[1] * m_Dims[2];

	switch ( m_VolBitsAllocated ) {
	case 8: {
			unsigned char* vol = new unsigned char[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;}
	case 16: {
			unsigned short* vol = new unsigned short[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;}
	}

	int pos1 = fileName.rfind( '\\');
	int pos2 = fileName.rfind( '/' );
	int pos = pos1 > pos2 ? pos1 : pos2;

	m_RawFileName = fileName.substr(0,  pos ) + std::string( volFileName );
	std::cout << " Raw file Name is : " << m_RawFileName << std::endl;

	std::ifstream volStream;

	volStream.open( m_RawFileName.c_str() , std::ios::binary );

	if (volStream == NULL)
		return false;

	volStream.ignore( ignoreBytes );
	volStream.read( (char*)m_VolData, (m_VolBitsAllocated / 8) * m_Dims[0] *
									m_Dims[1] * m_Dims[2] );

	if ( m_VolBitsAllocated > m_VolBitsStored ) {
		switch ( m_VolBitsAllocated ) {
		case 8: {
				unsigned char* vol = (unsigned char*) m_VolData;
				for ( int i = 0; i < m_Dims[0] * m_Dims[1] * m_Dims[2]; i++ ) {
					vol[i] <<= (m_VolBitsAllocated - m_VolBitsStored);
				}
				break;}
		case 16: {
				unsigned short* vol = (unsigned short*) m_VolData;
				for ( int i = 0; i < m_Dims[0] * m_Dims[1] * m_Dims[2]; i++ ) {
					vol[i] <<= (m_VolBitsAllocated - m_VolBitsStored);
				}
				break;}
		}
	}
	volStream.close();

	return true;
}

bool Volume::InitMhd(QString filename)
{
	char tmp[256];
	char tmp2[256];	/* for scanning the .mhd file and getting the volume
										 dimensions */
	char temp[256];	/* a dummy string */

	// std::cout << "Opening file: " << str << std::endl;

	FILE *fp_dat;	 /* pointer to the FILE structure for the corresponding .mhd
	*/
	if ((fp_dat = fopen(filename.toAscii(), "r")) == NULL) {
		// std::cout << "Could not open file" << std::endl;
		return false;
	}

	char fName[256];
	char c[8];
	while (!feof(fp_dat)) {
		fgets(temp, 255, fp_dat); 
		// std::cout << ">>: " << std::string(temp) << std::endl;
		if (ferror(fp_dat)) {
			//wxLogMessage(wxT("Error reading Volume header file"));
			return false;
		}
		if (strncmp(temp, "ElementDataFile =", strlen("ElementDataFile =")) ==
				0) {
			sscanf(temp, "%s %s %s", tmp2, c, fName);
		}
		if (strncmp(temp, "DimSize =", strlen("DimSize =")) == 0) {
			sscanf(temp, "%s %s %d %d %d", tmp2, c, &m_Dims[0], &m_Dims[1],
						 &m_Dims[2]);
		}
		if (strncmp(temp, "ElementSpacing =", strlen("ElementSpacing =")) == 0) {
			sscanf(temp, "%s %s %lf %lf %lf", tmp2, c, &m_VoxelSize[0],
						 &m_VoxelSize[1], &m_VoxelSize[2]);
		}
		if (strncmp(temp, "ElementType =", strlen("ElementType =")) == 0) {
			sscanf(temp, "%s %s %s", tmp2, c, tmp);
			if (!strncmp(tmp, "MET_UCHAR", strlen("MET_UCHAR"))) {
				m_VolBitsAllocated = 8; 
				m_VolBitsStored = 8;
			} else if (!strncmp(tmp, "MET_USHORT", strlen("MET_USHORT"))) {
				m_VolBitsAllocated = 16; 
				m_VolBitsStored = 12;
			} else if (!strncmp(tmp, "MET_FLOAT", strlen("MET_FLOAT"))) {
				m_VolBitsAllocated = 32; 
				m_VolBitsStored = 32;
			} else if (!strncmp(tmp, "MET_DOUBLE", strlen("MET_DOUBLE"))) {
				m_VolBitsAllocated = 64; 
				m_VolBitsStored = 64;
			}
		}
		if (strncmp(temp, "Position =", strlen("Position =")) == 0) {
			sscanf(temp, "%s %s %lf %lf %lf", tmp2, c, &m_Corner000[0],
						 &m_Corner000[1], &m_Corner000[2]);
		}
		if (strncmp(temp, "Orientation =", strlen("Orientation =")) == 0) {
			sscanf(temp, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf", tmp2,
						 c, &m_OrientationX[0], &m_OrientationX[1],
						 &m_OrientationX[2],
						 &m_OrientationY[0], &m_OrientationY[1],
						 &m_OrientationY[2], &m_OrientationZ[0], &m_OrientationZ[1],
						 &m_OrientationZ[2]);
		}
	}

	fclose(fp_dat);


	m_VolData = 0;
	int numOfEntries = m_Dims[0] * m_Dims[1] * m_Dims[2];

	switch ( m_VolBitsAllocated ) {
	case 8: {
			unsigned char* vol = new unsigned char[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;}
	case 16: {
			unsigned short* vol = new unsigned short[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;}
	case 32: {
			float* vol = new float[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;
		}
	case 64: {
			double* vol = new double[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;
		}
	}

	// strip the header name on the end and replace with image name.
	QString imgname = filename;
	int ii = imgname.lastIndexOf('/');
	imgname.truncate(ii);
	if (!imgname.isEmpty())
		imgname += QString("/");
	imgname += QString(fName);

	//std::cout << "Image file is " << qPrintable(imgname) << std::endl;
	std::ifstream volStream;
	volStream.open( imgname.toAscii() , std::ios::binary );

	if (volStream == NULL) {
		qDebug("Cannot open image file");
		return false;
	}

	volStream.read( (char*)m_VolData, (m_VolBitsAllocated / 8) * m_Dims[0] *
									m_Dims[1] * m_Dims[2] );

	if ( m_VolBitsAllocated > m_VolBitsStored ) {
		switch ( m_VolBitsAllocated ) {
		case 8: {
				unsigned char* vol = (unsigned char*) m_VolData;
				for ( int i = 0; i < m_Dims[0] * m_Dims[1] * m_Dims[2]; i++ ) {
					vol[i] <<= (m_VolBitsAllocated - m_VolBitsStored);
				}
				break;}
		case 16: {
				unsigned short* vol = (unsigned short*) m_VolData;
				for ( int i = 0; i < m_Dims[0] * m_Dims[1] * m_Dims[2]; i++ ) {
					vol[i] <<= (m_VolBitsAllocated - m_VolBitsStored);
				}
				break;}
		}
	}

	// Rescale data to unsigned char for ushort, float and double

	switch (m_VolBitsAllocated) {
	case 16: {
			unsigned short* vol = (unsigned short*) m_VolData;
			unsigned char* newvol = new unsigned char[numOfEntries] ;
			unsigned short _max =  0; 
			unsigned short _min =  5000; 
			for ( int i = 0; i < numOfEntries; i++ ) {
				if (vol[i] > _max) _max = vol[i];
				if (vol[i] < _min) _min = vol[i];
			}
			// Now rescale ...
			for (int i=0; i<numOfEntries; i++) {
				newvol[i] = static_cast<unsigned char>( (254*(vol[i] - _min))/(_max - _min) );
			}

			delete [] vol;
			m_VolData = newvol;
			m_VolBitsAllocated = m_VolBitsStored = 8;
			break;
		}
		case 32: {
			float* vol = (float*) m_VolData;
			unsigned char* newvol = new unsigned char[numOfEntries] ;
			float _max =  -5000; 
			float _min =  5000; 
			for ( int i = 0; i < numOfEntries; i++ ) {
				if (vol[i] > _max) _max = vol[i];
				if (vol[i] < _min) _min = vol[i];
			}
			// Now rescale ...
			for (int i=0; i<numOfEntries; i++) {
				newvol[i] = static_cast<unsigned char>( (254*(vol[i] - _min))/(_max - _min) );
			}

			delete [] vol;
			m_VolData = newvol;
			m_VolBitsAllocated = m_VolBitsStored = 8;
			break;
		}
		case 64: {
			double* vol = (double*) m_VolData;
			unsigned char* newvol = new unsigned char[numOfEntries] ;
			double _max =  -5000; 
			double _min =  5000; 
			for ( int i = 0; i < numOfEntries; i++ ) {
				if (vol[i] > _max) _max = vol[i];
				if (vol[i] < _min) _min = vol[i];
			}
			// Now rescale ...
			for (int i=0; i<numOfEntries; i++) {
				newvol[i] = static_cast<unsigned char>( (254*(vol[i] - _min))/(_max - _min) );
			}

			delete [] vol;
			m_VolData = newvol;
			m_VolBitsAllocated = m_VolBitsStored = 8;
			break;
		}
	}

	volStream.close();


	//std::cout << "Creating low and high res images" << std::endl;
	// create low and high Res images ...

	// determine sizes ... for now it is /4 and /2 ... no smart determination

	// allocate memory ...

	/*
	int szLow = (m_Dims[0]/4) * (m_Dims[1]/4) * (m_Dims[2]/4);
	int szHigh = (m_Dims[0]/2) * (m_Dims[1]/2) * (m_Dims[2]/2);

	m_ucpLowRes  =  new unsigned char [szLow];
	m_ucpHighRes =  new unsigned char [szHigh];


	//std::cout << "Dims are " << m_Dims[0] << " " << m_Dims[1] << " " << m_Dims[2] << std::endl;

	//std::cout << "Allocated memory for low and high res images" << std::endl;
	// Now set values ...

	switch ( m_VolBitsAllocated ) {
	case 8:{
			unsigned char* vol = (unsigned char*) m_VolData;
			// low-res first
			for (int k=0; k<m_Dims[2]; k+=4) {
				for (int j=0; j<m_Dims[1]; j+=4) {
					for (int i=0; i<m_Dims[0]; i+=4) {
						int idx = ( (k/4)*(m_Dims[1]/4) + (j/4))*(m_Dims[0]/4) + (i/4);
						if (idx < szLow) {
							m_ucpLowRes[idx] = ((vol[(k*m_Dims[1]+j)*m_Dims[0]+i]) >> 2);	// histograms bin size is 6
							if (m_ucpLowRes[idx] > 63) {
								std::cout << "LowRes is out of range " << idx << " = " << m_ucpLowRes[idx] << std::endl;
							}
						}
					}
				}
			}
			// now high-res
			for (int k=0; k<m_Dims[2]; k+=2) {
				for (int j=0; j<m_Dims[1]; j+=2) {
					for (int i=0; i<m_Dims[0]; i+=2) {
						int idx = ((k/2)*(m_Dims[1]/2) + (j/2))*(m_Dims[0]/2) + (i/2);
						if (idx < szHigh) {
							m_ucpHighRes[idx] = ((vol[(k*m_Dims[1]+j)*m_Dims[0]+i]) >> 2); // histograms bin size is 6
						}
					}
				}
			}

			break; 
		}
	case 16: {
			unsigned short* vol = (unsigned short*) m_VolData;
			// low-res first
			for (int k=0; k<m_Dims[2]; k+=4) {
				for (int j=0; j<m_Dims[1]; j+=4) {
					for (int i=0; i<m_Dims[0]; i+=4) {
						int idx = ( (k/4)*(m_Dims[1]/4) + (j/4))*(m_Dims[0]/4) + (i/4);
						if (idx < szLow) {
							m_ucpLowRes[idx] = ((vol[(k*m_Dims[1]+j)*m_Dims[0]+i]) >> 10); // histograms bin size is 6
						}
					}
				}
			}
			// now high-res
			for (int k=0; k<m_Dims[2]; k+=2) {
				for (int j=0; j<m_Dims[1]; j+=2) {
					for (int i=0; i<m_Dims[0]; i+=2) {
						int idx = ((k/2)*(m_Dims[1]/2) + (j/2))*(m_Dims[0]/2) + (i/2);
						if (idx < szHigh) {
							m_ucpHighRes[idx] = ((vol[(k*m_Dims[1]+j)*m_Dims[0]+i]) >> 10);	// histograms bin size is 6
						}
					}
				}
			}

			break;
		}
	}
	*/
	// std::cout << "Finished creating low and high res images" << std::endl;

	return true;
}

void Volume::ComputeHistogram(unsigned int *histogram, int nrBins)
{
	int totalSize = m_Dims[0]*m_Dims[1]*m_Dims[2];


	for (int i=0; i<nrBins; i++) {
		histogram[i] = 0;
	}

	for (int i=0; i< totalSize; i++) {
		int voxelIntensity = 0;

		if (m_VolBitsAllocated == 8)
			voxelIntensity=(int) *((unsigned char  *)m_VolData + i);
		else if (m_VolBitsAllocated == 16)
			voxelIntensity=(int) *((unsigned short *)m_VolData + i);

		if ( voxelIntensity < nrBins )
			histogram[voxelIntensity]++;
	}
}

bool Volume::InitAnalyze(QString filename)
{

	// @bug .. handle swap for the binary case too ...
	// first read in the image header
	struct dsr hdr; 
	FILE *fp;
	bool doSwap = false;

	if ((fp=fopen(filename.toAscii(),"r"))==NULL) {
		return false;
	}

	fread(&hdr,1, sizeof(struct dsr), fp);

	if (hdr.dime.dim[0] < 0 || hdr.dime.dim[0] > 15) {
		swap_hdr(&hdr);
		doSwap=true;
	}

	// now populate the Volume structure with the info ...

	// dimensions ...
	m_Dims[0] = hdr.dime.dim[1];
	m_Dims[1] = hdr.dime.dim[2];
	m_Dims[2] = hdr.dime.dim[3];

	// qDebug("Dims are %d,%d,%d", m_Dims[0], m_Dims[1], m_Dims[2]);
	// spacing   
	m_VoxelSize[0] = hdr.dime.pixdim[1];
	m_VoxelSize[1] = hdr.dime.pixdim[2];
	m_VoxelSize[2] = hdr.dime.pixdim[3];

	bool _signed = false;
	// data type ...
	if ( hdr.dime.datatype == DT_UNSIGNED_CHAR ) {
		m_VolBitsAllocated = 8; 
		m_VolBitsStored = 8;
	} else if ( hdr.dime.datatype == DT_SIGNED_SHORT) {
		m_VolBitsAllocated = 16; 
		m_VolBitsStored = 12;
		_signed = true;
	} else {
		return false;
	}

	// qDebug("Datatype: %d, %d", m_VolBitsAllocated, m_VolBitsStored);
	// Offset ... where is this ?
	m_Corner000[0] = 0.0;
	m_Corner000[1] = 0.0;
	m_Corner000[2] = 0.0;


	// orientation ..
	m_OrientationX[0] = 0.0; m_OrientationX[1] = 0.0; m_OrientationX[2] = 0.0;
	m_OrientationY[0] = 0.0; m_OrientationY[1] = 0.0; m_OrientationY[2] = 0.0;
	m_OrientationZ[0] = 0.0; m_OrientationZ[1] = 0.0; m_OrientationZ[2] = 0.0;

	fclose(fp);

	m_VolData = 0;
	int  numOfEntries  =  m_Dims[ 0]  *  m_Dims[ 1]  *  m_Dims[ 2];

	switch ( m_VolBitsAllocated ) {
	case 8: {
			unsigned char* vol = new unsigned char[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;}
	case 16: {
			unsigned short* vol = new unsigned short[ numOfEntries ];
			for ( int i = 0; i < numOfEntries; i++ ) vol[i] = 0;
			m_VolData = vol;
			break;}
	}

	// strip the header name on the end and replace with image name.
	QString imgname = filename;
	int ii = imgname.lastIndexOf('.');
	imgname.truncate(ii);
	imgname += QString(".img");

	// qDebug(imgname.toLatin1());

	std::ifstream volStream;
	volStream.open( imgname.toAscii() , std::ios::binary );

	if (volStream == NULL) {
		qDebug("failed to open image file");
		return false;
	}

	volStream.read( (char*)m_VolData, (m_VolBitsAllocated / 8) * m_Dims[0] *
									m_Dims[1] * m_Dims[2] );

	if ( m_VolBitsAllocated > m_VolBitsStored ) {
		switch ( m_VolBitsAllocated ) {
		case 8: {
				unsigned char* vol = (unsigned char*) m_VolData;
				for ( int i = 0; i < m_Dims[0] * m_Dims[1] * m_Dims[2]; i++ ) {
					vol[i] <<= (m_VolBitsAllocated - m_VolBitsStored);
				}
				break;}
		case 16: {
				unsigned short* vol = (unsigned short*) m_VolData;
				for ( int i = 0; i < m_Dims[0] * m_Dims[1] * m_Dims[2]; i++ ) {
					vol[i] <<= (m_VolBitsAllocated - m_VolBitsStored);
				}
				break;}
		}
	}
	volStream.close();

	return true;
}

void Volume::ComputeLinearOctree() {
	// Need to have working list of nodes ...
	std::vector<Node>   workingList;
	unsigned int octSz = 1 << m_uiMaxDepth;

	// insert root into the working list ...
	workingList.push_back( Node() );

	for ( unsigned int i=0; i < workingList.size(); i++ ) {
		Node nd = workingList[i];
		Point _off( nd.x, nd.y, nd.z );
		_off /= (double)octSz;
		int sz = 1 << nd.d; 
		// qDebug("Processing node %d, %d, added %d", i, nd.d, m_Octree.size());
		Point _sz(1./sz, 1./sz, 1./sz);

		if ( isUnique(_off, _sz) ) {
			// is unique ... add it to the linear Octree ...
			m_Octree.push_back(nd);
		} else {
			// is not unique .. split it ...
			if (! nd.addChildren(workingList) )
				nd.addChildren(workingList);
			// m_Octree.push_back(nd);
		}
		// workingList.erase(workingList.begin());
	}
	workingList.clear();
	std::cout << "Added " << m_Octree.size() << " octants" << std::endl;

	std::sort (m_Octree.begin(), m_Octree.end());
	// for now lets write it to a file ...
	qDebug("Writing out oct file");
	std::ofstream out("test.oct");
	out << "3 " << maxDepth << std::endl;
	out << m_Octree.size() << std::endl;
	for (unsigned int i=0; i<m_Octree.size(); i++)
		out << m_Octree[i].x << " " << m_Octree[i].y << " " << m_Octree[i].z << " " << m_Octree[i].d << std::endl;
	out.close();

}

// might want to change this to take a node as input ...
bool Volume::isUnique(Point _off, Point _sz) {
	// determine the lower and upper limts ...
	Point lower, upper, mid;

	lower.x() = _off.x() * m_Dims[0];  
	lower.y() = _off.y() * m_Dims[1];  
	lower.z() = _off.z() * m_Dims[2];

	upper.x() = lower.x() + _sz.x()*m_Dims[0];
	upper.y() = lower.y() + _sz.y()*m_Dims[1];
	upper.z() = lower.z() + _sz.z()*m_Dims[2];

	mid = (lower + upper)/2;

	unsigned char val = *((unsigned char *)GetAt(mid.xint(), mid.yint(), mid.zint()));
	// unsigned char val = *((unsigned char *)GetAt((int)lower.x(), (int)lower.y(), (int)lower.z()));
	// qDebug ("mid val is %d", val);
	// qDebug ("low is %f, %f, %f\n", floor(lower.x()), floor(lower.y()), floor(lower.z()));
	// qDebug ("high is %f, %f, %f\n", ceil(upper.x()), ceil(upper.y()), ceil(upper.z()));

	// std::cout << "val is " << (int)val << std::endl;
	for (int k=(int)(floor(lower.z())); k<(int)(ceil(upper.z())); k++)
		for (int j=(int)(floor(lower.y())); j<(int)(ceil(upper.y())); j++)
			for (int i=(int)(floor(lower.x())); i<(int)(ceil(upper.x())); i++) {
				unsigned char val2 = *((unsigned char *)GetAt(i, j, k));
				// std::cout << val << ", val2 is " << val2 << std::endl;
				//if ( abs(val2 - val) > 16 )
				if (val2 != val)
					return false;
			}

	return true;
}

