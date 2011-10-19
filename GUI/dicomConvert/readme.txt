dicomConvert v0.1

author:  	Hari Sundar
date:		29 Jan 2007

dicomConvert is a simple application that allows users to convert DICOM images into the Analyze and MetaIO formats. 

It provides a simple GUI to browse directories with DICOM datasets and displays all detected datasets. The user can then select datasets which need to be converted. Currently 3D and 4D (MR Cine) datasets are supported.

The 4D support is only for MR Cine Sequences, and it will in all probability not work with other forms of 4D sequences.  

Pre-Requisites
--------------
The application is written in C++ using the Qt toolkit from Trolltech. The dcmtk toolkit is used to parse dicom datasets. These two libraries need to be installed prior to installing dicomConvert. The version requirements for these libraries are,

Qt  > 4.x
	Can be obtained from http://www.trolltech.com 
	(Application will not work with Qt 3)

dcmtk	> 3.5.x 


Installation
------------

ALL OSes. Currently tested on Linux and MS Windows. It should work on all platforms where Qt 4 can be installed. Check Qt documentation for more information regarding your platform.

1. Make sure Qt and dcmtk are installed.
2. Make sure the Qt binary directory $QTDIR/bin is in the system path. 
   You will need to modify the systems PATH variable for this.
   Set environment variable DCMTK_DIR to point to the directory where dcmtk was installed. ( Recommended )
3. Check the dicomConvert.pro file and modify to suit your system. Check that,
	the variables LIBPATH and INCLUDEPATH point to the correct directories for dcmtk.
4. Invoke qmake in the root directory. Calling
   
	qmake dicomConvert.pro 

   will generate a makefile. This is usually good for *nix systems. On MS Windows, if using the MSVC compiler, it 
   might be easier to generate a VC project file. This can be generated using,

	qmake -tp vc dicomConvert.pro

5. Build the project using 'make' on *nix systems and if using the MSVC compiler (without the project) using 'nmake'.
6. The executible should be ready.

Usage
-----