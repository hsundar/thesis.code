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

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>

class MyExaminerViewer : public SoQtExaminerViewer {
	public:
		MyExaminerViewer(QWidget * parent);
  
                void setMaxDepth(int depth) { m_iDepth = depth; };
                void updateOct();
                int getDepth () { return m_iDepth; };
	protected:
		virtual void createViewerButtons(QWidget * parent, SbPList * buttonlist);
		virtual void leftWheelMotion(float value);

		virtual void bottomWheelMotion(float value);
		
		// variables ... used to threshold the depth's to be displayed.
                int		m_iDepth;
		int m_iLow;
		int m_iHigh;
};
  

