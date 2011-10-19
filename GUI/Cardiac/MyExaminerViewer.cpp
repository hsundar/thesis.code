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
#include "MyExaminerViewer.h"
#include <Inventor/nodes/SoSwitch.h>
#include <iostream>

MyExaminerViewer::MyExaminerViewer(QWidget * parent)
	: SoQtExaminerViewer(parent, NULL, TRUE,
						 SoQtFullViewer::BUILD_ALL, SoQtFullViewer::BROWSER,
                            // build == FALSE, to delay creation of decorations
						 FALSE)
{
    // Explicitly trigger the construction of viewer decorations.
	QWidget * widget = this->buildWidget(this->getParentWidget());

	this->setBaseWidget(widget);
        this->setDecoration(false);
	
        // default values ... will be changed ...
	m_iLow = 0;
	m_iHigh = 8;
        m_iDepth = 8;
	
	// Need to set the values of the wheels too ...
	
	this->setLeftWheelString("High (8)");
	this->setBottomWheelString("Low (0)");
}
  
void
MyExaminerViewer::createViewerButtons(QWidget * parent, SbPList * buttonlist)
{
	SoQtExaminerViewer::createViewerButtons(parent, buttonlist);
}

void
MyExaminerViewer::leftWheelMotion(float value) {
  // Use this to adjust the high threshold for the display of octree
  // depths.

  int _win = (int)(-value*10);

  if (_win < 0) _win = 0;
  if (_win >= m_iDepth) _win = m_iDepth -1;

  char str[30];
  sprintf(str, "High (%d)", _win);
  this->setLeftWheelString(str);
  
  m_iHigh = _win;

  // call the update function to turn off display.
  updateOct();
}

void
MyExaminerViewer::bottomWheelMotion(float value) {	
  // Use this to adjust the low threshold for the display of octree
  // depths.

  int _win = (int)(value*10);

  if (_win < 0) _win = 0;
  if (_win >= m_iDepth) _win = m_iDepth -1;

  char str[30];
  sprintf(str, "Low (%d)", _win);
  this->setBottomWheelString(str);

  m_iLow = _win;
  updateOct();
}

void
MyExaminerViewer::updateOct() {
  for (int i=0; i<m_iDepth; i++) {
    // get the right switch ...
    /*
	  char str[32];
    sprintf(str, "switch%d", i);
    SoSwitch * sw = (SoSwitch *)SoNode::getByName(str);

    if ( (i >= m_iLow) && (i <= m_iHigh) )
      sw->whichChild = SO_SWITCH_ALL;
    else 
      sw->whichChild = SO_SWITCH_NONE;
  */
  }
}

