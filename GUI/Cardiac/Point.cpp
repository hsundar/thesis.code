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
#include "Point.h"

Point::Point()
{
	_x=0.0;
	_y=0.0;
	_z=0.0;
}

Point::~Point()
{

}

Point::Point(const Point &newposition) {
	_x = newposition._x;
	_y = newposition._y;
	_z = newposition._z;
}

Point Point::operator - () const {
  return Point(-_x, -_y, -_z);
}

void Point::operator *= (const int factor){
	_x*=factor;
	_y*=factor;
	_z*=factor;
}

void Point::operator *= (const double factor){
	_x*=factor;
	_y*=factor;
	_z*=factor;
}

void Point::operator /= (const int divisor){
	if(divisor == 0) return;
	_x/=(double) divisor;
	_y/=(double) divisor;
	_z/=(double) divisor;
}

void Point::operator /= (const double divisor){
	if(divisor == 0) return;
	_x/= divisor;
	_y/= divisor;
	_z/= divisor;
}

void Point::operator += (const Point& other){
	_x += other._x;
	_y += other._y;
	_z += other._z;
}

Point Point::operator - (const Point &other){
	return Point(_x-other._x,_y-other._y, _z-other._z);
}

Point Point::operator - (const Point &other) const {
	return Point(_x-other._x,_y-other._y, _z-other._z);
}

Point Point::operator + (const Point &other){
	return Point(_x+other._x,_y+other._y, _z+other._z);
}

Point& Point::operator=(const Point &other){
	_x = other._x;
	_y = other._y;
	_z = other._z;
	return *this;
}

Point Point::operator /(const double divisor)
{
	return Point(_x/divisor,_y/divisor, _z/divisor);
}

Point Point::operator *(const double factor)
{
	return Point(_x*factor,_y*factor, _z*factor);
}

Point Point::TransMatMultiply(double *transMat, Point inPoint)
{
	Point outPoint;

	outPoint._x = transMat[ 0]*inPoint._x +transMat[ 4]*inPoint._y +transMat[8]
			*inPoint._z +transMat[12];
	outPoint._y = transMat[ 1]*inPoint._x +transMat[ 5]*inPoint._y +transMat[9]
			*inPoint._z +transMat[13];
	outPoint._z = transMat[ 2]*inPoint._x +transMat[ 6]*inPoint._y
			+transMat[10]*inPoint._z +transMat[14];

	return outPoint;
}

void Point::normalize() {
  double _abs = sqrt(_x*_x + _y*_y + _z*_z);
  if ( _abs > 1e-3) {
    _x /= _abs; _y /= _abs; _z /= _abs;
  } else {
    _x=0.; _y=0.; _z=0.;
  }
}

