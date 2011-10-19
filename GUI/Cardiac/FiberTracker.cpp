#include "FiberTracker.h"
#include <algorithm>
#include <iostream>

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoVertexProperty.h>
#include <Inventor/nodes/SoLineSet.h>

#include <QtGui>

FiberTracker::FiberTracker() {
	m_pFA = NULL;
	m_pPD = NULL;
	m_iColor = Gray;
}

FiberTracker::~FiberTracker() {

}

SoSeparator* FiberTracker::track() {
	// QMessageBox::information(0, "Fiber Track", "in Track");
	// Inventor stuff ...
	SoSeparator *root = new SoSeparator;
	root->ref();
	
	// compute the transformation ...
	
	SoTransform *fibTrans = new SoTransform;
	float tx1 = -m_pFA->getSize(0)*m_pFA->getVoxelSize(0)*0.5;
	float ty1 = -m_pFA->getSize(1)*m_pFA->getVoxelSize(1)*0.5;
	float tz1 = -m_pFA->getSize(2)*m_pFA->getVoxelSize(2)*0.5;		
    // isoTrans->scaleFactor.setValue(sx1, sy1, sz1);
    fibTrans->translation.setValue(tx1, ty1, tz1);
    root->addChild(fibTrans);
        
	SoVertexProperty *verts = new SoVertexProperty;
	verts->ref();
	verts->materialBinding = SoVertexProperty::PER_VERTEX;
	// root->addChild(verts);
	SoLineSet *lines = new SoLineSet;
	lines->ref();
	lines->vertexProperty = verts;
	root->addChild(lines);

	// global cnt of pts. ..
	long cnt =0, pcnt=0;
	// preliminaries ...
	int x,y,z;
        double sx, sy, sz;
	x = m_pFA->getSize(0);
	y = m_pFA->getSize(1);
	z = m_pFA->getSize(2);
	
        sx = m_pFA->getVoxelSize(0);
	sy = m_pFA->getVoxelSize(1);
	sz = m_pFA->getVoxelSize(2);

	int smp = (int)(1/m_fSeedDensity);

	// QMessageBox::information(0, "Fiber Track", "Generating seeds");
	// First generate seeds ...
        for (int k=0; k<z; k+=smp)
          for (int j=0; j<y; j+=smp)
            for (int i=0; i<x; i+=smp) {
              if ( ( *(m_pFA->getAt(i,j,k)) > FA_THRESHOLD )  ) 
                seeds.push_back(Point(i, j, k));
            }
	//	QMessageBox::information(0, "Fiber Track", "finished Generating seeds");

        // std::cout << seeds.size() << std::endl;
	// Now grow Fibers ...
	// Grow fibers at each step. Stop if FA is violated. If length is longer than
	// a certain threshold, then add it as an inverntor indexed line set ...
	QProgressDialog progress("Tracking Fibers...", "Abort", 0, seeds.size(), 0);
	progress.setValue(0);

	for (unsigned int p=0; p<seeds.size(); p++) {
		fiber.clear();
		norms.clear();
		int i;
		progress.setValue(p);
		qApp->processEvents();

                if (progress.wasCanceled())
                  break;

		// track for each seed
		Point lastV(0.,0.,0.);
		Point curr = seeds[p];
		Point last = seeds[p];
		Point next = seeds[p];
		// grow in negative direction ...
		for (i=0; i<max_iter; i++) {
			Point vel = getVelocity(curr, lastV);
			if (!i)
				vel = -vel;

			next = curr + vel*step_sz;
			lastV = vel;
			fiber.push_back(curr);
			norms.push_back(vel);

			if (lastV.abs() < 0.01) {
				// std::cout << "Velocity Failure " <<  fiber.size() << std::endl;
				break;
			}
			// Test to see if the FA is all right. Later on will add additional
			// conditions at this stage like testing for curvature etc. Lets keep it
			// simple for the time being.
			if ( *(m_pFA->getAt((int)curr.x(), (int)curr.y(), (int)curr.z())) < FA_THRESHOLD) {
				// std::cout << "FA failure" << " " << fiber.size() <<  std::endl;
				break;
			}

			// move to the new location if all is well :)
			last = curr;
			curr = next;
		}
		// reverse the list, since push_front doesn't work
		std::reverse(fiber.begin(), fiber.end());
		std::reverse(norms.begin(), norms.end());
		//norms.reserve(0);

		// track for each seed
		lastV = Point(0.,0.,0.);
		curr = seeds[p];
		last = seeds[p];
		next = seeds[p];
		// grow in positive direction ...
		for (i=0; i<max_iter; i++) {
			Point vel = getVelocity(curr, lastV);
			// std::cout << tc(DIM, Y, BL) << "Velocity: " << vel << tc(RST, W, BL) << std::endl;
			next = curr + vel*step_sz;
			lastV = vel;
			fiber.push_back(curr);
			norms.push_back(vel);

			if (lastV.abs() < 0.01) {
				// std::cout << "Velocity Failure " << fiber.size() << std::endl;
				break;
			}
			// Test to see if the FA is all right. Later on will add additional
			// conditions at this stage like testing for curvature etc. Lets keep it
			// simple for the time being.
			if ( *(m_pFA->getAt((int)curr.x(), (int)curr.y(), (int)curr.z())) < FA_THRESHOLD) {
				// std::cout << "FA failure "  << fiber.size() << std::endl;
				break;
			}

			// move to the new location if all is well :)
			last = curr;
			curr = next;
		}

		// now render it ....
		// add it to the root ...
		if (fiber.size() > 3000) {
			for (unsigned int ii=1; ii<fiber.size()-1; ii++){
				// loop through all pts in fiber ... and add the points,
				if ( (fiber[ii].x() > x)  || (fiber[ii].y() >y) || (fiber[ii].z() > z) || (fiber[ii].x() <0)  || (fiber[ii].y() <0) || (fiber[ii].z() < 0) )
					QMessageBox::information(0, "Fiber Track", QString("Wrong fiber point %1").arg(fiber[ii].abs()));
				verts->vertex.set1Value(cnt, sx*fiber[ii].x(), sy*fiber[ii].y(), sz*fiber[ii].z());
				// Point tang, norm;
				// tang = norms[ii];
				// norm = tang.cross(Point(0.,0.,1.0));
				// verts->normal.set1Value(cnt, norm.x(), norm.y(), norm.z());
				// verts->normal.set1Value(cnt, norms[ii].x(), norms[ii].y(), norms[ii].z());
				if (m_iColor == Normal) {
					unsigned int col = pt2Col(norms[ii], fiber[ii]);
					verts->orderedRGBA.set1Value(cnt, col);
				} else if (m_iColor == Chiral) {
                                  Point tang = norms[ii];
                                  tang.normalize();
                                  Point rad = Point((double)x/2, (double)y/2, fiber[ii].z());
                                  rad.normalize();
                                  float pitch = tang.z() * (rad.cross(tang)).z() ;
                                  
                                  // std::cout << pitch << std::endl;
                                  // color based on the pitch ...
                                  unsigned int col = float2Col(pitch);
                                  verts->orderedRGBA.set1Value(cnt, col);
                                }
				cnt++;
			}
			// now add one line set corresponding to this ...
			lines->numVertices.set1Value(pcnt, fiber.size()-2);
			pcnt++;
		}
	}
	progress.setValue(seeds.size());
	QMessageBox::information(0, "Fiber Track", QString("Added total of %1 fibers with %2 points").arg(pcnt).arg(cnt));
	return root;
}

Point FiberTracker::getVelocity ( int i, int j, int k) {
	float *vel = m_pPD->getAt(i,j,k);
	return Point(vel[0], vel[1], vel[2]);
}

// Susumo's version of getVelocity at a point.
Point FiberTracker::getVelocity (Point p, Point lastV) {
	// Here we will always use the nearest neighbour to get the principal
	// direction. The vel will be scaled to hit the next boundary, and also we
	// need to test if the eigenvector needs to be inverted. This probably means
	// that we will need to send in the previous velocity estimate.

	int i, j, k;
	i = (int)p.x(); j = (int)p.y(); k = (int)p.z();

	// get PD from the floored value of p
	Point pd = getVelocity( i,j,k );

	// Handle the case for the first step. Need to scale the PD before returning
	// it.
	if (lastV.abs() < 0.001) {
		pd.normalize();
		return pd;
	}

	// Compare it with the last estimate of velocity. Enforce Curvature here. We
	// can change the test in the main loop to one that checks the velocity.
	// Invert the eigenvector if necessary.
	float curv = (float)(pd.dot(lastV)/(pd.abs()*lastV.abs()));
	// first check if curvature is ok. If the abs is less than 0.8 then we should
	// set the velocity to zero and return. CHANGE 0.8 to a define or a static
	// variable that we read in from some config file.
	if ( fabs(curv) < 0.8 )
		return Point(0, 0, 0);

	// Since the vectors are reasonable aligned, lets flip the PD if necessary.
	if (curv < 0 ) {
		// The vectors are misaligned. Flip the PD.
		pd = - pd; // Add negation operator to the Point class.
	}

	// Scale the eigenvector to reach the next boundary. Simple first run can be
	// to simply normalize it to a unit vector which along with a step size of 1
	// should produce similar results.
	pd.normalize();

	return pd;
}

Point FiberTracker::getVelocityInterpolated (Point p ) {
	// Function takes a float index and returns velocity at that point after
	// interpolation. Currently, linear interpolation is used. Might need to
	// change that eventually.
	int i, j, k;
	float x, y, z;

	i = (int)p.x(); j = (int)p.y(); k = (int)p.z();
	x = p.x() - i;  y = p.y() - j;  z = p.z() - k;

	// Simple tri-linear Interpolation ....
	Point val = getVelocity( i,j,k )*( (1-x)*(1-y)*(1-z) ) +
		getVelocity( i+1,j,k )*( x*(1-y)*(1-z) ) +
		getVelocity( i,j+1,k )*( (1-x)*y*(1-z) ) +
		getVelocity( i,j,k+1 )*( (1-x)*(1-y)*z ) +
		getVelocity( i+1,j,k+1 )*( x*(1-y)*z ) +
		getVelocity( i,j+1,k+1 )*( (1-x)*y*z ) +
		getVelocity( i+1,j+1,k )*( x*y*(1-z) ) +
		getVelocity( i+1,j+1,k+1 )*( x*y*z ) ;

	return val;
}

unsigned int FiberTracker::pt2Col(Point p, Point q) {
	Point cvec = p;
	cvec.normalize();
	cvec *= 255;
	int x = (int)q.x(); int y = (int)q.y(); int z = (int)q.z(); 

	unsigned int r, g, b, a, clr;
	r = (int)fabs(cvec.x()); g = (int)fabs(cvec.y()); b = (int)fabs(cvec.z());
	a = (unsigned int)(255*(1.0 - *(m_pFA->getAt(x,y,z))));
	
	clr = (r << 24) + (g << 16) + (b << 8) + a; // opacity is zero for now ...
	return clr;
}

unsigned int FiberTracker::float2Col(float p) {
  // basically an interpolation between two colors ...

  float q = (p + 1.0)/2;
  unsigned int cw = 0xff000000; // Purple ...
  unsigned int ccw = 0x0000ff00; // Cyan ...

  return (unsigned int)(q*cw + (1.-q)*ccw);

}
