/**********************************************
  Hari Sundar - 8 Nov 2005
  Myocyte Tracking Program

  --* Need to add anyoption     *--
  --* Need to handle exceptions *--
  *********************************************/

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>

#include "Point.h"

// VTK Includes ...
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyLine.h"
#include "vtkDataSetMapper.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkCubeSource.h"
#include "vtkFloatArray.h"
#include "vtkCubeAxesActor2D.h"
#include "vtkTubeFilter.h"
#include "vtkLookupTable.h"
#include "vtkScalarBarActor.h"

class triplet {
public:
  int i;
  int j;
  int k;

  triplet(int x, int y, int z) {i=x; j=y; k=z;}
};

// Other Functions...
Point getVelocity(triplet ind, float *vel, triplet sz);
Point getVelocity(Point p, float *vel, triplet sz, Point lastV);
Point getVelocityInterpolated(Point p, float *vel, triplet sz);

int main(int argc, char **argv) {

  if (argc < 5) {
    std::cout << "Usage: " <<  argv[0] ; 
    // The args    1  2  3 4 5   6
    std::cout << " PD FA x y z sampling" << std::endl;  
    return -1;
  }

  std::cout << "Fiber Tracking Program" <<  std::endl;
 
  // Simple variables ...
  int x, y, z, smp;
  x = atoi(argv[3]); y = atoi(argv[4]); z = atoi(argv[5]); smp = atoi(argv[6]);
  int max_iter = 500000;
  float step_sz = 1.0;
  

  std::cout <<  "Allocating Memory:" << x << " " << y << " " << z; //<< tc(RST, W) << std::endl;
  // Allocate memory ...
  float * fa = new float [x*y*z];
  float * pd = new float [3*x*y*z];
  std::cout <<  " Success"  << std::endl;
  
  // Containers ...
  std::vector<Point> seeds;
  std::vector<Point> fiber;
  std::vector<Point> norms;
    
  std::cout <<  "Reading in anisotropy file "  << argv[2] << std::endl;
  // Read in the FA map
  std::ifstream fa_in (argv[2], std::ios_base::binary);
  fa_in.read((char *)fa, x*y*z*sizeof(float) );
  fa_in.close();

  std::cout << "Reading in principal direction "  <<  argv[1] <<  std::endl;
  // Read in the PD data
  std::ifstream pd_in (argv[1], std::ios_base::binary);
  pd_in.read((char *)pd, 3*x*y*z*sizeof(float));
  pd_in.close();
	
#ifdef _DEBUG
  float _mx=-100.0, _mn=100.0;
  for (int i=0; i<x*y*z; i++){
	if (fa[i] > _mx) _mx = fa[i];
	if (fa[i] < _mn) _mn = fa[i];
  }
  std::cout << "FA Range is " << _mn << " " << _mx << std::endl;
#endif

  // VTK Init. for the polylines ...
  vtkPolyData *data = vtkPolyData::New();
  vtkFloatArray *vecs = vtkFloatArray::New();
  vecs->SetNumberOfComponents(3);
  vtkFloatArray *scalars = vtkFloatArray::New();
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *lines = vtkCellArray::New();
  
  // Compute the set of seed points based on FA ...
  for (int k=45; k<46; k+=smp)
    for (int j=0; j<y; j+=smp)
      for (int i=0; i<x; i+=smp) {
        if ( (fa[k*x*y + j*x +i] > 0.1 )  ) // && (fa[k*x*y + j*x + i] < 0.5))
          seeds.push_back(Point(i, j, k));
      }

  std::cout <<  "Computed " << seeds.size() << " Points based on anisotropy"  << std::endl; 
 

  std::ofstream out("fibers.txt");
  out << seeds.size() << std::endl;
  
  // Grow fibers at each step. Stop if FA is violated. If length is longer than
  // a certain threshold, then add it as a VTK polyline ...
  for (int p=0; p<seeds.size(); p++) {
    //std::cout << p << " of " << seeds.size() << std::endl;
    fiber.clear();
	norms.clear();
	int i;

	// track for each seed
	Point lastV(0,0,0);
    Point curr = seeds[p];
    Point last = seeds[p];
    Point next = seeds[p];
	// grow in negative direction ...
	for (i=0; i<max_iter; i++) {
		Point vel = getVelocity(curr, pd, triplet(x,y,z), lastV);
		if (!i)
			vel = -vel;
		// std::cout << tc(DIM, Y, BL) << "Velocity: " << vel << tc(RST, W, BL) << std::endl;
		next = curr + vel*step_sz;
		lastV = vel;
		fiber.push_back(curr);
		norms.push_back(vel);
		
		if (lastV.abs() < 0.01) {
			// std::cout << tc(RST, R, BL) << "Velocity Failure" << lastV << " " << fiber.size() << tc(RST, W, BL) << std::endl;
			break;
		}
		// Test to see if the FA is all right. Later on will add additional
		// conditions at this stage like testing for curvature etc. Lets keep it
		// simple for the time being.
		if (fa[curr.getIndex(x,y)] < 0.1) {
			//std::cout << tc(RST, R, BL) << "FA failure" << fa[curr.getIndex(x,y)] << " " << fiber.size() << tc(RST, W, BL) <<  std::endl;
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
    lastV = Point(0,0,0);
    curr = seeds[p];
    last = seeds[p];
    next = seeds[p];
	// grow in positive direction ...
    for (i=0; i<max_iter; i++) {
      Point vel = getVelocity(curr, pd, triplet(x,y,z), lastV);
      // std::cout << tc(DIM, Y, BL) << "Velocity: " << vel << tc(RST, W, BL) << std::endl;
      next = curr + vel*step_sz;
      lastV = vel;
      fiber.push_back(curr);
	  norms.push_back(vel);

	if (lastV.abs() < 0.01) {
	  // std::cout << tc(RST, R, BL) << "Velocity Failure" << lastV << " " << fiber.size() << tc(RST, W, BL) << std::endl;
          break;
	}
      // Test to see if the FA is all right. Later on will add additional
      // conditions at this stage like testing for curvature etc. Lets keep it
      // simple for the time being.
      if (fa[curr.getIndex(x,y)] < 0.1) {
        //std::cout << tc(RST, R, BL) << "FA failure" << fa[curr.getIndex(x,y)] << " " << fiber.size() << tc(RST, W, BL) <<  std::endl;
        break;
      }
        
      // move to the new location if all is well :)
      last = curr;
      curr = next;
    }
    

    
	if (fiber.size() > 30) {
      // std::cout << tc(RST, G, BL) << "Point " << p << ": Computed a fiber of length " << fiber.size() << " points " <<tc(RST, W, BL) <<std::endl;
      // if of sufficient length ... decide to render it ...

      out << fiber.size() << std::endl; 
          
      int num_p =  points->GetNumberOfPoints();
      vtkIdList *ids = vtkIdList::New();
      
	  // Color by x,y,z
/*
	  for (unsigned int ii=0; ii<fiber.size(); ii++){
                  out << fiber[ii];
		  double v[3];
		  v[0] = norms[ii].x(); v[1] = norms[ii].y(); v[2] = norms[ii].z(); 
		  
		  double sc;
		  if ( (fabs(v[0]) > fabs(v[1])) && (fabs(v[0]) > fabs(v[2])) )
			  sc = 0.66+ fabs(v[0])/3;
		  else if ( fabs(v[1]) > fabs(v[2]) ) 
			  sc = fabs(v[1])/3;
		  else
			  sc = 0.33 + fabs(v[2])/3;

		  // std::cout << sc << std::endl;
		  vtkIdType _id = points->InsertNextPoint(fiber[ii].x(), fiber[ii].y(), 3*fiber[ii].z());
		  ids->InsertId(ii,_id);
		  vecs->InsertTuple(_id, v);
		  scalars->InsertTuple(_id, &sc);		
      }*/

	  // color by rising / falling ...
	  double sc =  fiber.size(); //[fiber.size()].z() - fiber[0].z();
          sc = seeds[p].x();
	  // std::cout << sc << std::endl;
	  for (unsigned int ii=0; ii<fiber.size()-1; ii++){
                  out << fiber[ii];
		  float v[3];
		  v[0] = norms[ii].x(); v[1] = norms[ii].y(); v[2] = norms[ii].z(); 
		  
		  double sd = sc;
		  //std::cout << sc << std::endl;
		  vtkIdType _id = points->InsertNextPoint(fiber[ii].x(), fiber[ii].y(), 3*fiber[ii].z());
		  ids->InsertId(ii,_id);
		  vecs->InsertTuple(_id, v);
		  scalars->InsertTuple(_id, &sd);
		  
      }

      lines->InsertNextCell(ids);
      ids->Delete();
    }
  }

  out.close();
  // Finished with all seed points ...
  
  std::cout <<  "Rendering " << lines->GetNumberOfCells() << " fibers"  << std::endl; 
  
  // Clean Up ...
  std::cout << "Cleaning Up ... " <<  std::endl;
  delete [] pd;
  delete [] fa;
  
  // Render usage text ...
  std::cout << std::endl << "Now Rendering. Please use the following keys to navigate." << std::endl << std::endl;
  std::cout << "L-Mouse: \t\t"  << "Rotate" << std::endl;
  std::cout << "R-Mouse/Scroll Wheel: \t"  << "Zoom" << std::endl;
  std::cout << "M-Mouse: \t\t"  << "Move" << std::endl;
  std::cout << std::endl  << "q/Q \t\t\t"  << "Quit" <<  std::endl;
  
  
  // Finish and Render in VTK.

  data->SetPoints(points);
  int idd = data->GetPointData()->AddArray(scalars);
  data->GetPointData()->SetActiveAttribute(idd, vtkDataSetAttributes::SCALARS);
  data->SetLines(lines);
  data->GetPointData()->SetVectors(vecs);
  data->Update();

  // Set up the lookup tables ...

  double *rg = data->GetPointData()->GetScalars()->GetRange();

  // Start by creatin a black/white lookup table.
  vtkLookupTable *bwLut = vtkLookupTable::New();
  bwLut->SetTableRange (rg);
  bwLut->SetSaturationRange (0, 0);
  bwLut->SetHueRange (0, 0);
  bwLut->SetValueRange (0, 1);
  
  // Now create a lookup table that consists of the full hue circle
  // (from HSV).
  vtkLookupTable *hueLut = vtkLookupTable::New();
  hueLut->SetTableRange (rg);
  hueLut->SetHueRange (0, 1);
  hueLut->SetSaturationRange (1, 1);
  hueLut->SetValueRange (1, 1);
  
  // Finally, create a lookup table with a single hue but having a range
  // in the saturation of the hue.
  vtkLookupTable *satLut = vtkLookupTable::New();
  satLut->SetTableRange (rg);
  satLut->SetHueRange (.6, .6);
  satLut->SetSaturationRange (0, 1);
  satLut->SetValueRange (1, 1);
 
  vtkDataSetMapper *aPolyLineMapper = vtkDataSetMapper::New();
  aPolyLineMapper->SetInput(data);
  aPolyLineMapper->SetScalarRange(data->GetPointData()->GetScalars()->GetRange());
  aPolyLineMapper->SetLookupTable(hueLut);  

  // Set Up a Scalar bar ...
  vtkScalarBarActor *colorbar = vtkScalarBarActor::New();
  colorbar->SetLookupTable(aPolyLineMapper->GetLookupTable());
  colorbar->SetTitle("Fibers");
  colorbar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  colorbar->GetPositionCoordinate()->SetValue(0.88, 0.1);
  colorbar->SetOrientationToVertical();
  colorbar->SetWidth(0.1);
  colorbar->SetHeight(0.8);

  vtkActor *aPolyLineActor = vtkActor::New();
  aPolyLineActor->SetMapper(aPolyLineMapper);
  aPolyLineActor->AddPosition(-x/2, -y/2, -1.5*z);

  // The Bounding box ... just to get an idea ...
  vtkCubeSource *cubeSrc = vtkCubeSource::New();
  cubeSrc->SetBounds(-x/2, x/2, -y/2, y/2, -1.5*z, 1.5*z);
  vtkPolyDataMapper *cubeMap = vtkPolyDataMapper::New();
  cubeMap->SetInput(cubeSrc->GetOutput());
  
  vtkActor *cube = vtkActor::New();
  cube->SetMapper(cubeMap);
  cube->GetProperty()->SetRepresentationToWireframe();
  
  vtkRenderer *ren= vtkRenderer::New();
  ren->SetBackground( 0.2, 0.4, 0.8 );
  // ren->SetBackground( 0.8, 0.8, 0.8 );
  ren->AddActor(aPolyLineActor);
  ren->AddActor(cube);
  ren->AddActor2D(colorbar);

  // Axes Labeling ...
  vtkCubeAxesActor2D *axes = vtkCubeAxesActor2D::New();
  axes->SetInput(cubeSrc->GetOutput());
  axes->SetCamera(ren->GetActiveCamera());
  axes->SetLabelFormat("%6.4g");
  axes->SetFlyModeToOuterEdges();
  axes->SetFontFactor(0.8);
 // ren->AddProp(axes); 
  
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->SetSize(600,600);
  renWin->AddRenderer( ren );  
  ren->ResetCamera();
  
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  
  vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
  iren->SetInteractorStyle(style);
  
  iren->Initialize();
  iren->Start();

  return 0;
}


Point getVelocity (triplet ind, float *vel, triplet sz) {
  // std::cout << "getVel(int) " << ind.i << " " << ind.j << " " << ind.k << std::endl; 
  int index = 3*(ind.k*sz.i*sz.j + ind.j*sz.i + ind.i);
  return Point(vel[index], vel[index+1], vel[index+2]);
}

// Susumo's version of getVelocity at a point.
Point getVelocity (Point p, float *vel, triplet sz, Point lastV) {
  // Here we will always use the nearest neighbour to get the principal
  // direction. The vel will be scaled to hit the next boundary, and also we
  // need to test if the eigenvector needs to be inverted. This probably means
  // that we will need to send in the previous velocity estimate.

  int i, j, k;
  i = (int)p.x(); j = (int)p.y(); k = (int)p.z();
  
  // get PD from the floored value of p
  Point pd = getVelocity( triplet(i,j,k), vel, sz);
  
  // Handle the case for the first step. Need to scale the PD before returning
  // it.
  if (lastV.abs() < 0.001) {
    pd.normalize();
    return pd;
  }
  
  // Compare it with the last estimate of velocity. Enforce Curvature here. We
  // can change the test in the main loop to one that checks the velocity.
  // Invert the eigenvector if necessary.
  float curv = pd.dot(lastV)/(pd.abs()*lastV.abs());
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

Point getVelocityInterpolated (Point p, float *vel, triplet sz) {
  // Function takes a float index and returns velocity at that point after
  // interpolation. Currently, linear interpolation is used. Might need to
  // change that eventually.
  int i, j, k;
  float x, y, z;

  i = (int)p.x(); j = (int)p.y(); k = (int)p.z();
  x = p.x() - i;  y = p.y() - j;  z = p.z() - k;

  // Simple tri-linear Interpolation ....
  Point val = getVelocity( triplet(i,j,k), vel, sz )*( (1-x)*(1-y)*(1-z) ) +
    getVelocity( triplet(i+1,j,k), vel, sz )*( x*(1-y)*(1-z) ) +
    getVelocity( triplet(i,j+1,k), vel, sz )*( (1-x)*y*(1-z) ) +
    getVelocity( triplet(i,j,k+1), vel, sz )*( (1-x)*(1-y)*z ) +
    getVelocity( triplet(i+1,j,k+1), vel, sz )*( x*(1-y)*z ) +
    getVelocity( triplet(i,j+1,k+1), vel, sz )*( (1-x)*y*z ) +
    getVelocity( triplet(i+1,j+1,k), vel, sz )*( x*y*(1-z) ) +
    getVelocity( triplet(i+1,j+1,k+1), vel, sz )*( x*y*z ) ;

  return val;
}
