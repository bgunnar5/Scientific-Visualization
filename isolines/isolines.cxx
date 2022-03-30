/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <math.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
//  Function: interpolate
//
//  Purpose: 
//     A helper function that does interpolation for us using the equations
//     provided in class:
//                      F(X) = F(A) + t*(F(B) - F(A))
//                      t = (X - A) / (B - A)
//
//  Arguments:
//       X (input):     a float representing the value we're interpolating
//                      for
//       A (input):     a float representing the A value in our equation
//       B (input):     a float representing the B value in our equation
//       FA (input):    a float representing the F(A) value in our equation
//       FB (input):    a float representing the F(B) value in our equation
//
//   Returns: the interpolated value.
//
// ****************************************************************************
float interpolate(float X, float A, float B, float FA, float FB) {
    float t = (X - A) / (B - A);
    float interpolatedValue = FA + t * (FB - FA);
    return interpolatedValue;
}


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

// ****************************************************************************
//  Function: GenerateNumSegments
//
//  Purpose: 
//     For each case 0-15, provide the number of edges in each case.
//
//  Arguments:
//       numSegments (output): an array we're storing edge cases in
//
// ****************************************************************************
void GenerateNumSegments(int *numSegments) {
  // 0b0000; No lines; Opposite coloring of index 15; All field values at vertices are greater than isovalue
  numSegments[0] = 0;

  // 0b0001; Line from edge 0 to edge 1; Opposite coloring of index 14
  numSegments[1] = 1;

  // 0b0010; Line from edge 1 to edge 2; Opposite coloring of index 13
  numSegments[2] = 1;

  // 0b0011; Line from edge 0 to edge 2; Opposite coloring of index 12
  numSegments[3] = 1;

  // 0b0100; Line from edge 0 to edge 3; Opposite coloring of index 1
  numSegments[4] = 1;

  // 0b0101; Line from edge 1 to edge 3; Opposite coloring of index 10
  numSegments[5] = 1;

  // 0b0110; Line from edge 0 to edge 3 and line from edge 1 to edge 2
  numSegments[6] = 2;

  // 0b0111; Line from edge 2 to edge 3; Opposite coloring of index 8
  numSegments[7] = 1;

  // 0b1000; Line from edge 2 to edge 3; Opposite coloring of index 7
  numSegments[8] = 1;

  // 0b1001; Line from edge 0 to edge 1 and line from edge 2 to edge 3
  numSegments[9] = 2;

  // 0b1010; Line from edge 1 to edge 3; Opposite coloring of index 5
  numSegments[10] = 1;

  // 0b1011; Line from edge 0 to edge 3; Opposite coloring of index 4
  numSegments[11] = 1;

  // 0b1100; Line from edge 0 to edge 2; Opposite coloring of index 3
  numSegments[12] = 1;

  // 0b1101; Line from edge 1 to edge 2; Opposite coloring of index 2
  numSegments[13] = 1;
  
  // 0b1110; Line from edge 0 to edge 1; Opposite coloring of index 1
  numSegments[14] = 1;

  // 0b1111; No lines; Opposite coloring of index 0; All field values at vertices are greater than isovalue
  numSegments[15] = 0;
}

// ****************************************************************************
//  Function: GenerateLookupTable
//
//  Purpose: 
//     Populate the lookup table with where line segments are for each case.
//
//  Arguments:
//       lup (output):                a 2D array storing information about what edges each line segment in the case touches (if any)
//       numSegments (input):         an int array storing the number of edges per case
//       sizeNumSegments (input):     an int representing the size of the numSegments array
//
// ****************************************************************************
void GenerateLookupTable(int lup[16][4], int *numSegments, int sizeNumSegments) {
    // Handling the -1 entries first
    for (int i = 0; i < sizeNumSegments; i++) {
        // No line segments so every entry in lup at this index should be -1
        if (numSegments[i] == 0) {
            for (int j = 0; j < 4; j++) {
                lup[i][j] = -1;
            }
        }
        // One line segment so last 2 entries in lup at this index should be -1
        else if (numSegments[i] == 1) {
            for (int j = 2; j < 4; j++) {
                lup[i][j] = -1;
            }
        }
    }

    // Indices 1 and 14 have a line segment from edge 0 to edge 1
    lup[1][0] = lup[14][0] = 0;
    lup[1][1] = lup[14][1] = 1;

    // Indices 2 and 13 have a line segment from edge 1 to edge 2
    lup[2][0] = lup[13][0] = 1;
    lup[2][1] = lup[13][1] = 2;

    // Indices 3 and 12 have a line segment from edge 0 to edge 2
    lup[3][0] = lup[12][0] = 0;
    lup[3][1] = lup[12][1] = 2;

    // Indices 4 and 11 have a line segment from edge 0 to edge 3
    lup[4][0] = lup[11][0] = 0;
    lup[4][1] = lup[11][1] = 3;

    // Indices 5 and 10 have a line segment from edge 1 to edge 3
    lup[5][0] = lup[10][0] = 1;
    lup[5][1] = lup[10][1] = 3;

    // Index 6 has two line segments from edge 0 to edge 3 and edge 1 to edge 2
    lup[6][0] = 0;
    lup[6][1] = 3;
    lup[6][2] = 1;
    lup[6][3] = 2;

    // Indices 7 and 8 have a line segment from edge 2 to edge 3
    lup[7][0] = lup[8][0] = 2;
    lup[7][1] = lup[8][1] = 3;

    // Index 9 has two line segments from edge 0 to edge 1 and edge 2 to edge 3
    lup[9][0] = 0;
    lup[9][1] = 1;
    lup[9][2] = 2;
    lup[9][3] = 3;
}

// ****************************************************************************
//  Function: IdentifyCase
//
//  Purpose: 
//     Figure out which case (0 - 15) we're dealing with for a specific cell.
//
//  Arguments:
//       isovalue (input):  the isovalue we're comparing field values to (provided for us in the project handout)
//       F (input):         a float representing the A value in our equation
//       v0 (input):        the point index for vertex 0
//       v1 (input):        the point index for vertex 1
//       v2 (input):        the point index for vertex 2
//       v3 (input):        the point index for vertex 3
//
//   Returns: the case number.
//
// ****************************************************************************
int IdentifyCase(float isovalue, float *F, int v0, int v1, int v2, int v3) {
    // Initialize case number to 0
    int icase = 0;

    // If the field value at vertex 0 is greater than the isovalue, add 1 to the case (2^0 spot in the binary number)
    if (F[v0] > isovalue) {
        icase++;
    }

    // If the field value at vertex 1 is greater than the isovalue, add 2 to the case (2^1 spot in the binary number)
    if (F[v1] > isovalue) {
        icase += 2;
    }

    // If the field value at vertex 2 is greater than the isovalue, add 4 to the case (2^2 spot in the binary number)
    if (F[v2] > isovalue) {
        icase += 4;
    }

    // If the field value at vertex 3 is greater than the isovalue, add 8 to the case (2^3 spot in the binary number)
    if (F[v3] > isovalue) {
        icase += 8;
    }

    return icase;
}

// ****************************************************************************
//  Function: interpolatePoint
//
//  Purpose: 
//     Interpolate to find the x and y coordinates of the end of a line segment based on the edge.
//
//  Arguments:
//       pt (output):       a float array to store the x and y coordinates that we interpolated for
//       isovalue (input):  the isovalue we're comparing field values to (provided for us in the project handout)
//       F (input):         a float representing the A value in our equation
//       v0 (input):        the point index for vertex 0
//       v1 (input):        the point index for vertex 1
//       v2 (input):        the point index for vertex 2
//       v3 (input):        the point index for vertex 3
//       X:                 an array containing the X locations of a rectilinear mesh.
//       Y:                 an array containing the Y locations of a rectilinear mesh.
//       vertices (input):  a 2D array storing the logical point indices of each vertex.
//       edge (input):      an int representing the edge we need to interpolate on
//
// ****************************************************************************
void interpolatePoint(float *pt, float isovalue, float *F, int v0, int v1, int v2, int v3, float *X, float *Y, int vertices[4][2], int edge) {
    // Based on the edge that has the point, interpolate the point's x and y coordinates
    switch (edge) {
        // Edge between v0 and v2
        case 0:
            pt[0] = interpolate(isovalue, F[v0], F[v2], X[vertices[0][0]], X[vertices[2][0]]);
            pt[1] = Y[vertices[0][1]];
            break;
        // Edge between v0 and v1
        case 1:
            pt[0] = X[vertices[0][0]];
            pt[1] = interpolate(isovalue, F[v0], F[v1], Y[vertices[0][1]], Y[vertices[1][1]]);
            break;
        // Edge between v1 and v3
        case 2:
            pt[0] = interpolate(isovalue, F[v1], F[v3], X[vertices[1][0]], X[vertices[3][0]]);
            pt[1] = Y[vertices[1][1]];
            break;
        // Edge between v2 and v3
        case 3:
            pt[0] = X[vertices[2][0]];
            pt[1] = interpolate(isovalue, F[v2], F[v3], Y[vertices[2][1]], Y[vertices[3][1]]);
            break;
        // Should never be here
        default:
            break;
    }   
}

int main()
{
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("isolines.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    SegmentList sl;
    sl.AddSegment(-10, -10, +10, -10); // Add segment (-10,-10) -> (+10, -10)
    sl.AddSegment(-10, +10, +10, +10);
    sl.AddSegment(-10, -10, -10, +10);
    sl.AddSegment(+10, -10, +10, +10);

    // Create a variable to represent the isovalue so it's easy to change if we want to
    float isovalue = 3.2;

    // Generate an array to store the number of line segments going through the cell for each of the 16 cases
    int sizeNumSegments = 16;
    int numSegments[sizeNumSegments];
    GenerateNumSegments(numSegments);

    // Generate the lookup table to determine which edges in a cell have lines connecting them
    int lup[sizeNumSegments][4];
    GenerateLookupTable(lup, numSegments, sizeNumSegments);

    // Get the number of cells for looping
    int numCells = GetNumberOfCells(dims);

    for (int i = 0; i < numCells; i++) {
        // Find the logical index of the cell (this corresponds to the point index at the lower left vertex of the cell)
        int logicalCellIdx[2];
        GetLogicalCellIndex(logicalCellIdx, i, dims);
        
        // Establish the 4 vertices in the order v0, v1, v2, v3
        // v0 = bottom left; v1 = top left; v2 = bottom right; v3 = top right
        int vertices[4][2] = {{logicalCellIdx[0], logicalCellIdx[1]}, {logicalCellIdx[0], logicalCellIdx[1] + 1}, {logicalCellIdx[0] + 1, logicalCellIdx[1]}, {logicalCellIdx[0] + 1, logicalCellIdx[1] + 1}};

        // Get the point index for every vertex in the cell
        int v0PtIdx = GetPointIndex(vertices[0], dims);
        int v1PtIdx = GetPointIndex(vertices[1], dims);
        int v2PtIdx = GetPointIndex(vertices[2], dims);
        int v3PtIdx = GetPointIndex(vertices[3], dims);

        // Get the case we're looking for
        int icase = IdentifyCase(isovalue, F, v0PtIdx, v1PtIdx, v2PtIdx, v3PtIdx);

        // Get the number of line segments for this case
        int nsegments = numSegments[icase];

        // Loop through each line segment
        for (int j = 0; j < nsegments; j++) {
            // Find what edges the line segments are on
            int edge1 = lup[icase][2*j];
            int edge2 = lup[icase][2*j+1];

            // Initialize variables to store the points
            float pt1[2];
            float pt2[2];

            // Interpolate to find the correct points to add a line between
            interpolatePoint(pt1, isovalue, F, v0PtIdx, v1PtIdx, v2PtIdx, v3PtIdx, X, Y, vertices, edge1);
            interpolatePoint(pt2, isovalue, F, v0PtIdx, v1PtIdx, v2PtIdx, v3PtIdx, X, Y, vertices, edge2);

            sl.AddSegment(pt1[0], pt1[1], pt2[0], pt2[1]);
        }
    }

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!

    vtkPolyData *pd = sl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
