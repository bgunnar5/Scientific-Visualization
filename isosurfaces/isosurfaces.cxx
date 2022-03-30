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
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include "LUT.h"
#include <math.h>

// Constant used for debugging; 1 means debug is on, 0 means it's off
const int DEBUG = 0;

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
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
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
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
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
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
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
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
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
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    // idx[0] = pointId%dims[0];
    // idx[1] = pointId/dims[0];
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
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
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


class TriangleList
{
   public:
                   TriangleList() { maxTriangles = 1000000; triangleIdx = 0; pts = new float[9*maxTriangles]; };
     virtual      ~TriangleList() { delete [] pts; };

     void          AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxTriangles;
     int           triangleIdx;
};

void
TriangleList::AddTriangle(float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3)
{
    pts[9*triangleIdx+0] = X1;
    pts[9*triangleIdx+1] = Y1;
    pts[9*triangleIdx+2] = Z1;
    pts[9*triangleIdx+3] = X2;
    pts[9*triangleIdx+4] = Y2;
    pts[9*triangleIdx+5] = Z2;
    pts[9*triangleIdx+6] = X3;
    pts[9*triangleIdx+7] = Y3;
    pts[9*triangleIdx+8] = Z3;
    triangleIdx++;
}

vtkPolyData *
TriangleList::MakePolyData(void)
{
    int ntriangles = triangleIdx;
    int numPoints = 3*(ntriangles);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *tris = vtkCellArray::New();
    tris->EstimateSize(numPoints,4);
    for (int i = 0 ; i < ntriangles ; i++)
    {
        double pt[3];
        pt[0] = pts[9*i];
        pt[1] = pts[9*i+1];
        pt[2] = pts[9*i+2];
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[9*i+3];
        pt[1] = pts[9*i+4];
        pt[2] = pts[9*i+5];
        vtk_pts->SetPoint(ptIdx+1, pt);
        pt[0] = pts[9*i+6];
        pt[1] = pts[9*i+7];
        pt[2] = pts[9*i+8];
        vtk_pts->SetPoint(ptIdx+2, pt);
        vtkIdType ids[3] = { ptIdx, ptIdx+1, ptIdx+2 };
        tris->InsertNextCell(3, ids);
        ptIdx += 3;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetPolys(tris);
    tris->Delete();
    vtk_pts->Delete();

    return pd;
}

// ****************************************************************************
//  Function: IdentifyCase
//
//  Purpose: 
//     Figure out which case (0 - 255) we're dealing with for a specific cell.
//
//  Arguments:
//       isovalue (input):          the isovalue we're comparing field values to (provided for us in the project handout)
//       fieldValues (input):       an array of floats with 8 entries representing the field values at each vertex of our cell
//
//   Returns: the case number.
//
// ****************************************************************************
int IdentifyCase(float isovalue, float *fieldValues) {
    // Initialize case number to 0
    int icase = 0;

    // If the field value at vertex i is greater than the isovalue, add 2^i to the case number
    for (int i = 0; i < 8; i++) {
        if (fieldValues[i] > isovalue) {
            icase += (int) pow(2.0, i);
        }
    }

    return icase;
}

// ****************************************************************************
//  Function: interpolatePoints
//
//  Purpose: 
//     Interpolate to find the x, y, and z coordinates of all 3 vertices of a triangle
//
//  Arguments:
//       Xt (output):               a float array to store the x coordinates that we interpolated for
//       Yt (output):               a float array to store the y coordinates that we interpolated for
//       Zt (output):               a float array to store the z coordinates that we interpolated for
//       isovalue (input):          the isovalue we're using for interpolation
//       fieldValues (input):       an array of floats with 8 entries representing the field values at each vertex of our cell
//       vertices (input):          a 2D array storing the logical point indices of each vertex.
//       edge (input):      an int representing the edge we need to interpolate on
//       pointOnTriangle:   an int representing which point on the triangle we're working on interpolating for; helps index into Xt, Yt, and Zt appropriately
//       X:                 an array containing the X locations of our mesh
//       Y:                 an array containing the Y locations of our mesh
//       Z:                 an array containing the Z locations of our mesh
//
// ****************************************************************************
void interpolatePoints(float *Xt, float *Yt, float *Zt, float isovalue, float *fieldValues, int vertices[8][3], int edge, int pointOnTriangle, float *X, float *Y, float *Z) {
    switch (edge) {
        // Edge between v0 and v1
        case 0:
            Xt[pointOnTriangle] = interpolate(isovalue, fieldValues[0], fieldValues[1], X[vertices[0][0]], X[vertices[1][0]]);
            Yt[pointOnTriangle] = Y[vertices[0][1]];
            Zt[pointOnTriangle] = Z[vertices[0][2]];
            break;
        // Edge between v1 and v2
        case 1:
            Xt[pointOnTriangle] = X[vertices[1][0]];
            Yt[pointOnTriangle] = Y[vertices[1][1]];
            Zt[pointOnTriangle] = interpolate(isovalue, fieldValues[1], fieldValues[2], Z[vertices[1][2]], Z[vertices[2][2]]);
            break;
        // Edge between v2 and v3
        case 2:
            Xt[pointOnTriangle] = interpolate(isovalue, fieldValues[2], fieldValues[3], X[vertices[2][0]], X[vertices[3][0]]);
            Yt[pointOnTriangle] = Y[vertices[2][1]];
            Zt[pointOnTriangle] = Z[vertices[2][2]];
            break;
        // Edge between v0 and v3
        case 3:
            Xt[pointOnTriangle] = X[vertices[0][0]];
            Yt[pointOnTriangle] = Y[vertices[0][1]];
            Zt[pointOnTriangle] = interpolate(isovalue, fieldValues[0], fieldValues[3], Z[vertices[0][2]], Z[vertices[3][2]]);
            break;
        // Edge between v4 and v5
        case 4:
            Xt[pointOnTriangle] = interpolate(isovalue, fieldValues[4], fieldValues[5], X[vertices[4][0]], X[vertices[5][0]]);
            Yt[pointOnTriangle] = Y[vertices[4][1]];
            Zt[pointOnTriangle] = Z[vertices[4][2]];
            break;
        // Edge between v5 and v6
        case 5:
            Xt[pointOnTriangle] = X[vertices[5][0]];
            Yt[pointOnTriangle] = Y[vertices[5][1]];
            Zt[pointOnTriangle] = interpolate(isovalue, fieldValues[5], fieldValues[6], Z[vertices[5][2]], Z[vertices[6][2]]);
            break;
        // Edge between v6 and v7
        case 6:
            Xt[pointOnTriangle] = interpolate(isovalue, fieldValues[6], fieldValues[7], X[vertices[6][0]], X[vertices[7][0]]);
            Yt[pointOnTriangle] = Y[vertices[6][1]];
            Zt[pointOnTriangle] = Z[vertices[6][2]];
            break;
        // Edge between v4 and v7
        case 7:
            Xt[pointOnTriangle] = X[vertices[4][0]];
            Yt[pointOnTriangle] = Y[vertices[4][1]];
            Zt[pointOnTriangle] = interpolate(isovalue, fieldValues[4], fieldValues[7], Z[vertices[4][2]], Z[vertices[7][2]]);
            break;
        // Edge between v0 and v4
        case 8:
            Xt[pointOnTriangle] = X[vertices[0][0]];
            Yt[pointOnTriangle] = interpolate(isovalue, fieldValues[0], fieldValues[4], Y[vertices[0][1]], Y[vertices[4][1]]);
            Zt[pointOnTriangle] = Z[vertices[0][2]];
            break;
        // Edge between v1 and v5
        case 9:
            Xt[pointOnTriangle] = X[vertices[1][0]];
            Yt[pointOnTriangle] = interpolate(isovalue, fieldValues[1], fieldValues[5], Y[vertices[1][1]], Y[vertices[5][1]]);
            Zt[pointOnTriangle] = Z[vertices[1][2]];
            break;
        // Edge between v3 and v7
        case 10:
            Xt[pointOnTriangle] = X[vertices[3][0]];
            Yt[pointOnTriangle] = interpolate(isovalue, fieldValues[3], fieldValues[7], Y[vertices[3][1]], Y[vertices[7][1]]);
            Zt[pointOnTriangle] = Z[vertices[3][2]];
            break;
        // Edge between v2 and v6
        case 11:
            Xt[pointOnTriangle] = X[vertices[2][0]];
            Yt[pointOnTriangle] = interpolate(isovalue, fieldValues[2], fieldValues[6], Y[vertices[2][1]], Y[vertices[6][1]]);
            Zt[pointOnTriangle] = Z[vertices[2][2]];
            break;
        // Should never be here
        default:
            break;
    }
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("isosurfaces.vtk");
    rdr->Update();
    if (rdr->GetOutput() == NULL || rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Could not find input file." << endl;
        exit(EXIT_FAILURE);
    }

    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    int dims[3];
    rgrid->GetDimensions(dims);
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    // These were useful to me
    int edgeToVertex[12][2] =
        {
            {  0,  1 },
            {  2,  1 },
            {  2,  3 },
            {  0,  3 },
            {  4,  5 },
            {  5,  6 },
            {  6,  7 },
            {  4,  7 },
            {  0,  4 },
            {  1,  5 },
            {  3,  7 },
            {  2,  6 }
         };
    // This follows the convention in Lecture 11 slides (and project 6)
    // X is left-to-right, Y is up-and-down, Z is front-and-back.
    int offsetsI[8] = { 0, 1, 1, 0, 0, 1, 1, 0 };
    int offsetsJ[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
    int offsetsK[8] = { 0, 0, 1, 1, 0, 0, 1, 1 };

    TriangleList tl;
    int ncells = rgrid->GetNumberOfCells();
    cerr << "Number of cells to isosurface is " << ncells << endl;
    float isoval = 3.2;
    for (int i = 0 ; i < ncells ; i++)
    {
        // Debug
        if (DEBUG && (i != 4771)) continue;

        /*
            STEP 1:
            (1) collect all of the info about a cell (8 vertex locations, 8 field values) and put them in arrays you can use
        */

        // Get the logical cell index so we can find all the vertices
        int logicalCellIdx[3];
        GetLogicalCellIndex(logicalCellIdx, i, dims);

        // 8 vertices, 3 coordinate points each
        int vertices[8][3] = 
        {
            {logicalCellIdx[0], logicalCellIdx[1], logicalCellIdx[2]}, // v0
            {logicalCellIdx[0] + 1, logicalCellIdx[1], logicalCellIdx[2]}, // v1
            {logicalCellIdx[0] + 1, logicalCellIdx[1], logicalCellIdx[2] + 1}, // v2
            {logicalCellIdx[0], logicalCellIdx[1], logicalCellIdx[2] + 1}, // v3
            {logicalCellIdx[0], logicalCellIdx[1] + 1, logicalCellIdx[2]}, // v4
            {logicalCellIdx[0] + 1, logicalCellIdx[1] + 1, logicalCellIdx[2]}, // v5
            {logicalCellIdx[0] + 1, logicalCellIdx[1] + 1, logicalCellIdx[2] + 1}, // v6
            {logicalCellIdx[0], logicalCellIdx[1] + 1, logicalCellIdx[2] + 1} // v7
        };

        // Store the field values at every vertex
        float fieldValues[8];
        for (int j = 0; j < 8; j++) {
            // Get the point index first so we can index into F
            int ptIdx = GetPointIndex(vertices[j], dims);
            fieldValues[j] = F[ptIdx];
        }

        // Debug
        if (DEBUG) {
            cout << "Cell's 8 points are" << endl;
            for (int j = 0; j < 8; j++) {
                cout << "\tPt[" << j << "] = " << vertices[j][0] << ", " << vertices[j][1] << ", " << vertices[j][2] << endl;
                int ptIdx = GetPointIndex(vertices[j], dims);
                cout << "\t ptIdx = " << ptIdx << ", field = " << fieldValues[j] << ", loc = (" << X[vertices[j][0]] << ", " << Y[vertices[j][1]] << ", " << Z[vertices[j][2]] << ")" << endl;
            }
        }

        /*
            STEP 2:
            (2) figure out which case in the lookup table this cell goes to 
        */

        // Get the case we're looking for
        int icase = IdentifyCase(isoval, fieldValues);

        // Debug
        if (DEBUG) cout << "Triangle case index is " << icase << endl;

        /*
            STEP 3:
            (3) determine the triangles within the cell using the lookup table and info about the cell
                 -- hints:
                     check to make sure T is between 0 and 1.  call abort() if not and debug.
                     you will need to calculate the exact position of each triangle vertex along an edge.
                     I put my three vertices in "float Xt[3], Yt[3], Zt[3];"
                     And then called:             tl.AddTriangle(Xt[0], Yt[0], Zt[0], Xt[1], Yt[1], Zt[1], Xt[2], Yt[2], Zt[2]);
        */

        // Loop through the possible edges for our case in the lookup table; increment by 3 since 3 edges go together to make a triangle
        for (int j = 0; j < 16; j += 3) {
            // If the value we're looking at in the lookup table is -1 then there are no more triangles to add for our case
            if (lookupTable[icase][j] == -1) break;

            // Initialize variables to store our x, y, and z values for the triangle we're adding
            float Xt[3];
            float Yt[3];
            float Zt[3];

            // We'll use this as an index into Xt, Yt, and Zt so we know which point we're storing data for
            int trianglePt1 = j % 3;
            int trianglePt2 = (j + 1) % 3;
            int trianglePt3 = (j + 2) % 3;

            // Find what edges the triangle touches
            int edge1 = lookupTable[icase][j];
            int edge2 = lookupTable[icase][j + 1];
            int edge3 = lookupTable[icase][j + 2];

            if (DEBUG) {
                cout << "\t\tWorking on triangle with vertices along edges " << edge1 << ", " << edge2 << ", " << edge3 << endl;
            }

            // Interpolate the x, y, z coordinates for all 3 vertices of the triangle
            interpolatePoints(Xt, Yt, Zt, isoval, fieldValues, vertices, edge1, trianglePt1, X, Y, Z);
            interpolatePoints(Xt, Yt, Zt, isoval, fieldValues, vertices, edge2, trianglePt2, X, Y, Z);
            interpolatePoints(Xt, Yt, Zt, isoval, fieldValues, vertices, edge3, trianglePt3, X, Y, Z);

            // Add the triangle to our image
            tl.AddTriangle(Xt[0],  Yt[0], Zt[0], Xt[1], Yt[1], Zt[1], Xt[2], Yt[2], Zt[2]);
        }

        // My print statements for cell 4771:
/*
	Cell Log Idx = 18, 48, 1
	Cell's 8 points are
		Pt[0] = 18, 48, 1
		 ptIdx = 4918, field = 3.29052, loc = (-2.65306, 9.59184, -9.59184)
		Pt[1] = 19, 48, 1
		 ptIdx = 4919, field = 3.18607, loc = (-2.2449, 9.59184, -9.59184)
		Pt[2] = 19, 48, 2
		 ptIdx = 7419, field = 3.16074, loc = (-2.2449, 9.59184, -9.18367)
		Pt[3] = 18, 48, 2
		 ptIdx = 7418, field = 3.25718, loc = (-2.65306, 9.59184, -9.18367)
		Pt[4] = 18, 49, 1
		 ptIdx = 4968, field = 3.10792, loc = (-2.65306, 10, -9.59184)
		Pt[5] = 19, 49, 1
		 ptIdx = 4969, field = 3.02942, loc = (-2.2449, 10, -9.59184)
		Pt[6] = 19, 49, 2
		 ptIdx = 7469, field = 3.02574, loc = (-2.2449, 10, -9.18367)
		Pt[7] = 18, 49, 2
		 ptIdx = 7468, field = 3.10095, loc = (-2.65306, 10, -9.18367)
	Triangle case index is 9
		 Working on triangle with vertices along edges 0, 10, 8
			 Edge 0 has vertices at 0 and 1, and the new vertex should be placed at t=0.866634
			 Interpolated to get point -2.29933, 9.59184, -9.59184
			 Edge 10 has vertices at 3 and 7, and the new vertex should be placed at t=0.365998
			 Interpolated to get point -2.65306, 9.74123, -9.18367
			 Edge 8 has vertices at 0 and 4, and the new vertex should be placed at t=0.495728
			 Interpolated to get point -2.65306, 9.79418, -9.59184
		 Working on triangle with vertices along edges 10, 2, 0
			 Edge 10 has vertices at 3 and 7, and the new vertex should be placed at t=0.365998
			 Interpolated to get point -2.65306, 9.74123, -9.18367
			 Edge 2 has vertices at 2 and 3, and the new vertex should be placed at t=0.407094
			 Interpolated to get point -2.41106, 9.59184, -9.18367
			 Edge 0 has vertices at 0 and 1, and the new vertex should be placed at t=0.866634
			 Interpolated to get point -2.29933, 9.59184, -9.59184
 */
    }

    vtkPolyData *pd = tl.MakePolyData();

/*
    //This can be useful for debugging
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("proj6_out.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkCleanPolyData *cpd = vtkCleanPolyData::New();
    cpd->SetInputData(pd);
    cpd->SetAbsoluteTolerance(0);
    cpd->PointMergingOn();
    cpd->Update();
    vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
    pdn->SetInputData(cpd->GetOutput());
    //pdn->SetInputData(pd);
    pdn->Update();

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pdn->GetOutput());
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

    ren1->GetActiveCamera()->SetFocalPoint(0., 0., 0.);
    ren1->GetActiveCamera()->SetPosition(0,0,-62);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(1, 100);
    //ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
