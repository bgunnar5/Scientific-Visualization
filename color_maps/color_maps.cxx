#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


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

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F)
{
    // Check for invalid point
    if (pt[0] > X[dims[0] - 1] || pt[0] < X[0]) {
        return 0;
    }
    if (pt[1] > Y[dims[1] - 1] || pt[1] < Y[1]) {
        return 0;
    }

    // Placeholders for the bottom left indices of the cell that the point is located at
    int logicalPointIdx[2];

    // Locate which cell the x value is in
    for (int i = 1; i < dims[0]; i++) {
        if (pt[0] > X[i-1] && pt[0] < X[i]) {
            // Set x of the logical point index to be the idx at the bottom left of the cell
            logicalPointIdx[0] = i-1;
            break;
        }
    }

    // Locate which cell the y value is in
    for (int i = 1; i < dims[1]; i++) {
        if (pt[1] > Y[i-1] && pt[1] < Y[i]) {
            // Set y of the logical point inndex to be the idx at the bottom left of the cell
            logicalPointIdx[1] = i-1;
            break;
        }
    }

    int pointIdx = GetPointIndex(logicalPointIdx, dims);
    
    // Grab the scalar field values at each vertex
    float bottomLeft = F[pointIdx];
    float bottomRight = F[pointIdx + 1];
    // Adding by dims[0] is the same as going up a column
    float topLeft = F[pointIdx + dims[0]];
    float topRight = F[pointIdx + dims[0] + 1];

    // NOTE: I factored out interpolation into it's own function so I could use it multiple times
    // Interpolate to find the top and bottom values that are in line with the point we're looking for
    float bottomInterpolatedValue = interpolate(pt[0], X[logicalPointIdx[0]], X[logicalPointIdx[0] + 1], bottomLeft, bottomRight);
    float topInterpolatedValue = interpolate(pt[0], X[logicalPointIdx[0]], X[logicalPointIdx[0] + 1], topLeft, topRight);

    // Interpolate using the newfound top/bottom values and find the value of the point we're looking for
    float finalInterpolatedValue = interpolate(pt[1], Y[logicalPointIdx[1]], Y[logicalPointIdx[1] + 1], bottomInterpolatedValue, topInterpolatedValue);

    return finalInterpolatedValue;
}


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,0.5) 
//        F=1: (1.0,1.0,1.0) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
//  Note: color map is in float (0 to 1) and output is in unsigned char
//        (0 to 255), so the floating point values you calculate for 
//        the three color channels need to be multiplied by 255 when
//        being assigned to RGB.
//
// ****************************************************************************

void ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    // FA & FB for interpolation are:
    // R is between 0.0 and 1.0
    // G is betweenn 0.0 and 1.0
    // B is between 0.5 and 1.0

    // Colors are between 0 and 1 so use that for your A and B values in the interpolation
    float r = interpolate(F, 0.0, 1.0, 0.0, 1.0);
    float g = interpolate(F, 0.0, 1.0, 0.0, 1.0);
    float b = interpolate(F, 0.0, 1.0, 0.5, 1.0);

    RGB[0] = r * 255;
    RGB[1] = g * 255;
    RGB[2] = b * 255;
}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color 
//     using a difference colormap.
//
//     The difference color map has:
//        F=0:   (0.0, 0.0, 0.5) 
//        F=0.5: (1.0, 1.0, 1.0) 
//        F=1:   (0.5, 0.0, 0.0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
//  Note: color map is in float (0 to 1) and output is in unsigned char
//        (0 to 255), so the floating point values you calculate for 
//        the three color channels need to be multiplied by 255 when
//        being assigned to RGB.
//
// ****************************************************************************

void ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    float r, g, b;
    // FA & FB values are as follows:
    // R is between 0.0 and 1.0 or 1.0 and 0.5 respectively
    // G is betweenn 0.0 and 1.0 or 1.0 and 0.0 respectively
    // B is between 0.5 and 1.0 or 1.0 and 0.0 respectively

    // A and B values are 0 & 0.5 or 0.5 & 1 since that's where F's color values change
    if (F <= 0.5) {
        r = interpolate(F, 0.0, 0.5, 0.0, 1.0);
        g = interpolate(F, 0.0, 0.5, 0.0, 1.0);
        b = interpolate(F, 0.0, 0.5, 0.5, 1.0);
    }
    else {
        r = interpolate(F, 0.5, 1.0, 1.0, 0.5);
        g = interpolate(F, 0.5, 1.0, 1.0, 0.0);
        b = interpolate(F, 0.5, 1.0, 1.0, 0.0);
    }

    RGB[0] = r * 255;
    RGB[1] = g * 255;
    RGB[2] = b * 255;
}

// ****************************************************************************
//  Function: hsvToRGB
//
//  Purpose: 
//     Converts from an HSV color space to an RGB color space
//
//  Arguments:
//       hue (input):           a float between 0 and 360 degrees
//       saturation (input):    a float between 0.0 and 1.0
//       value (input):         a float between 0.0 and 1.0
//       rgb (output):          the location to store the color
//      
//  Note: rgb float values are not converted to unsigned chars yet since this
//        is just a helper function
//
// ****************************************************************************
void hsvToRGB(float hue, float saturation, float value, float *rgb) {
    // achromatic (grey)
    if (saturation == 0) {
        rgb[0] = rgb[1] = rgb[2] = value;
    }
    else {
        hue /= 60.f;
        // sector 0 to 5
        int i = floor(hue);
        float f = hue - i;
        // factorial part of h
        float p = value * (1 - saturation);
        float q = value * (1 - saturation * f);
        float t = value * (1 - saturation * (1 - f));
        switch (i) {
            case 0: rgb[0] = value; rgb[1] = t;     rgb[2] = p;     break;
            case 1: rgb[0] = q;     rgb[1] = value; rgb[2] = p;     break;
            case 2: rgb[0] = p;     rgb[1] = value; rgb[2] = t;     break;
            case 3: rgb[0] = p;     rgb[1] = q;     rgb[2] = value; break;
            case 4: rgb[0] = t;     rgb[1] = p;     rgb[2] = value; break;
            case 5: rgb[0] = value; rgb[1] = p;     rgb[2] = q;     break;
        }
    }
}

// ****************************************************************************
//  Function: ApplyHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0 <= F <= 1) to a color using 
//     an HSV rainbow colormap.
//
//     The rainbow colormap uses a saturation = 1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees.
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1  
//       RGB (output):  the location to store the color
//      
//  Note: as with the first two functions, make sure to multiple by 255 
//        when converting floats to unsigned chars.
//
// ****************************************************************************

void ApplyHSVColorMap(float F, unsigned char *RGB)  
{
    // F value is between 0 and 1 so that's our A and B values
    // Hue is between 0 and 360 degrees so that's our FA and FB values
    float hue = interpolate(F, 0.0, 1.0, 0, 360);

    // Convert the hsv value to RGB
    float rgb[3];
    hsvToRGB(hue, 1.0, 1.0, rgb);

    RGB[0] = rgb[0] * 255;
    RGB[1] = rgb[1] * 255;
    RGB[2] = rgb[2] * 255;
}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("color_maps_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3 ; i++)
       for (j = 0 ; j < 3*nx*ny ; j++)
            buffer[i][j] = 0;

    for (i = 0 ; i < nx ; i++) {
        // m value in y=mx+b equation; 18 since we go from x=-9 to x=9 and nx-1 for 0 based indexing
        float m_for_x = 18.0/(nx-1);

        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            // m value in y=mx+b equation; 18 since we go from x=-9 to x=9 and nx-1 for 0 based indexing
            float m_for_y = 18.0/(ny-1);

            pt[0] = m_for_x * i - 9;
            pt[1] = m_for_y * j - 9;

            // Get the field value at the determined point
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);

            // Normalize the field value with the normalization equation:
            // x_norm = (x - x_min) / (x_max - x_min)
            // Here our range is from 1.2 to 5.02 so x_min = 1.2 and x_max = 5.02
            float normalizedF = (f - 1.2) / (5.02 - 1.2);
            
            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }
    }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
