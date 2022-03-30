/*
Volume Rendering
Author: Brian Gunnarson
Date Created: March 1, 2022

This is my final project for CIS 410 Scientific Visualization at the University of Oregon.
In this project I am implementing volume rendering on a data set that was provided by Hank Childs.
This volume rendering is done through ray casting. Early ray termination is implemented and defaults
to terminate when the opacity value of a pixel hits 99.9%.

To obtain a low resolution image, change the HIGH_RES value to 0.
To change the size of the image, change the WIDTH & HEIGHT constants.
To change the number of samples per ray, change the NUMSAMPLES constant.
To change the ray termination criteria (i.e. the max opaqueness percentage of a pixel), change the TERMINATION_CRITERIA constant.

Note: Any code that is not mine I have assumed to be Prof. Hank Childs and is credited as such.
*/
#include <iostream>
#include <math.h>
#include <time.h>

#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>

using namespace std;

#define HIGH_RES 1
#define DEBUG 0

const int WIDTH = 1000;
const int HEIGHT = 1000;
const int NUMSAMPLES = 1024;
const double TERMINATION_CRITERIA = 0.999;

// Camera struct and setup function provided by Hank Childs
struct Camera
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};

// Transfer function struct and setup function provided by Hank Childs
struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    double         *opacities; // size is numBins

    // Take in a value and applies the transfer function.
    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
    {
        // cout << "Applying transfer function to: " << value << endl;
        // Figure out which bin our value is in
        int bin = GetBin(value);

        // Invalid bin
        if (bin == -1)
        {
            RGB[0] = (unsigned char) 0;
            RGB[1] = (unsigned char) 0;
            RGB[2] = (unsigned char) 0;
            opacity = 0.0;
        }
        // Valid bin so get the color and opacity
        else
        {
            RGB[0] = colors[3*bin+0];
            RGB[1] = colors[3*bin+1];
            RGB[2] = colors[3*bin+2];
            opacity = opacities[bin];
        }

        // cout << "RGB: " << (int) RGB[0] << ", " << (int) RGB[1] << ", " << (int) RGB[2] << endl;
        // cout << "Opacity: " << opacity << endl;
    }

    // Figure out which bin "value" lies in.
    // If "min" is 2 and "max" is 4, and there are 10 bins, then
    //   bin 0 = 2->2.2
    //   bin 1 = 2.2->2.4
    //   bin 2 = 2.4->2.6
    //   bin 3 = 2.6->2.8
    //   bin 4 = 2.8->3.0
    //   bin 5 = 3.0->3.2
    //   bin 6 = 3.2->3.4
    //   bin 7 = 3.4->3.6
    //   bin 8 = 3.6->3.8
    //   bin 9 = 3.8->4.0
    // and, for example, a "value" of 3.15 would return the color in bin 5
    // and the opacity at "opacities[5]".
    int GetBin(double value)
    {
        // Find the size of each bin
        int diff = max - min;
        double binSize = (double) diff / numBins;

        // Start at the min
        double currValue = (double) min;
        for (int i = 0; i < numBins; i++)
        {
            // If value is between the min and max of the current bin then that's the bin we return
            if ((value >= currValue) && (value < (currValue + binSize)))
            {
                // cout << "Mapped to bin: " << i << endl;
                return i;
            }

            currValue += binSize;
        }

        // cout << "Out of range: no valid bin" << endl;
        // Value is not in a valid bin
        return -1;
    }
};

TransferFunction
SetupTransferFunction(void)
{
    int  i;

    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    rv.opacities = new double[256];
    unsigned char charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 
            4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
            13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 
            11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 
            2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 
            11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 
            17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 
            9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 
            4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 
            18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 
            44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 
            67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 
            64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 
            38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 
            28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 
            105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 
            205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 
            208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 
            109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 
            30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 
            4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0, 0, 0, 0, 0, 0
        };

    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = charOpacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = { 
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76 
       };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0 ; i < numControlPoints-1 ; i++)
    {
        int start = controlPointPositions[i]*rv.numBins;
        int end   = controlPointPositions[i+1]*rv.numBins+1;
        // cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << endl;
        if (end >= rv.numBins)
            end = rv.numBins-1;
        for (int j = start ; j <= end ; j++)
        {
            double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (int k = 0 ; k < 3 ; k++)
                rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                                 + controlPointColors[3*i+k];
        }
    }    

    return rv;
}

Camera
SetupCamera(void)
{
    Camera rv;
    rv.focus[0] = 0;
    rv.focus[1] = 0;
    rv.focus[2] = 0;
    rv.up[0] = 0;
    rv.up[1] = -1;
    rv.up[2] = 0;
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;

    return rv;
}

// ****************************************************************************
//  Function: crossProduct
//
//  Purpose: 
//     A helper function that takes the cross product of 2 vectors
//
//  Arguments:
//       a:     an array of doubles representing a vector
//                      i.e. index 0 = x, index 1 = y, index 2 = z
//       b:     an array of doubles representing a vector
//                      i.e. index 0 = x, index 1 = y, index 2 = z
//       arr (output):   an array of doubles to store a x b
//
// ****************************************************************************
void crossProduct(double *a, double *b, double (&arr)[3])
{
    // A.y * B.z - A.z * B.y
    arr[0] = (a[1] * b[2]) - (a[2] * b[1]);
    // A.z * B.x - A.x * B.z
    arr[1] = (a[2] * b[0]) - (a[0] * b[2]);
    // A.x * B.y - A.y * B.x
    arr[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

// ****************************************************************************
//  Function: magnitude
//
//  Arguments:
//       vector:    an array of doubles representing a vector
//                      i.e. index 0 = x, index 1 = y, index 2 = z
//
//  Returns: the magnitude value of a vector.
//
// ****************************************************************************
double magnitude(double *vector)
{
    double x_squared = vector[0] * vector[0];
    double y_squared = vector[1] * vector[1];
    double z_squared = vector[2] * vector[2];
    double result = sqrt(x_squared + y_squared + z_squared);
    return result;
}

// ****************************************************************************
//  Function: GetPointIndex
//  Author: Hank Childs
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
    // return idx[1]*dims[0]+idx[0];
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
//     pt: a three-dimensional location
//     dims: an array of size three.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y, and the third the size of Z.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     Z: an array (size is specified by dims).  
//              This contains the Z locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1]*dims[2].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************
double EvaluateFieldAtLocation(const double *pt, const int *dims, const float *X, 
                              const float *Y, const float *Z, const float *F)
{
    // Check for invalid point
    if (pt[0] > X[dims[0] - 1] || pt[0] < X[0]) {
        return 0;
    }
    if (pt[1] > Y[dims[1] - 1] || pt[1] < Y[0]) {
        return 0;
    }
    if (pt[2] > Y[dims[2] - 1] || pt[2] < Z[0]) {
        return 0;
    }

    int logicalPointIdxV0[3];

    // Locate which cell the x value is in
    for (int i = 1; i < dims[0]; i++) {
        if (pt[0] > X[i-1] && pt[0] < X[i]) {
            // Set x of the logical point index to be the idx at the bottom left of the cell
            logicalPointIdxV0[0] = i-1;
            break;
        }
    }

    // Locate which cell the y value is in
    for (int i = 1; i < dims[1]; i++) {
        if (pt[1] > Y[i-1] && pt[1] < Y[i]) {
            // Set y of the logical point index to be the idx at the bottom left of the cell
            logicalPointIdxV0[1] = i-1;
            break;
        }
    }

    // Locate which cell the z value is in
    for (int i = 1; i < dims[2]; i++) {
        if (pt[2] > Z[i-1] && pt[2] < Z[i]) {
            // Set z of the logical point index to be the idx at the bottom left of the cell
            logicalPointIdxV0[2] = i-1;
            break;
        }
    }

    // From the logical point index at vertex 0 we can obtain the logical point indices of the remaining vertices
    int logicalPointIdxV1[3] = {logicalPointIdxV0[0] + 1, logicalPointIdxV0[1]    , logicalPointIdxV0[2]    };
    int logicalPointIdxV2[3] = {logicalPointIdxV0[0] + 1, logicalPointIdxV0[1]    , logicalPointIdxV0[2] + 1};
    int logicalPointIdxV3[3] = {logicalPointIdxV0[0]    , logicalPointIdxV0[1]    , logicalPointIdxV0[2] + 1};
    int logicalPointIdxV4[3] = {logicalPointIdxV0[0]    , logicalPointIdxV0[1] + 1, logicalPointIdxV0[2]    };
    int logicalPointIdxV5[3] = {logicalPointIdxV0[0] + 1, logicalPointIdxV0[1] + 1, logicalPointIdxV0[2]    };
    int logicalPointIdxV6[3] = {logicalPointIdxV0[0] + 1, logicalPointIdxV0[1] + 1, logicalPointIdxV0[2] + 1};
    int logicalPointIdxV7[3] = {logicalPointIdxV0[0]    , logicalPointIdxV0[1] + 1, logicalPointIdxV0[2] + 1};

    // Use the logical point indices to get the point index of each vertex
    int v0PointIdx = GetPointIndex(logicalPointIdxV0, dims);
    int v1PointIdx = GetPointIndex(logicalPointIdxV1, dims);
    int v2PointIdx = GetPointIndex(logicalPointIdxV2, dims);
    int v3PointIdx = GetPointIndex(logicalPointIdxV3, dims);
    int v4PointIdx = GetPointIndex(logicalPointIdxV4, dims);
    int v5PointIdx = GetPointIndex(logicalPointIdxV5, dims);
    int v6PointIdx = GetPointIndex(logicalPointIdxV6, dims);
    int v7PointIdx = GetPointIndex(logicalPointIdxV7, dims);

    // Interpolate in x; need to interpolate between v0 & v1, v2 & v3, v4 & v5, and v6 & v7
    double val01 = interpolate(pt[0], X[logicalPointIdxV0[0]], X[logicalPointIdxV1[0]], F[v0PointIdx], F[v1PointIdx]);
    double val23 = interpolate(pt[0], X[logicalPointIdxV2[0]], X[logicalPointIdxV3[0]], F[v2PointIdx], F[v3PointIdx]);
    double val45 = interpolate(pt[0], X[logicalPointIdxV4[0]], X[logicalPointIdxV5[0]], F[v4PointIdx], F[v5PointIdx]);
    double val67 = interpolate(pt[0], X[logicalPointIdxV6[0]], X[logicalPointIdxV7[0]], F[v6PointIdx], F[v7PointIdx]);

    // Interpolate in z; need to interpolate between v0/v1 & v2/v3 and v4/v5 & v6/v7
    double val0123 = interpolate(pt[2], Z[logicalPointIdxV0[2]], Z[logicalPointIdxV3[2]], val01, val23);
    double val4567 = interpolate(pt[2], Z[logicalPointIdxV4[2]], Z[logicalPointIdxV7[2]], val45, val67);

    // Interpolate in y; put it all together
    double finalInterpolatedValue = interpolate(pt[1], Y[logicalPointIdxV0[1]], Y[logicalPointIdxV7[1]], val0123, val4567);

    return finalInterpolatedValue;
}

// ****************************************************************************
//  Function: compositeColors
//
//  Purpose: 
//     A helper function that composites 2 samples
//
//  Arguments:
//       r:             a double representing the current red color of a pixel
//       g:             a double representing the current green color of a pixel
//       b:             a double representing the current blue color of a pixel
//       BRGB:          an unsigned char array consisting of the back RGB values 
//                          that we're compositing with our current RGB values
//       opacity:       a double representing the current opacity of a pixel
//       bOpacity:      a double representing the opacity of the back sample
//                          we're compositing with
//
// ****************************************************************************
void compositeColors(double &r, double &g, double &b, unsigned char *BRGB, double &fOpacity, double bOpacity)
{   
    // This equation expects values to be between 0 and 1 so we divide the back RGB (BRGB) by 255 to get that
    // Removed the initial multiplication by front opacity since it'd be redundant with how I loop and calculate RGB
    r = r + (1 - fOpacity) * bOpacity * ((double) BRGB[0] / 255.0);
    g = g + (1 - fOpacity) * bOpacity * ((double) BRGB[1] / 255.0);
    b = b + (1 - fOpacity) * bOpacity * ((double) BRGB[2] / 255.0);
    fOpacity = fOpacity + (1 - fOpacity) * bOpacity;
}

// ****************************************************************************
//  Function: WriteImage
//  Author: Hank Childs
//
//  Purpose: 
//     A function that writes data to a vtk image.
//
//  Arguments:
//       img:       a vtkImageData object; the actual image we're making
//       filename:  a string representing the name of the file we're creating
//
//  Note: I obtained this code from a project I did in CIS 441 with Prof. Childs
//
// ****************************************************************************
void WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

int main()
{
    // Get the data from the provided vtk file
    vtkDataSetReader *reader = vtkDataSetReader::New();
    if (HIGH_RES)
        reader->SetFileName("astro512.vtk");
    else
        reader->SetFileName("astro64.vtk");
    reader->Update();

    // Set up the rectilinear grid and get the dimensions
    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) reader->GetOutput();
    rgrid->GetDimensions(dims);

    // Get the X, Y, and Z coordinates as well as the scalar data
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    // Initialize the camera and the transfer function
    Camera c = SetupCamera();
    TransferFunction tf = SetupTransferFunction();

    // Calculate the step size between each sample
    double stepSize = (c.far - c.near) / (NUMSAMPLES - 1.0);

    // Create an image specified by the number of pixels at the top of the file
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(WIDTH, HEIGHT, 1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    // Create the color buffer for this image (RGB value for each pixel so size = num pixels * 3); initialize every entry to be 0
    unsigned char *colorBuffer = (unsigned char *) image->GetScalarPointer(0, 0, 0);
    int bufferSize = WIDTH * HEIGHT * 3;
    for (int i = 0; i < bufferSize; i++)
        colorBuffer[i] = 0;

    // Start time for the rendering
    clock_t tStart = clock();

    /*
    *
    * STEP 1: FIND THE RAY FOR A GIVEN PIXEL
    *
    */
    double look[3], temp[3], u[3], v[3], x[3], y[3], tempNormalized, numerator, mulFactor;

    // Calculate look vector using equation: focus - position
    for (int k = 0; k < 3; k++)
        look[k] = c.focus[k] - c.position[k];

    // Calculate u vector for the ray finding equation
    crossProduct(look, c.up, temp);
    tempNormalized = magnitude(temp);
    for (int k = 0; k < 3; k++)
        u[k] = temp[k] / tempNormalized;

    // Calculate v vector for the ray finding equation
    crossProduct(look, u, temp);
    tempNormalized = magnitude(temp);
    for (int k = 0; k < 3; k++)
        v[k] = temp[k] / tempNormalized;

    // Calculate x and y vector for the ray finding equation
    numerator = 2 * tan((c.angle * (M_PI/180)) / 2);

    mulFactor = numerator / WIDTH;
    for (int k = 0; k < 3; k++)
        x[k] = mulFactor * u[k];
    
    mulFactor = numerator / HEIGHT;
    for (int k = 0; k < 3; k++)
        y[k] = mulFactor * v[k];

    // 1. look/||look||
    tempNormalized = magnitude(look);
    for (int k = 0; k < 3; k++)
        look[k] = look[k] / tempNormalized;

    for (int i = 0; i < WIDTH; i++) 
    {
        if (DEBUG && i != 50) continue;
        for (int j = 0; j < HEIGHT; j++) 
        {
            if (DEBUG && j != 50) continue;
            // Equation to find ray for a pixel: d(i, j) = look/||look|| + (((2i + 1 - W)/2)x) + ((2j + 1 - H)/2)y)
            double ray[3], xForThisPixel[3], yForThisPixel[3], coefficient;

            // 2. ((2i + 1 - W)/2)x
            coefficient = (2 * i + 1 - WIDTH) / 2.0;
            for (int k = 0; k < 3; k++)
                xForThisPixel[k] = coefficient * x[k];

            // 3. ((2j + 1 - H)/2)y
            coefficient = ((2 * j) + 1 - HEIGHT) / 2.0;
            for (int k = 0; k < 3; k++)
                yForThisPixel[k] = coefficient * y[k];

            // 4. d(i, j) <-- Combine 1, 2, and 3 from above
            for (int k = 0; k < 3; k++)
                ray[k] = look[k] + xForThisPixel[k] + yForThisPixel[k];
            
            /*
             * STEP 2: INTERSECT VOLUME WITH RAY
             * & STEP 3: CALCULATE COLOR FROM INTERSECTION
             * Used to be 2 separate steps and implemented as a 2 pass algorithm but I switched it to 1 pass to try to improve performance
             */
            // Initialize variables to store the field value at specific locations in the ray, distance from camera, current position, front opacity, back opacity, and the final RGB values
            double fieldVal, distance, currPosition[3], opacity, bOpacity, r, g, b;

            // Initialize array to store the RGB value at a certain field value
            unsigned char RGB[3];

            // Start at the near plane with our ray
            distance = c.near;

            // Calculate the offset so we know exactly how far to move in each direction
            double offset[3] = {ray[0] * distance, ray[1] * distance, ray[2] * distance};

            // Apply the offset to the position of the camera to get our point; think of it like adding the offset vector to the origin
            for (int l = 0; l < 3; l++)
                currPosition[l] = c.position[l] + offset[l];
            
            // Calculate the field value at the current position
            fieldVal = EvaluateFieldAtLocation(currPosition, dims, X, Y, Z, F);
            // cout << "Value at that position is " << fieldVal << endl;

            // Apply transfer function to the first value we get so we can get initial RGBA values to start with
            tf.ApplyTransferFunction(fieldVal, RGB, opacity);

            // Compositing equation expects variables to be between 0 and 1 so do that conversion here; also switching to doubles bc unsigned chars drive me insane
            r = (double) RGB[0] / 255.0;
            g = (double) RGB[1] / 255.0;
            b = (double) RGB[2] / 255.0;

            // Increase the distance we travel by our step size
            distance += stepSize;

            // Loop through the rest of the samples
            for (int k = 1; k < NUMSAMPLES; k++)
            {
                /*
                 *
                 * EARLY RAY TERMINATION: stop doing calculations when our opacity becomes opaque or (close to) 1
                 *
                 */
                if (opacity >= TERMINATION_CRITERIA) break;

                // Calculate the offset so we know exactly how far to move in each direction
                double offset[3] = {ray[0] * distance, ray[1] * distance, ray[2] * distance};

                // Apply the offset to the position of the camera to get our point; think of it like adding the offset vector to the origin
                for (int l = 0; l < 3; l++)
                    currPosition[l] = c.position[l] + offset[l];
                // cout << "Position for sample[" << k << "]: " << currPosition[0] << ", " << currPosition[1] << ", " << currPosition[2] << endl;

                // Calculate the value at the current position
                fieldVal = EvaluateFieldAtLocation(currPosition, dims, X, Y, Z, F);
                // cout << "Value at that position is " << fieldVal << endl;

                // Apply the transfer function to get the back RGB for the compositing equation
                tf.ApplyTransferFunction(fieldVal, RGB, bOpacity);

                // Optimization to only composite samples if volume intersects with ray (otherwise compositing won't change our r, g, b values at all)
                // Also, if the volume intersects but the back opacity is 0 then there will be no change to r,g,b either
                if ((RGB[0] != 0 || RGB[1] != 0 || RGB[2] != 0) && (bOpacity != 0)) 
                {
                    // Opacity correction
                    bOpacity = 1 - pow((1 - bOpacity), 500 / (double) NUMSAMPLES);

                    // Composite 2 samples to update our rgb values
                    compositeColors(r, g, b, RGB, opacity, bOpacity);
                }

                // Increase the distance we travel by our step size
                distance += stepSize;
            }

            /*
             *
             * STEP 4: ASSIGN COLOR TO PIXEL
             *
             */
            // Get the correct index into the buffer
            int pixelIndex = (j * WIDTH + i) * 3;

            // Assign the color at the specified index
            colorBuffer[pixelIndex + 0] = floor(r*255);
            colorBuffer[pixelIndex + 1] = floor(g*255);
            colorBuffer[pixelIndex + 2] = floor(b*255);

            if (DEBUG) cout << "Final color for the pixel is " << (int) colorBuffer[pixelIndex] << ", " << (int) colorBuffer[pixelIndex + 1] << ", " << (int) colorBuffer[pixelIndex + 2] << endl;
        }
    }

    // End time for the rendering
    cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << " seconds" << endl;

    // Create the image and delete the reader object
    WriteImage(image, "astro");
    reader->Delete();
}

