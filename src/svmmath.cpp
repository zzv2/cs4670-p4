/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #1:
 * svmmath.cpp
 *		a routine for intersecting >2 lines (for vanishing point
 *		computation);
 *		routines for computing the homography for the reference
 *		plane and arbitrary polygons
 **************************************************************/

#pragma warning(disable : 4996)

#include "Eigen/Core"
#include "MinEig.h"

#include "svmmath.h"
#include "jacob.h"
#include "vec.h"
#include <cstring>
#include <cstdio>
#include <assert.h>
#include <iostream>

using namespace Eigen;
using namespace std;

//
// TODO 1: BestFitIntersect()
//		Given lines, the list of 3 or more lines to be intersected,
//		find the best fit intersection point.
//		See http://www.cs.cornell.edu/courses/cs4670/2013fa/projects/p4/vanishing.txt
//      Note that the "numerical conditioning" part of this description is optional, but recommended
//	
SVMPoint BestFitIntersect(const std::list<SVMLine> &lines, int imgWidth, int imgHeight)
{
    // check
    if (lines.size() < 2)
	{
            fprintf(stderr, "Not enough lines to compute the best fit.");
            abort();
	}

    SVMPoint bestfit;
    list<SVMLine>::const_iterator iter;

    // To accumulate stuff
    typedef Matrix<double, Dynamic, 3, RowMajor> Matrix3;

    int numLines = (int) lines.size();
    Matrix3 A = Matrix3::Zero(numLines, 3);	

    // Transformation for numerical stability

    // Note: iterate through the lines list as follows:
    //		for (iter = lines.begin(); iter != lines.end(); iter++) {
    //			...iter is the pointer to the current line...
    //		}
    // Note: Function to find eigenvector with smallest eigenvalue of A^TA is MinEig(A, eval, evec)
    // This minimum eigenvector of A^T A is also the same as the minimum singular vector of A.

    //TODO-BLOCK-BEGIN
    int rowCount = 0;

    for (iter = lines.begin(); iter != lines.end(); iter++) {
        // ...iter is the pointer to the current line...

        // specify each line's endpoints e1 and e2 in homogeneous coordinates
        // e1 = (x1_i , y1_i, w)
        // e2 = (x2_i , y2_i, w)
        Vec3d e1 = Vec3d(iter->pnt1->u,iter->pnt1->v,1);
        Vec3d e2 = Vec3d(iter->pnt2->u,iter->pnt2->v,1);

        // compute a homogenous coordinate vector representing the line
        // as the cross product of its two endpoints
        // (a_i,b_i,c_i) = e1  X  e2 (there is a built in cross product available in vec.h)
        // note that this resulting vector is just the parameters of
        // the equation a_i x + b_i y + c_i = 0 of the 2D infinite line
        // passing through the two endpoints
        Vec3d l = cross(e1,e2);


        // arrange all of the lines into a nx3 matrix:
        A.row(rowCount) << l[0], l[1], l[2];

        rowCount++;
    }

    double eval, evec[3];
    MinEig(A, eval, evec);
    bestfit = SVMPoint(evec[0]/evec[2], evec[1]/evec[2]);

    //TODO-BLOCK-END
    /******** END TODO ********/
	
    return bestfit;
}


//
// TODO 2: ConvertToPlaneCoordinate()
//		Given a plane defined by points, converts their coordinates into
//		plane coordinates. See the following document for more detail.
//		http://www.cs.cornell.edu/courses/cs4670/2013fa/projects/p4/homography.pdf.
//      
//      points contains the three co-planer points which will define the desired plane
//      basisPts should be filled with the corresponding points, now in plane coordinates
//      The final divisors you apply to the u and v coordinates should be saved uScale and vScale
//
void ConvertToPlaneCoordinate(const vector<SVMPoint>& points, vector<Vec3d>& basisPts, double &uScale, double &vScale)
{
    int numPoints = points.size();

    /******** BEGIN TODO ********/
    //TODO-BLOCK-BEGIN
    Vec3d p = Vec3d(points[0].X, points[0].Y, points[0].Z);
    Vec3d q = Vec3d(points[1].X, points[1].Y, points[1].Z);
    Vec3d r = Vec3d(points[2].X, points[2].Y, points[2].Z);

    Vec3d ex = (p - r);
    ex.normalize();

    Vec3d s = ((q - r) * ex) * ex;
    Vec3d ey = (q - r) - s;
    ey.normalize();

    double umin, umax, vmin, vmax;

    for (int i = 0; i < numPoints; i++)
    {
        Vec3d origPt = Vec3d(points[i].X, points[i].Y, points[i].Z);
        Vec3d pt = Vec3d((origPt - r) * ex, (origPt - r) * ey, 1);
        basisPts.push_back(pt);

        umin = min(umin, pt[0]);
        umax = max(umax, pt[0]);
        vmin = min(vmin, pt[1]);
        vmax = max(vmax, pt[1]);
    }

    for (int i = 0; i < numPoints; i++)
    {
        basisPts[i][0] = (basisPts[i][0] - umin)/(umax - umin);
        basisPts[i][1] = (basisPts[i][1] - vmin)/(vmax - vmin);
    }

    uScale = 1.0/(umax - umin);
    vScale = 1.0/(vmax - vmin);
    //TODO-BLOCK-END
    /******** END TODO ********/
}



//
// TODO 3: ComputeHomography()
//		Computes the homography H from the plane specified by "points" to the image plane,
//		and its inverse Hinv.
//		If the plane is the reference plane (isRefPlane == true), don't convert the
//		coordinate system to the plane. Only do this for polygon patches where
//		texture mapping is necessary.
//		Coordinate system conversion is to be implemented in a separate routine
//		ConvertToPlaneCoordinate.
//		For more detailed explaination, see
//		http://www.cs.cornell.edu/courses/cs4670/2013fa/projects/p4/homography.pdf.
//
void ComputeHomography(CTransform3x3 &H, CTransform3x3 &Hinv, const vector<SVMPoint> &points, vector<Vec3d> &basisPts, bool isRefPlane)
{
    int i;
    int numPoints = (int) points.size();
    assert( numPoints >= 4 );

    basisPts.clear();
    if (isRefPlane) // reference plane
    {
        for (i=0; i < numPoints; i++) {
            Vec3d tmp = Vec3d(points[i].X, points[i].Y, points[i].W); // was Z, not W
            basisPts.push_back(tmp);
        }
    } 
    else // arbitrary polygon
    {
        double uScale, vScale; // unused in this function
        ConvertToPlaneCoordinate(points, basisPts, uScale, vScale);
    }

    // A: 2n x 9 matrix where n is the number of points on the plane
    //    as discussed in lecture
    int numRows = 2 * numPoints;
    const int numCols = 9;

    typedef Matrix<double, Dynamic, 9, RowMajor> MatrixType;
    MatrixType A = MatrixType::Zero(numRows, numCols);

    /******** BEGIN TODO ********/
    /* Fill in the A matrix for the call to MinEig */
    //TODO-BLOCK-BEGIN
    for (int i = 0; i < numPoints; i++)
    {
        SVMPoint origPt = points[i];
        Vec3d basisPt = basisPts[i];

        A.row(2*i) << basisPt[0], basisPt[1], basisPt[2], 0, 0, 0, -origPt.u * basisPt[0], -origPt.u * basisPt[1], -origPt.u * basisPt[2];
        A.row(2*i+1) << 0, 0, 0, basisPt[0], basisPt[1], basisPt[2], -origPt.v * basisPt[0], -origPt.v * basisPt[1], -origPt.v * basisPt[2];
    }
    //TODO-BLOCK-END

    double eval, h[9];
    MinEig(A, eval, h);

    H[0][0] = h[0];
    H[0][1] = h[1];
    H[0][2] = h[2];

    H[1][0] = h[3];
    H[1][1] = h[4];
    H[1][2] = h[5];

    H[2][0] = h[6];
    H[2][1] = h[7];
    H[2][2] = h[8];

    /******** END TODO ********/

    // compute inverse of H
    if (H.Determinant() == 0)
        fl_alert("Computed homography matrix is uninvertible \n");
    else
        Hinv = H.Inverse();

    int ii;
    printf("\nH=[\n");
    for (ii=0; ii<3; ii++)
        printf("%e\t%e\t%e;\n", H[ii][0]/H[2][2], H[ii][1]/H[2][2], H[ii][2]/H[2][2]);
    printf("]\nHinv=[\n");

    for (ii=0; ii<3; ii++)
        printf("%e\t%e\t%e;\n", Hinv[ii][0]/Hinv[2][2], Hinv[ii][1]/Hinv[2][2], Hinv[ii][2]/Hinv[2][2]);

    printf("]\n\n");
    fflush(stdout);
}

