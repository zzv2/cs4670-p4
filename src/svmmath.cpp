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


    double evec[3];
    double eval;

    Matrix3 ata = A.transpose() * A;
    cout << "A^t * A =" << endl << ata << endl;

    MinEig(ata, eval, evec);

    // the singular vector associated with the smallest singular value
    // is the best-fit vanishing point vector V (evec)
    evec[0] /= evec[2];
    evec[1] /= evec[2];

    bestfit = SVMPoint(evec[0],evec[1]);

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
    printf("TODO: %s:%d\n", __FILE__, __LINE__);
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
    printf("TODO: %s:%d\n", __FILE__, __LINE__);
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

