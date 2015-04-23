/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #3:
 * ImgViewBox.cpp
 *		routines for finding the corners of an axis-aligned
 *      box (in 2D and in 3D)
 **************************************************************/
#pragma warning(disable : 4996)

#include "ImgView.h"

//
// TODO 6: solveForOppositeCorners()
//     Given the 2D positions of two corners of a rectangular face parallel to the XZ plane, compute
//     the 2D positions of the other two corners
void ImgView::solveForOppositeCorners(double u0, double v0, double u2, double v2,
                                      double &u1, double &v1, double &u3, double &v3)
{
    /* Vanishing points must be known */
    assert(xVanish.known() && yVanish.known() && zVanish.known());    


    /******** BEGIN TODO ********/ 
    // Given the 2D positions of corners p0 and p2 of the face, compute the 2D positions of p1 and p3
    // Remember that this face is on a plane perpendicular to the plane x=0
    // Store the results in variables 'u1, v1' and 'u3, v3'

    //TODO-BLOCK-BEGIN

    // Get the vanishing points in homogenous image coordinates
    Vec3d vx = Vec3d(xVanish.u, xVanish.v, 1.0);
    Vec3d vy = Vec3d(yVanish.u, yVanish.v, 1.0);
    Vec3d vz = Vec3d(zVanish.u, zVanish.v, 1.0);

    // Get the image coords for the known corners
    Vec3d p0 = Vec3d(u0, v0, 1.0);
    Vec3d p2 = Vec3d(u2, v2, 1.0);

    // Line [p0, vz]
    Vec3d l1 = cross(p0, vz);
    l1[0] /= l1[2];
    l1[1] /= l1[2];
    l1[2] = 0.0;

    // Line [p2, vx]
    Vec3d l2 = cross(p2, vx);
    l2[0] /= l2[2];
    l2[1] /= l2[2];
    l2[2] = 0.0;

    // Line [p0, vx]
    Vec3d l3 = cross(p0, vx);
    l3[0] /= l3[2];
    l3[1] /= l3[2];
    l3[2] = 0.0;

    // Line [p2, vz]
    Vec3d l4 = cross(p2, vz);
    l4[0] /= l4[2];
    l4[1] /= l4[2];
    l4[2] = 0.0;

    // Find p3
    Vec3d p3 = cross(l1, l2);
    p3[0] /= p3[2];
    p3[2] /= p3[2];
    p3[3] = 0.0;

    // Find p1
    Vec3d p1 = cross(l3, l4);
    p1[0] /= p1[2];
    p1[1] /= p1[2];
    p1[2] = 0.0;

    // Set the output variables
    u1 = p1[0];
    v1 = p1[1];
    u3 = p3[0];
    v3 = p3[1];
    //TODO-BLOCK-END
    /********* END TODO ********/
}

//
// TODO 7: solveForOppositeFace()
//     Given the 2D positions of one rectangular face parallel to the XZ plane, 
//     compute the 2D positions of a parallel face being swept out from it.
//     The mouse position is given; one of the lines on the parallel face should pass
//     through the mouse position
void ImgView::solveForOppositeFace(SVMSweep *sweep, double imgX, double imgY,
                                   Vec3d &p4_out, Vec3d &p5_out, Vec3d &p6_out, Vec3d &p7_out)
{
    SVMPolygon *poly = sweep->poly;

    if (poly == NULL)
        return;

    // Get the four existing points
    SVMPoint *n0, *n1, *n2, *n3;
    poly->getFourPoints(&n0, &n1, &n2, &n3);

    Vec3d p0(n0->u, n0->v, n0->w);
    Vec3d p1(n1->u, n1->v, n1->w);
    Vec3d p2(n2->u, n2->v, n2->w);
    Vec3d p3(n3->u, n3->v, n3->w);

    Vec3d pMouse(imgX, imgY, 1.0);

    /******** BEGIN TODO ********/
    // Find the 2D image positions of box corners p4, p5, p6, p7, as described on the webpage.  
    // You will compute these positions using the known corners of the box (p0, p1, p2, p3, defined above)
    // and the vanishing points.
    // The line through points p4 and p5 will go through the mouse position, pMouse
    // Store the results in variables p4, p5, p6, and p7.
	Vec3d p4, p5, p6, p7;

    //TODO-BLOCK-BEGIN
    Vec3d xV = Vec3d(xVanish.u, xVanish.v, 1.0);
    Vec3d yV = Vec3d(yVanish.u, yVanish.v, 1.0);
    Vec3d zV = Vec3d(zVanish.u, zVanish.v, 1.0);

    Vec3d lm = cross(pMouse,xV);
    lm[0] /= lm[2];
    lm[1] /= lm[2];
    lm[2] = 0.0;

    Vec3d l1 = cross(p0, yV);
    l1[0] /= l1[2];
    l1[1] /= l1[2];
    l1[2] = 0.0;

    p4 = cross(l1,lm);
    p4[0] /= p4[2];
    p4[1] /= p4[2];
    p4[2] = 0.0;

    Vec3d l2 = cross(p1,yV);
    l2[0] /= l2[2];
    l2[1] /= l2[2];
    l2[2] = 0.0;

    p5 = cross(l2,lm);
    p5[0] /= p5[2];
    p5[1] /= p5[2];
    p5[0] = 0.0;

    Vec3d l3 = cross(p3,yV);
    l3[0] /= l3[2];
    l3[1] /= l3[2];
    l3[2] = 0.0;

    Vec3d l4 = cross(p4,zV);
    l4[0] /= l4[2];
    l4[1] /= l4[2];
    l4[2] = 0.0;

    p7 = cross(l3,l4);
    p7[0] /= p7[2];
    p7[1] /= p7[2];
    p7[2] = 0.0;

    Vec3d l5 = cross(p2,yV);
    l5[0] /= l5[2];
    l5[1] /= l5[2];
    l5[2] = 0.0;

    Vec3d l6 = cross(p5,zV);
    l6[0] /= l6[2];
    l6[1] /= l6[2];
    l6[2] = 0.0;

    p6 = cross(l5,l6);
    p6[0] /= p6[2];
    p6[1] /= p6[2];
    p6[2] = 0.0;
    //TODO-BLOCK-END

    /******** END TODO ********/

    p4_out = p4;
    p5_out = p5;
    p6_out = p6;
    p7_out = p7;
}

//
// TODO 8: find3DPositionsBox()
//    Find the 3D positions of the 8 corners of the box.  The 3D position of points[0] is known.
void ImgView::find3DPositionsBox(SVMPoint *points[8]) 
{
    /******** BEGIN TODO ********/
    // Implement this function.  You will compute the 3D positions of the corners of the box
    // using their 2D positions and the known 3D position of points[0].  Store the results in 
    // points[1] through points[7].  You can and should use the sameXY and sameZ routines that 
	// you need to implement.  For that to work, you will need to push and pop points from
	// pntSelStack.  There are multiple ways to implement this function.

    //TODO-BLOCK-BEGIN
    
    // point 1
    pntSelStack.push_back(points[0]);
    pntSelStack.push_back(points[1]);
    sameZPlane();
    pntSelStack.pop_back();
    pntSelStack.pop_back();
    
    // point 2
    pntSelStack.push_back(points[1]);
    pntSelStack.push_back(points[2]);
    sameXY();
    pntSelStack.pop_back();
    pntSelStack.pop_back();
    
    // point 3
    pntSelStack.push_back(points[0]);
    pntSelStack.push_back(points[3]);
    sameXY();
    pntSelStack.pop_back();
    pntSelStack.pop_back();
    
    // point 4
    pntSelStack.push_back(points[0]);
    pntSelStack.push_back(points[4]);
    sameZPlane();
    pntSelStack.pop_back();
    pntSelStack.pop_back();
    
    // point 5
    pntSelStack.push_back(points[4]);
    pntSelStack.push_back(points[5]);
    sameZPlane();
    pntSelStack.pop_back();
    pntSelStack.pop_back();
    
    // point 6
    pntSelStack.push_back(points[5]);
    pntSelStack.push_back(points[6]);
    sameXY();
    pntSelStack.pop_back();
    pntSelStack.pop_back();
    
    // point 7
    pntSelStack.push_back(points[4]);
    pntSelStack.push_back(points[7]);
    sameXY();
    pntSelStack.pop_back();
    pntSelStack.pop_back();

    //TODO-BLOCK-END

	/********* END TODO ********/
}


