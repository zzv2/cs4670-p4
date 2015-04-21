/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #2:
 * ImgView.inl (included from ImgView.cpp)
 *		contains routines for computing the 3D position of points
 ***************************************************************/

//
// TODO 4: sameXY()
//		Computes the 3D position of newPoint using knownPoint
//		that has the same X and Y coordinate, i.e. is directly
//		below or above newPoint.
//		See lecture slide on measuring heights.
//
// HINT1: make sure to dehomogenize points when necessary
// HINT2: there is a degeneracy that you should look out for involving points already in line with the reference
// HINT3: make sure to get the sign of the result right, i.e. whether it is above or below ground
// HINT4: check for some helpful/necessary variables which are listed in ImgView.h such as the vanishing points and homography H
void ImgView::sameXY()
{
	if (pntSelStack.size() < 2)
	{
		fl_alert("Not enough points on the stack.");
		return;
	}

	SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
	SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];

	if( !knownPoint.known() )
	{
		fl_alert("Can't compute relative values for unknown point.");
		return;
	}

	if( refPointOffPlane == NULL )
	{
		fl_alert("Need to specify the reference height first.");
		return;
	}

	/******** BEGIN TODO ********/

	// See the lecture note on measuring heights
	// using a known point directly below the new point.

	// printf("sameXY() to be implemented!\n");


    //TODO-BLOCK-BEGIN
	// todo - sign?
	// todo - degenerative case?

	Vec3d vx = Vec3d(xVanish.u, xVanish.v, 1.0);
	Vec3d vy = Vec3d(yVanish.u, yVanish.v, 1.0);
	Vec3d vz = Vec3d(zVanish.u, zVanish.v, 1.0);

    double bu, bv;
    ApplyHomography(bu, bv, H, refPointOffPlane->X, refPointOffPlane->Y, 1.0);
    Vec3d b = Vec3d(bu, bv, 1.0);
    Vec3d r = Vec3d(refPointOffPlane->u, refPointOffPlane->v, 1.0);

	Vec3d b0 = Vec3d(knownPoint.u, knownPoint.v, 1.0);
	Vec3d t0 = Vec3d(newPoint.u, newPoint.v, 1.0);

    Vec3d v = cross(cross(b, b0), cross(vx, vy));
    v[0] = v[0]/v[2];
    v[1] = v[1]/v[2];
    v[2] = 1.0;

    Vec3d t = cross(cross(v, t0), cross(r, b));
    t[0] = t[0]/t[2];
    t[1] = t[1]/t[2];
    t[2] = 1.0;

    Vec3d t_b = t - b;
    Vec3d vz_r = vz - r;
    Vec3d r_b = r - b;
    Vec3d vz_t = vz - t;
    double cross_ratio = t_b.length() * vz_r.length()/(r_b.length() * vz_t.length());

    double height = cross_ratio * referenceHeight;
    newPoint.X = knownPoint.X;
    newPoint.Y = knownPoint.Y;
    newPoint.Z = knownPoint.Z + height;
    //TODO-BLOCK-END
	/******** END TODO ********/

	newPoint.known(true);

	printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );

	redraw();
}



//
// TODO 5: sameZPlane()
//		Compute the 3D position of newPoint using knownPoint
//		that lies on the same plane and whose 3D position is known.
//		See the man on the box lecture slide.
//		If newPoint is on the reference plane (Z==0), use homography (this->H, or simply H) directly.
//
// HINT: For this function, you will only need to use the three vanishing points and the reference homography 
//       (in addition to the known 3D location of knownPoint, and the 2D location of newPoint)
void ImgView::sameZPlane()
{
	if (pntSelStack.size() < 2)
	{
		fl_alert("Not enough points on the stack.");
		return;
	}

	SVMPoint &newPoint = *pntSelStack[pntSelStack.size() - 1];
	SVMPoint &knownPoint = *pntSelStack[pntSelStack.size() - 2];

	if( !knownPoint.known() )
	{
		fl_alert("Can't compute relative values for unknown point.");
		return;
	}

	/******** BEGIN TODO ********/
    //TODO-BLOCK-BEGIN
    if (knownPoint.Z == 0)
    {
    	double X, Y;
    	ApplyHomography(X, Y, Hinv, newPoint.u, newPoint.v, 1.0);
    	newPoint.X = X;
    	newPoint.Y = Y;
    	newPoint.Z = 0.0;
    }
     else
    {
    	// Get the vanishing points
    	Vec3d vx = Vec3d(xVanish.u, xVanish.v, 1.0);
    	Vec3d vy = Vec3d(yVanish.u, yVanish.v, 1.0);
    	Vec3d vz = Vec3d(zVanish.u, zVanish.v, 1.0);

    	// Get t1 in image coords
    	Vec3d t1 = Vec3d(knownPoint.u, knownPoint.v, 1.0);

    	// Get m0 in image coords
    	Vec3d m0 = Vec3d(newPoint.u, newPoint.v, 1.0);

    	// Find b1 from t1
    	double _u, _v;
    	ApplyHomography(_u, _v, H, knownPoint.X, knownPoint.Y, 1.0);
    	Vec3d b1 = Vec3d(_u, _v, 1.0);

    	// Find v
    	Vec3d v = cross(cross(t1, m0), cross(vx, vy));
    	v[0] = v[0]/v[2];
    	v[1] = v[0]/v[2];
    	v[2] = 1.0;

    	// Find b0
    	Vec3d b0 = cross(cross(m0, vz), cross(b1, v));
    	b0[0] = b0[0]/b0[2];
    	b0[1] = b0[1]/b0[2];
    	b0[2] = 1.0;

    	// Finally find the new point's 3d coords
    	double x, y;
    	ApplyHomography(x, y, Hinv, b0[0], b0[1], b0[2]);
    	newPoint.X = x;
    	newPoint.Y = y;
    	newPoint.Z = knownPoint.Z;

    }
    //TODO-BLOCK-END
	/******** END TODO ********/

	newPoint.known(true);

	printf( "Calculated new coordinates for point: (%e, %e, %e)\n", newPoint.X, newPoint.Y, newPoint.Z );

	redraw();
}


