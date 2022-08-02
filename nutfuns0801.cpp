
#include "stdafx.h"

#include "ImageProc.h"

using namespace ImageProc;



void Cutter::find_key_points_with_offset(const std::vector<Point2i> nutcontour, Point2f* keypoints, const int offset)
{
	// the most left and right points are found for further process, then
	// the two side points around them are obtained with x- offset by
	// forward and backward searching.

	for (int ix = 1; ix < 20; ix++)
	{
		RotatedRect rRec = minAreaRect(Mat(nutcontour)); // test 
	}

	int idx[2] = { 0,0 };
	for (int id = 0; id < (int)nutcontour.size(); id++)
	{
		if (nutcontour[idx[0]].x > nutcontour[id].x)
		{ // find left point
			idx[0] = id;
		}
		if (nutcontour[idx[1]].x < nutcontour[id].x)
		{ // find right point
			idx[1] = id;
		}
	}
	
	int len = (int)nutcontour.size(), id0, dx = 0;
	int sidepts[2][2] = {};

	for (int jd = 0; jd < 2; jd++)
	{
		for (int id = 0; id < (len / 2); id++)
		{
			// find two side points around the most left/right point
			id0 = (idx[jd] - id + len) % len;

			dx = nutcontour[id0].x - nutcontour[idx[jd]].x;
			dx = (jd == 0 ? dx : -dx);

			if ((dx <= offset) && (dx + 1 > offset))
			{
				sidepts[jd][0] = id0;
			}

			id0 = (idx[jd] + id) % len;

			dx = nutcontour[id0].x - nutcontour[idx[jd]].x;
			dx = (jd == 0 ? dx : -dx);

			if ((dx <= offset) && (dx + 1 > offset))
			{
				sidepts[jd][1] = id0;
			}
		}
	}
	
	keypoints[0].x=  ((nutcontour[sidepts[0][0]].x + nutcontour[sidepts[0][1]].x) / 2.0f);
	keypoints[0].y = ((nutcontour[sidepts[0][0]].y + nutcontour[sidepts[0][1]].y) / 2.0f);
	keypoints[1].x = ((nutcontour[sidepts[1][0]].x + nutcontour[sidepts[1][1]].x) / 2.0f);
	keypoints[1].y = ((nutcontour[sidepts[1][0]].y + nutcontour[sidepts[1][1]].y) / 2.0f);
	//printf("shoulders = (%f, %f) and (%f, %f)\n", keypoints[0].x, keypoints[0].y, keypoints[1].x, keypoints[1].y);

}

RotatedRect Cutter::FindRect(const std::vector<Point2i> nutcontour, const Point2f* keypoints, const int dwidth)
{
	// the main idea:
	// (1) get the coordinates of the center point of the left/right two offseted points;
	// (2) with these two center points, get the unit vector along them and its vertical vector
	// (3) other vectors, from the left center point, to arbitray point on the nut contour, 
	//     can be obtained;
	// (4) with the projection of those vectors onto the first unit vector, we can get the 
	//     maximum and minimum projection along the two center points, then the sum of the
	//     absolute projection value, namely the width of the rectangle, can be obtained;
	// (5) with similiar process, the height of the rectangle can be derived;
	// (6) with the first two projections, the center point along width direction can be obtained,
	//     with the second vertical projections, the center point can shift along vertical vector,
	//     therefore, the center point of the rectangle can be derived;
	// (7) the rotation angle is the angle of the first unit vector;
	// dwidth: unit: pixel.

	// get the unit vector along the two center points
	Vec2f vech(keypoints[1].x - keypoints[0].x, keypoints[1].y - keypoints[0].y);
	float vechn = (float)norm(vech);
	vech[0] = vech[0] / vechn;
	vech[1] = vech[1] / vechn;
	// get the vertical vector
	Vec2f vecv(vech[1], -vech[0]), vecm(0.0f, 0.0f);
	//printf("vech(%f, %f), vecv(%f, %f)\n", vech[0], vech[1], vecv[0], vecv[1]);

	float disth[] = { 0.0f,0.0f }, distv[] = { 0.0f,0.0f }, projh, projv;

	for (int id = 0; id < (int)nutcontour.size(); id++)
	{
		// the vector from left center point to arbitray point of nut contour
		vecm[0] = (float)(nutcontour[id].x - keypoints[0].x);
		vecm[1] = (float)(nutcontour[id].y - keypoints[0].y);

		// the projection along the first unit vector
		projh = vecm.dot(vech);
		if (disth[0] > projh)
		{
			disth[0] = projh;
		}
		if (disth[1] < projh)
		{
			disth[1] = projh;
		}

		// the projection along the vetical vector
		projv = vecm.dot(vecv);
		if (distv[0] > projv)
		{
			distv[0] = projv;
		}
		if (distv[1] < projv)
		{
			distv[1] = projv;
		}
	}

	float width = disth[1] - disth[0] - 2 * dwidth, height = distv[1] - distv[0]; // width, height of the rectangle
	
	Vec2f increment = (disth[0] + disth[1]) / 2 * vech + (distv[0] + distv[1]) / 2 * vecv; // shift vector of the center point
	
	Point2f pcenter;
	pcenter.x = keypoints[0].x + increment[0];
	pcenter.y = keypoints[0].y + increment[1];

	float pi = 2 * cos(0.0f);
	float angle = (atan2(vech[1], vech[0]) * 180 / pi);
	return RotatedRect(pcenter, Size2f(width, height), 90 - angle);
}


void Cutter::rotate_contour_by_angle(std::vector<Point2f>& nutcontour, RotatedRect rRect, float angle)
{
	Point2f pcenter = rRect.center;
	float dx, dy;

	float pi = 2 * acos(0.0f);
	float cosx = cos(angle * pi / 180), sinx = sin(angle * pi / 180);

	for (int id = 0; id < (int)nutcontour.size(); id++)
	{
		dx = nutcontour[id].x - pcenter.x;
		dy = nutcontour[id].y - pcenter.y;

		nutcontour[id].x = dx * cosx + dy * sinx + pcenter.x;
		nutcontour[id].y = dy * cosx - dx * sinx + pcenter.y;
	}
}

	
void Cutter::trim_contour_with_offset(std::vector<Point2i>& nutcontour, const int dwidth)
{
	// a new contour must be copied before using this function, 
	// otherwise the original contour coordinates might be changed since
	// conversion from float to int and vice versa. though the coordinates 
	// are changed(slightly), the point order has not been changed.
	//
	Point2f keypoints[2];
	find_key_points_with_offset(nutcontour, keypoints, 30);
	RotatedRect rRect = FindRect(nutcontour, keypoints, 0);

	std::vector<Point2f> contourx = convert_contour2i_to_contour2f(nutcontour);
	rotate_contour_by_angle(contourx, rRect, -rRect.angle);

	float half_width = rRect.size.width / 2 - (float)dwidth;

	for (int id = 0; id < (int)contourx.size(); id++)
	{
		if (contourx[id].y < rRect.center.y - half_width)
		{
			contourx[id].y = rRect.center.y - half_width;
		}
		if (contourx[id].y > rRect.center.y + half_width)
		{
			contourx[id].y = rRect.center.y + half_width;
		}
	}
	rotate_contour_by_angle(contourx, rRect, +rRect.angle);
	nutcontour = convert_contour2f_to_contour2i(contourx);
}


std::vector<Point2i> Cutter::convert_contour2f_to_contour2i(std::vector<Point2f> contour2f)
{
	std::vector<Point2i> contour2i(contour2f.size());

	for (int id = 0; id < (int)contour2i.size(); id++)
	{
		contour2i[id].x = (int)(contour2f[id].x);
		contour2i[id].y = (int)(contour2f[id].y);
	}
	return contour2i;
}


std::vector<Point2f> Cutter::convert_contour2i_to_contour2f(std::vector<Point2i> contour2i)
{
	std::vector<Point2f> contour2f(contour2i.size());

	for (int id = 0; id < (int)contour2f.size(); id++)
	{
		contour2f[id].x = (float)contour2i[id].x;
		contour2f[id].y = (float)contour2i[id].y;
	}
	return contour2f;
}


bool Cutter::remove_branches(Mat1b& fstMat, Mat1b& ctrMat, std::vector<Point2i>& contour, int bThresh)
{
	// Input: srx, grayscale image;
	// generate a black image with a white nut in, and return this image and contour.
	bThresh = 50;
	Mat1b tmpMat = fstMat.clone(); // nut without branches.
	Mat se = getStructuringElement(MORPH_RECT, cv::Size(bThresh, bThresh));
	erode(tmpMat, tmpMat, se);
	dilate(tmpMat, tmpMat, se);

	std::vector <std::vector<Point2i>> contours;
	findContours(tmpMat, contours, CV_RETR_TREE, CV_CHAIN_APPROX_NONE);

	int idx = -1;
	for (int id = 0; id < (int)contours.size(); id++)
	{
		int xsize = (int)contours[id].size();
		if (xsize > 500 && xsize < 1500)
		{
			idx = id;
		}
	}

	if (idx < 0) { return false; }

	contour = contours[idx]; // the size of nut contour should within [500, 1500] pixels.
	contour[1].y=contour[1].y+40;
	contour[50].y=contour[50].y+80;
	contour[100].y=contour[100].y+120;
	ctrMat = Mat1b::zeros(fstMat.size()); // only nut, without small pots and branches.
	drawContours(ctrMat, contours, idx, Scalar(1), CV_FILLED);
	return true;
}


void Cutter::FindFoot(std::vector<Point2i> &contour, Point2i &headpoint, Point2i &footpoint)
{
	// trim contours with width offset, 6.16
	//trim_contour_with_offset(contour, 30);
	// trim end
	float distmax = 0, dist = 0;
	for (int id=0; id < (int)contour.size();id++)
	{
		dist = (float)norm(headpoint-contour[id]);
		if (dist >= distmax)
		{
			distmax = dist;
			footpoint = contour[id];
		}
	}
}


// overload function of FindHead, used for detect2 and detect3.
bool Cutter::FindHead(Mat1b grayMat, std::vector<Point2i>& contour, Point2i &head, Point2i &foot)
{
	// trim contours with width offset, 6.16
	///trim_contour_with_offset(contour, 20);
	// trim end

	// obtained trimed image with white nut in.
	Mat1b rData = Mat1b::zeros(grayMat.size());
	std::vector<std::vector<Point2i>> contours;
	contours.push_back(contour);
	drawContours(rData, contours, 0, Scalar(1), CV_FILLED);

	cv::Rect rect = boundingRect(Mat(contour));

	// find two farest points in contour
	float distmax = 0, dist = 0;
	Point2i pt1, pt2;

	for (int id = 0; id < (int)contour.size() - 1; id++)
	{
		for (int jd =id + 1;jd < (int)contour.size();jd++) 
		{
			dist = (float)norm(contour[id]-contour[jd]);
			if (dist >= distmax)
			{
				distmax = dist;
				pt1 = contour[id];
				pt2 = contour[jd];
			}
		}
	}

	// find head
	bool bhead = FindHead(rData, contour, rect, pt1, pt2, head);
	if (!bhead) {return false;}

	// find foot
	distmax = 0, dist = 0;
	for (int id = 0; id < (int)contour.size(); id++)
	{
		dist = (float)norm(head-contour[id]);
		if (dist >= distmax)
		{
			distmax = dist;
			foot = contour[id];
		}
	}
	return true;
}


