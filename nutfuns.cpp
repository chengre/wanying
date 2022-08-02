
#include "stdafx.h"

#include "ImageProc.h"

using namespace ImageProc;



void Cutter::find_the_left_right_points(const std::vector<Point2i>& cnt, int idx[4])
{
	for (int id = 0; id < (int)cnt.size(); id++)
	{
		if (cnt[idx[0]].x > cnt[id].x)
		{ // find leftmost point
			idx[0] = id;
		}
		if (cnt[idx[1]].x < cnt[id].x)
		{ // find rightmost point
			idx[1] = id;
		}
		if (cnt[idx[2]].y > cnt[id].y)
		{ // the topmost point
			idx[2] = id;
		}
		if (cnt[idx[3]].y < cnt[id].y)
		{ // the downmost point
			idx[3] = id;
		}
	}
}



void Cutter::find_the_left_right_centers(const std::vector<Point2i>& cnt, Point2i& ct1, Point2i& ct2)
{
	int dw = 10;
	find_the_left_right_centers(cnt, dw, dw, ct1, ct2);
}



void Cutter::find_the_left_right_centers(const std::vector<Point2i>& cnt, int indentleft, int indentright, Point2i& ct1, Point2i& ct2)
{
	int idx[4] = {};
	find_the_left_right_points(cnt, idx);

	Point2i pt = cnt[idx[1]] - cnt[idx[0]];

	int dx = pt.x, dy = pt.y;
	double dl = sqrt(double(dx * dx + dy * dy));

	double vech[2] = { dx / dl,dy / dl }; // from left to right

	int indents[2] = { indentleft, indentright };
	int pid[4] = {};
	int cntlen = (int)cnt.size();

	for (int kd = 0; kd < 2; kd++) // left and right side
	{
		for (int rd = 0; rd < 2; rd++) // top and down side
		{
			for (int id = 0; id < (cntlen / 6); id++)
			{
				// sgn = -1: backward; sgn = 1: forward.
				int sgn = (rd < 1 ? 1 : -1);

				int ip = (idx[kd] + sgn * id + cntlen) % cntlen;
				Point2i pti = cnt[ip] - cnt[idx[kd]];

				double dist = pti.x * vech[0] + pti.y * vech[1]; // projection
				if (kd > 0) { dist = -dist; }

				if (dist >= indents[kd] && dist < indents[kd] + 1)
				{
					pid[2 * kd + rd] = ip;
					break;
				}
			}
		}
	}

	//ct1 = (cnt[pid[0]] + cnt[pid[1]]) / 2; // left center
	//ct2 = (cnt[pid[2]] + cnt[pid[3]]) / 2; // right center
	ct1.x = (cnt[pid[0]].x + cnt[pid[1]].x) / 2;
	ct1.y = (cnt[pid[0]].y + cnt[pid[1]].y) / 2;
	ct2.x = (cnt[pid[2]].x + cnt[pid[3]].x) / 2;
	ct2.y = (cnt[pid[2]].y + cnt[pid[3]].y) / 2;
}



RotatedRect Cutter::find_nut_rectangle(const std::vector<Point2i>& cnt, int dwidth)
{
	// Input: dwidth: unit: pixel, offset from the left and right points to center.
	//

	if (dwidth < 5) { dwidth = 5; }
	if (dwidth > 90) { dwidth = 90; }

	Point2i ct1, ct2, ct, cpt;
	find_the_left_right_centers(cnt, dwidth, dwidth, ct1, ct2);
	ct = ct2 - ct1, cpt = (Point2i)((ct1 + ct2) * 0.5);

	double dl = norm(ct);
	double vech[2] = { ct.x / dl, ct.y / dl }, vecv[2] = { vech[1], -vech[0] };

	double disth[] = { 0.0,0.0 }, distv[] = { 0.0,0.0 };

	for (int id = 0; id < (int)cnt.size(); id++)
	{
		Point2i pti = cnt[id] - cpt;

		double projh = pti.x * vech[0] + pti.y * vech[1]; // projection on vech
		if (disth[0] > projh)
		{
			disth[0] = projh;
		}
		if (disth[1] < projh)
		{
			disth[1] = projh;
		}

		double projv = pti.x * vecv[0] + pti.y * vecv[1]; // projection on vecv
		if (distv[0] > projv)
		{
			distv[0] = projv;
		}
		if (distv[1] < projv)
		{
			distv[1] = projv;
		}
	}

	double width = disth[1] - disth[0], height = distv[1] - distv[0]; // width, height of the rectangle

	Point2d pcenter;
	double disthmean = (disth[0] + disth[1]) / 2, distvmean = (distv[0] + distv[1]) / 2;

	pcenter.x = cpt.x + disthmean * vech[0] + distvmean * vecv[0]; // shift vector of the center point
	pcenter.y = cpt.y + disthmean * vech[1] + distvmean * vecv[1]; // shift vector of the center point

	double pi = 2 * cos(0.0);
	double angle = (atan2(vech[1], vech[0]) * 180 / pi);
	return RotatedRect((Point2f)pcenter, Size2f((float)width, (float)height), 90 - (float)angle);
}



void Cutter::rotate_contour_by_angle(std::vector<Point2f>& cnt, RotatedRect rRect, float angle)
{
	Point2f pcenter = rRect.center;
	float dx, dy;

	float pi = 2 * acos(0.0f);
	float cosx = cos(angle * pi / 180), sinx = sin(angle * pi / 180);

	for (int id = 0; id < (int)cnt.size(); id++)
	{
		dx = cnt[id].x - pcenter.x;
		dy = cnt[id].y - pcenter.y;

		cnt[id].x = dx * cosx + dy * sinx + pcenter.x;
		cnt[id].y = dy * cosx - dx * sinx + pcenter.y;
	}
}



void Cutter::get_nut_bending_tendency(std::vector<Point2i>& cnt, int headside, double& indicator)
{
	// Purpose: to evaluate the similarity between the nut shell and the banana shape.
	// Output: indicator, the more it close to 1, the less similar to banana shape;
	//		   the more it close to 0, the more similar to banana shape.
	//
	Point2i ct1, ct2;
	find_the_left_right_centers(cnt, ct1, ct2);
	int dx = ct2.x - ct1.x, dy = ct2.y - ct1.y;

	Point2i cpt = (Point2i)((ct1 + ct2) * 0.5);
	Point2i cpv = cpt + Point2i(dy, -dx);

	int idx[2] = { -1,-1 };
	double xlength = 0;
	get_length_of_splited_contour(cnt, cpt, cpv, headside, 1000, 0, idx, xlength);

	double hleft = norm(cnt[idx[1]] - cpt);
	double hright = norm(cnt[idx[0]] - cpt);

	indicator = (hleft < hright ? hleft / hright : hright / hleft);
}



void Cutter::trim_contour_with_offset(std::vector<Point2i>& cnt, int dwidth)
{
	// a new contour must be copied before using this function, 
	// otherwise the original contour coordinates might be changed since
	// conversion from float to int and vice versa. though the coordinates 
	// are changed(slightly), the point order has not been changed.
	//

	Point2i ct1, ct2;
	find_the_left_right_centers(cnt, 30, 30, ct1, ct2);
	Point2i ct = ct2 - ct1;
	double dl = norm(ct);
	double vech[2] = { ct.x / dl,ct.y / dl };

	RotatedRect rRect = find_nut_rectangle(cnt, dwidth);

	std::vector<Point2f> contourx = convert_contour2i_to_contour2f(cnt);
	rotate_contour_by_angle(contourx, rRect, rRect.angle);

	float half_width = rRect.size.width / 2 - (float)dwidth;

	for (int id = 0; id < (int)contourx.size(); id++)
	{
		if (contourx[id].x < rRect.center.x - half_width)
		{
			contourx[id].x = rRect.center.x - half_width;
		}
		if (contourx[id].x > rRect.center.x + half_width)
		{
			contourx[id].x = rRect.center.x + half_width;
		}
	}
	rotate_contour_by_angle(contourx, rRect, -rRect.angle);
	cnt = convert_contour2f_to_contour2i(contourx);
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



void Cutter::sort_contours_by_area(std::vector <std::vector<Point2i>>& cnts)
{
	// sort contours, with the largest area first, and smallest the last.
	for (int id = 0; id < (int)cnts.size() - 1; id++)
	{
		int kd = id;
		for (int jd = id + 1; jd < (int)cnts.size(); jd++)
		{
			if (contourArea(Mat(cnts[kd])) < contourArea(Mat(cnts[jd])))
			{
				kd = jd;
			}
		}
		if (kd > id)
		{
			std::vector<Point2i> cnti = cnts[id];
			cnts[id] = cnts[kd];
			cnts[kd] = cnti;
		}
		if (id > 5) { break; } // only find the top 5 maximum area.
	}
}



void Cutter::sort_contours_by_length(std::vector <std::vector<Point2i>>& cnts)
{
	// sort contours, with the longest first, and the shortest the last.
	// usually used for sorting possible branches, the longest one maybe a branch
	//

	std::vector<double> maxlens;

	for (int id = 0; id < (int)cnts.size(); id++)
	{
		int idx[4] = {};
		find_the_left_right_points(cnts[id], idx);
		double maxlen = max(norm(cnts[id][idx[1]] - cnts[id][idx[0]]), norm(cnts[id][idx[3]] - cnts[id][idx[2]]));

		maxlens.push_back(maxlen);
	}

	for (int id = 0; id < (int)cnts.size() - 1; id++)
	{
		int kd = id;
		for (int jd = id + 1; jd < (int)cnts.size(); jd++)
		{
			if (maxlens[kd] < maxlens[jd]) { kd = jd; }
		}
		if (kd > id)
		{
			std::vector<Point2i> cnti = cnts[id];
			cnts[id] = cnts[kd];
			cnts[kd] = cnti;

			double maxlen = maxlens[id];
			maxlens[id] = maxlens[kd];
			maxlens[kd] = maxlen;
		}
		if (id > 4) { break; } // only find the top 5 longest contours.
	}
}



double Cutter::distance_from_point_to_line(Point2i pt, Point2i pt1, Point2i pt2)
{
	// get the distance from point pt to the line through pt1 and pt2, it should 
	// be noted that the distance can be positive(negative) if pt is at the left(right) 
	// half plane of vector {pt2 - pt1}.
	//
	int dx = pt2.x - pt1.x, dy = pt2.y - pt1.y;
	double dl = sqrt((double)(dx * dx + dy * dy));
	if (dl < 1.0e-10) { dl = 1.0e-10; } // dl should not be zero.
	// distance has direction. the vertical vector {dy, -dx} is 90 deg ahead of vector {pt2 - pt1},
	// so if pt is at the left(right) half plane of vector {pt2 - pt1}, distance is positive(negative).
	return (dy * (pt.x - pt1.x) - dx * (pt.y - pt1.y)) / dl; // see formula of distance from point to line.
}



bool Cutter::get_length_of_splited_contour(std::vector<Point2i>& cnt, Point2i pti, Point2i ptj, int headside, int floweroffset, int flowerindent, int idx[2], double& xlength)
{
	// Input: pti, ptj are the two points to generate base point and line direction.
	// the base point is pti, which must be an inner point, and the direction is (ptj - pti).
	// floweroffset: always positive, no matter at left or right of the nut.
	// flowerindent: the length of flower to be reduced from the line segment length.
	// xlength: the length of the line segment, after recognizing flower.
	// Output: idx[2] are the two intersection points on cnt. the vector {cnt[idx[1]] - cnt[idx[0]]}
	// is of the same direction as the vector {ptj - pti}.
	//

	int dx = ptj.x - pti.x, dy = ptj.y - pti.y;
	idx[0] = -1, idx[1] = -1; // initialization of the intersection points

	for (int id = 0; id < (int)cnt.size(); id++)
	{
		double dist = distance_from_point_to_line(cnt[id], pti, ptj);

		if (std::abs(dist) < 1) // cnt[id] is on the line through pti and ptj
		{
			int dxt = cnt[id].x - pti.x, dyt = cnt[id].y - pti.y;

			(dxt * dx + dyt * dy < 0 ? idx[0] = id : idx[1] = id);
		}
	}
	// if less than two points were found, return false.
	// return (idx[0] > 0 && idx[1] > 0) ? true : false;
	if (idx[0] < 0 || idx[1] < 0) { return false; }

	// to see if any point in idx is at the outside of the floweroffset, if true, change the line length
	xlength = norm(cnt[idx[0]] - cnt[idx[1]]); // length of the line segment

	int idt[4] = {};
	find_the_left_right_points(cnt, idt);

	Point2i cpt = (Point2i)((cnt[idt[0]] + cnt[idt[1]]) * 0.5);// now we get the coordinates of the center point
	Point2i hpt = (headside < 0 ? cnt[idt[0]] : cnt[idt[1]]);

	dx = hpt.x - cpt.x, dy = hpt.y - cpt.y; // vech[2] ={dx, dy}
	double dl = norm(Point2i(dx, dy));
	double vech[] = { dx / dl,dy / dl }; // the unit vector

	for (int id = 0; id < 2; id++)
	{
		int dxt = cnt[idx[id]].x - cpt.x, dyt = cnt[idx[id]].y - cpt.y; //

		if ((dxt * vech[0] + dyt * vech[1]) < -floweroffset) // projection of [cnt[id] - cpt] on [hpt - cpt]
		{
			xlength -= flowerindent; // change length of line segment
		}
	}
	return true;
}



bool Cutter::get_heights_of_splited_contour(std::vector<Point2i>& cnt, Point2i pt1, Point2i pt2, int headside, int floweroffset, double xheights[2])
{
	// Input: pt1, pt2 defines a vector {pt2 - pt1}, and the line to seperate cnt.
	// floweroffset: always positive, no matter at left or right of the nut.
	// Output: xheights[0], xheight[1], both positive, are the left, right, of above vector, maximum height of the half contour.
	//
	
	int idt[4] = {};
	find_the_left_right_points(cnt, idt);
	
	Point2i cpt = (Point2i)((cnt[idt[0]] + cnt[idt[1]]) * 0.5);// now we get the coordinates of the center point
	Point2i hpt = (headside < 0 ? cnt[idt[0]] : cnt[idt[1]]);

	int dx = hpt.x - cpt.x, dy = hpt.y - cpt.y; // vech[2] ={dx, dy}
	double dl = norm(Point2i(dx, dy));
	double vech[] = { dx / dl,dy / dl }; // the unit vector

	xheights[0] = 0, xheights[1] = 0; // initialization of heights.

	for (int id = 0; id < (int)cnt.size(); id++)
	{
		int dxt = cnt[id].x - cpt.x, dyt = cnt[id].y - cpt.y; //
		// ignore those points at the outer side of flower
		if ((dxt * vech[0] + dyt * vech[1]) < -floweroffset) { continue; } // projection of [cnt[id] - cpt] on [hpt - cpt]

		double dist = distance_from_point_to_line(cnt[id], pt1, pt2); // left half plane: positive; right half plane: negative.

		if (dist > xheights[0]) { xheights[0] = dist; }
		if (dist < xheights[1]) { xheights[1] = dist; }
	}

	xheights[1] = -xheights[1]; // to ensure heights are positive values.
	return (xheights[0] == 0 || xheights[1] == 0 ? false : true);
}



bool Cutter::remove_branches(Mat1b& fstMat, Mat1b& ctrMat, std::vector<Point2i>& contour, int branchwidththresh)
{
	// Input: srx, grayscale image;
	// generate a black image with a white nut in, and return this image and contour.
	branchwidththresh = 50;
	Mat1b tmpMat = fstMat.clone(); // nut without branches.
	Mat se = getStructuringElement(MORPH_RECT, cv::Size(branchwidththresh, branchwidththresh));
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
	//contour[1].y=contour[1].y+40;
	//contour[50].y=contour[50].y+80;
	//contour[100].y=contour[100].y+120;
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



// PART II: functions with image as input.

bool Cutter::receive_image_find_nut_contour(Mat1b& src, int nutside, int graythresh, std::vector<Point2i>& cnt, Mat1b& srx)
{
	// nutside = -1: find nut in left region; nutside = 1: right region.

	convert_colored_to_binary_image(src, srx, graythresh);

	// erase the margin
	int margin = 10; // margin are to be erased.

	// find nut contours in left and right regions
	std::vector<std::vector<Point2i>> cnts;
	if (!find_valid_nut_contours(srx, cnts)) { return false; }
	
	int width = srx.cols, height = srx.rows;
	int pos = 0;
	for (int id = 0; id < (int)cnts.size(); id++)
	{
		int idx[4] = {};
		find_the_left_right_points(cnts[id], idx);
		int xc = (cnts[id][idx[0]].x + cnts[id][idx[1]].x) / 2;

		if (nutside > 0) // find nut in right region
		{
			if (xc > width / 2 && (cnts[id][idx[1]].x < width - margin))
			{
				pos = id;
			}
		}
		else // find nut in left region
		{
			if (xc < width / 2 && (cnts[id][idx[0]].x > margin))
			{
				pos = id;
			}
		}
	}
	
	cnt = cnts[pos];
	srx = Mat1b::zeros(srx.size());
	drawContours(srx, cnts, pos, Scalar(255), CV_FILLED);

	return true;
}



void Cutter::convert_colored_to_binary_image(Mat1b& src, Mat1b& srx, int graythresh)
{
	// cvtColor(src, srx, CV_BGR2GRAY);
	// convert the background to black and nut to white.
	// since the command 'findContours' only find white(not black) pixels as contour,
	// and the command 'drawContours' set the contour and its inside as white or black.
	cv::threshold(src, srx, graythresh, 255, THRESH_BINARY+THRESH_BINARY_INV);

	// erase black pots on nut
	Mat sx = getStructuringElement(MORPH_ELLIPSE, cv::Size(5, 5));
	dilate(srx, srx, sx);
	erode(srx, srx, sx);
}



bool Cutter::find_valid_nut_contours(Mat1b& srx, std::vector<std::vector<Point2i>>& cnts)
{
	// binary image is default with black background and white foreground.
	// remove small black pots on nut.
	cnts.clear();
	std::vector <std::vector<Point2i>> cntx;
	findContours(srx, cntx, RETR_LIST, CHAIN_APPROX_NONE); // find white (not black) points as contour

	srx = Mat1b::zeros(srx.size());

	for (int id = 0; id < (int)cntx.size(); id++)
	{
		// no more than three nut contours will be found.
		if (id > 2) { break; }

		// find the maximum nut contour in the rest contours.
		for (int jd = id + 1; jd < (int)cntx.size(); jd++)
		{
			if (cntx[id].size() < cntx[jd].size())
			{
				std::vector<Point2i> cnti = cntx[id];
				cntx[id] = cntx[jd];
				cntx[jd] = cnti;
			}
		}
		
		int idx[4] = {};
		find_the_left_right_points(cntx[id], idx);
		int width = cntx[id][idx[1]].x - cntx[id][idx[0]].x, height = cntx[id][idx[3]].y - cntx[id][idx[2]].y;
		
		if ((width > minLength/dimPerPix) && (width < overMaxLength/dimPerPix) && (height > minHeight/dimPerPix))
		{
			cnts.push_back(cntx[id]);
			drawContours(srx, cntx, id, Scalar(255), CV_FILLED); // set the contour and its interior as white
		}
	}
	
	return cnts.size() < 1 ? false : true;
}



void Cutter::find_head_by_cluster_method(Mat1b& srx, int& headside)
{
	// assuming the nut contour has no branches at tail.
	std::vector<std::vector<Point2i>> cnts;
	if (!find_valid_nut_contours(srx, cnts)) { return; }

	std::vector<Point2i> cnt = cnts[0];
	int idx[4] = {};
	find_the_left_right_points(cnt, idx);

	int offset=10, dwidth = (int)((cnt[idx[1]].x - cnt[idx[0]].x) / 10); // 1/10 is the width scale of interest region.

	int whitecounts[2] = { 0,0 };

	for (int jd = (int)cnt[idx[2]].y; jd < (int)cnt[idx[3]].y; jd++)
	{
		for (int id = (int)cnt[idx[0]].x + offset; id < (int)cnt[idx[0]].x + offset + dwidth; id++)
		{
			if (srx.at<uchar>(jd, id) == 255)
			{
				whitecounts[0]++; // count nut pixels in left side
			}
		}

		for (int id = (int)cnt[idx[1]].x - offset; id > (int)cnt[idx[1]].x - offset - dwidth; id--)
		{
			if (srx.at<uchar>(jd, id) == 255)
			{
				whitecounts[1]++; // count nut pixels in right side
			}
		}
	}

	// -1: head is left; +1: head is right; 0: unknown.
	headside = (whitecounts[0] < whitecounts[1] ? -1 : 1); // -1: head is left; +1: head is right.
}



void Cutter::trim_branch_and_find_head(Mat1b& srx, int branchwidththresh, int branchlengththresh, int& headside)
{
	std::vector<std::vector<Point2i>> cnts;
	if (!find_valid_nut_contours(srx, cnts)) { return; }
	
	Mat1b srz = Mat1b::zeros(srx.size()); // srz storages all possible branches.

	// remove all posible branches, including the true branches, the flower branches and the nut head.
	int width = srx.cols, height = srx.rows;

	for (int id = 0; id < width - 1; id++)
	{
		for (int jd = 0; jd < height - 1; jd++)
		{
			if (srx.at<uchar>(jd, id) == 0 && srx.at<uchar>(jd + 1, id) == 255)
			{
				if (jd + branchwidththresh < height)
				{
					for (int kd = 2; kd < branchwidththresh; kd++)
					{
						if (srx.at<uchar>(jd + kd - 1, id) == 255 && srx.at<uchar>(jd + kd, id) == 0)
						{
							for (int hd = 1; hd < kd; hd++)
							{
								srz.at<uchar>(jd + hd, id) = srx.at<uchar>(jd + hd, id);
							}
						}
					}
				}
			}
		}
	}
	
	std::vector<std::vector<Point2i>> cntz; // cntz are trimed parts from nut, it possibly contains branches.
	findContours(srz, cntz, RETR_LIST, CHAIN_APPROX_NONE);
	//sort_contours_by_area(cntz); // sort the trimed contours, the larger the closer to first, so the branches might be the first.
	sort_contours_by_length(cntz);

	// if all the possible branches are trimed, then the head and tail may also be trimed unexpectedly.
	// headside < 0: head is left; > 0: head is right; == 0: unknown.
	headside = 0;
	
	if ((int)cntz.size() > 0)
	{
		// find the center of the nut
		int idx[4] = {};
		find_the_left_right_points(cnts[0], idx);
		int xc = (cnts[0][idx[0]].x + cnts[0][idx[1]].x) / 2;

		// find the center of the branch
		int idz[4] = {};
		find_the_left_right_points(cntz[0], idz);
		int xz = (cntz[0][idz[0]].x + cntz[0][idz[1]].x) / 2;

		double maxlen = max(norm(cntz[0][idz[1]] - cntz[0][idz[0]]), norm(cntz[0][idz[3]] - cntz[0][idz[2]]));

		// remove a branch if it is longer or higher than the predefined branch length thresh
		if (maxlen > branchlengththresh)
		{
			// remove branch from the orignal image.
			drawContours(srx, cntz, 0, Scalar(0), CV_FILLED);
			// check which side is the head.
			headside = (xz > xc ? -1 : +1);
		}
	}

	// find the nut contours, all possible branches have been removed from the image.
	cnts.clear();
	find_valid_nut_contours(srx, cnts);
}



void Cutter::find_flower_and_find_head(Mat1b& srx, int flowerwidththresh, int flowerheightthresh, double mblankscale, double mwidthscale, int& headside, int& floweroffset)
{
	// flowers will not be trimed, actually, nothing will be trimed.
	// triming branch would be used before this command.
	// output: floweroffset: unit(pixel), from the center of cpt to flower, is always positive, no matter flower is at left or right of nut.
	//		   when no flower exist, flowerdist = 1000 by default.
	//
	std::vector<std::vector<Point2i>> cnts;
	if (!find_valid_nut_contours(srx, cnts)) { return; }

	// the most left and right points will be found.
	std::vector<Point2i> cnt = cnts[0]; // trimed contour.
	int idx[4] = {};
	find_the_left_right_points(cnt, idx);

	// flowers exist in [pt1.x + x1, pt1.x + x2] and [pt2.x - x2, pt2.x - x1] regions.
	int nutwidth = cnt[idx[1]].x - cnt[idx[0]].x, xc = (cnt[idx[0]].x + cnt[idx[1]].x) / 2;
	int mblank = (int)(nutwidth * mblankscale), mwidth = (int)(nutwidth * mwidthscale);
	
	Mat1b srz = Mat1b::zeros(srx.size());
	int width = srx.cols, height = srx.rows;

	for (int jd = (int)cnt[idx[2]].y; jd < (int)cnt[idx[3]].y; jd++)
	{
		for (int id = (int)cnt[idx[0]].x + mblank + mwidth; id > (int)cnt[idx[0]].x + mblank; id--)
		{
			if (id - flowerwidththresh > 0)
			{
				if (srx.at<uchar>(jd, id) == 255 && srx.at<uchar>(jd, id + 1) == 0)
				{
					for (int kd = flowerwidththresh; kd > 0; kd--)
					{
						if (srx.at<uchar>(jd, id - kd) == 0)
						{
							for (int hd = 0; hd < kd; hd++)
							{
								srz.at<uchar>(jd, id - hd) = srx.at<uchar>(jd, id - hd);
							}
						}
						break;
					}
				}
			}
		}

		for (int id = (int)cnt[idx[1]].x - mblank - mwidth; id < (int)cnt[idx[1]].x - mblank; id++)
		{
			if (id + flowerwidththresh < width)
			{
				if (srx.at<uchar>(jd, id) == 255 && srx.at<uchar>(jd, id - 1) == 0)
				{
					for (int kd = flowerwidththresh; kd > 0; kd--)
					{
						if (srx.at<uchar>(jd, id + kd) == 0)
						{
							for (int hd = 0; hd < kd; hd++)
							{
								srz.at<uchar>(jd, id + hd) = srx.at<uchar>(jd, id + hd);
							}
						}
						break;
					}
				}
			}
		}
	}
	
	std::vector<std::vector<Point2i>> cntz; // cntz are trimed parts from nut, it possibly contains branches.
	findContours(srz, cntz, RETR_LIST, CHAIN_APPROX_NONE);
	sort_contours_by_area(cntz); // sort the trimed contours, the larger the closer to first, so the branches might be the first.

	// headside < 0: head is left; > 0: head is right; == 0: unknown.
	headside = 0; floweroffset = 1000; // floweroffset is initialized to a value a normal nut flower can never reach.

	if ((int)cntz.size() > 0)
	{
		// find the center of the nut
		int idx[4] = {};
		find_the_left_right_points(cnts[0], idx);
		int xc = (cnts[0][idx[0]].x + cnts[0][idx[1]].x) / 2;

		// find the center of the flower
		int idz[4] = {};
		find_the_left_right_points(cntz[0], idz);
		int xz = (cntz[0][idz[0]].x + cntz[0][idz[1]].x) / 2;

		int dy = cntz[0][idz[3]].y - cntz[0][idz[2]].y;

		// if the flower height is larger than the threshold, we take it as a flower.
		if (dy > flowerheightthresh)
		{
			// check which side is the head.
			headside = (xz > xc ? -1 : +1);

			// to get a line, through the center of hpt and fpt, verticle to the line from hpt to fpt
			Point2i hpt = (headside < 0 ? cnt[idx[0]] : cnt[idx[1]]);
			Point2i fpt = (headside < 0 ? cnt[idx[1]] : cnt[idx[0]]);

			int dx = hpt.x - fpt.x, dy = hpt.y - fpt.y;
			Point2i cpt = (Point2i)((hpt + fpt) * 0.5);
			Point2i cpv = cpt + Point2i(dy, -dx);

			Point2i flower = (headside < 0 ? cntz[0][idz[0]] : cntz[0][idz[1]]); // the inner point of the flower contour
			floweroffset = (int)distance_from_point_to_line(flower, cpt, cpv); // floweroffset is always positive.
		}
	}

}



void Cutter::find_head_by_flower_neck(std::vector<Point2i>& cnt, int flowerheightthresh, double mblankscale, double mwidthscale, int& headside, int& floweroffset)
{
	// prerequisite: all the branches have been trimed, with or without flower.
	// output: headside: -1: head is left; +1: head is right; 0: unknown.
	// output: floweroffset: unit(pixel), from the center of cpt to flower, is always positive, no matter flower is at left or right of nut.
	//		   when no flower exist, flowerdist = 1000 by default.
	//
	int idx[4] = {};
	find_the_left_right_points(cnt, idx);

	int nutwidth = cnt[idx[1]].x - cnt[idx[0]].x, xc = (cnt[idx[0]].x + cnt[idx[1]].x) / 2;
	int mblank = (int)(nutwidth * mblankscale), mwidth = (int)(nutwidth * mwidthscale);

	int cntlen = (int)cnt.size(), neckinfo[4][3] = {};

	// from the most left and right point, searching forward and backward, 
	// from left to right, in the region [mblank, mblank+mwidth], find the highest
	// vertical distance of two points in 4 contour segments.
	for (int kd = 0; kd < 2; kd++) // left and right side
	{
		for (int rd = 0; rd < 2; rd++) // top and down side
		{
			for (int id = 0; id < (cntlen / 6); id++)
			{
				for (int jd = id + 1; jd < (cntlen / 4); jd++)
				{
					// sgn = -1: backward; sgn = 1: forward.
					int sgn = (rd == 0 ? 1 : -1);

					int ip = (idx[kd] + sgn * id + cntlen) % cntlen;
					int jp = (idx[kd] + sgn * jd + cntlen) % cntlen;

					int ix = cnt[ip].x - cnt[idx[kd]].x;
					int jx = cnt[jp].x - cnt[idx[kd]].x;

					ix = (kd == 0 ? ix : -ix);
					jx = (kd == 0 ? jx : -jx);

					if ((ix > mblank) && (jx < mblank + mwidth))
					{
						int dy = cnt[ip].y - cnt[jp].y;
						dy = (kd == 0 ? dy : -dy);
						dy = (rd == 0 ? dy : -dy);

						if (dy > neckinfo[2 * kd + rd][0])
						{
						neckinfo[2 * kd + rd][0] = dy; // the 1st column is flower height
						neckinfo[2 * kd + rd][1] = ip; // the 2nd column is outer point id;
						neckinfo[2 * kd + rd][2] = jp; // the 2nd column is inner point id;
						}
					}
				}
			}
		}
	}

	// the highest segment is the possible flower.
	int ix = 0;
	for (int id = 0; id < 4; id++)
	{
		if (neckinfo[id][0] > neckinfo[ix][0]) { ix = id; }
	}

	// if the possible segment is higher than the height thresh, then it is flower.
	headside = 0; floweroffset = 1000; // floweroffset is initialized to a value a normal nut flower can never reach.
	if (neckinfo[ix][0] > flowerheightthresh)
	{
		// -1: head is left; +1: head is right; 0: unknown.
		headside = (cnt[neckinfo[ix][1]].x > xc ? -1 : +1);

		// to get a line, through the center of hpt and fpt, verticle to the line from hpt to fpt
		Point2i hpt = (headside < 0 ? cnt[idx[0]] : cnt[idx[1]]);
		Point2i fpt = (headside < 0 ? cnt[idx[1]] : cnt[idx[0]]);

		int dx = hpt.x - fpt.x, dy = hpt.y - fpt.y;
		Point2i cpt = (Point2i)((hpt + fpt) * 0.5);
		Point2i cpv = cpt + Point2i(dy, -dx);

		// to get the distance from flower to the center vertical line of nut
		int ip = neckinfo[ix][1], jp = neckinfo[ix][2];
		
		floweroffset = (int)distance_from_point_to_line(cnt[jp], cpt, cpv); // floweroffset is always positive.
	}
}





// PART IV: detect images

bool Cutter::Detect2(Bitmap^ bgrImage)
{
	// the first process when the whole body of nut has came into sight.
	Mat1b grayMat = BgrBmp2GrayMat(bgrImage), binMat;

	std::vector<Point2i> cnt;
	if (!receive_image_find_nut_contour(grayMat, -1, (int)(grayThresh), cnt, binMat))
	{
		failureCode = 1;
		return false;
	}
	
	// 1. is nut lying between jig?
	cv::Rect rect = boundingRect(Mat(cnt));
	topY = rect.tl().y;
	btmY = rect.br().y;

	if (topY<=JigTopY || btmY>=JigBtmY)
	{
		failureCode = 4;
		return false;
	}
	
	//Point2i fp1, fp2;
	//Point2i center = FindCenter(cnt, fp1, fp2);
	//Center = ToPoint(center);

	// 2. finding head with 3 methods, and triming the branches.
	int headsides[] = { 0,0,0 };
	int floweroffset = 0;
	find_head_by_cluster_method(binMat, headsides[0]); 
	trim_branch_and_find_head(binMat, branchwidththresh, branchlengththresh, headsides[1]); // if branch exists, it has been removed, thus the contour should be updated.

	std::vector<std::vector<Point2i>> cntz;
	find_valid_nut_contours(binMat, cntz);
	cnt = cntz[0]; // cnt has been updated.
	find_head_by_flower_neck(cnt, flowerheightthresh, mblankscale, mwidthscale, headsides[2], floweroffset);

	headside = headsides[0]; // headside is a public variable.
	if (headsides[1] != 0) { headside = headsides[1]; }
	if (headsides[2] != 0) { headside = headsides[2]; }

	// 3. finding head and foot points
	Point2i ct1, ct2, head, foot;
	find_the_left_right_centers(cnt, ct1, ct2);
	head = (headside < 0 ? ct1 : ct2);
	foot = (headside < 0 ? ct2 : ct1);

	// 4. finding rectangle
	RotatedRect rRect = find_nut_rectangle(cnt, 20); // this command should wait until cnt has been updated.

	RectPnts = RectForShow(rRect); // show rectangle

	//=============================================
	Head = ToPoint(head);
	Center = ToPoint(rRect.center);
	Contour = ContourForShow(cnt);
	//=============================================

	// 5. checking moon
	length = max(rRect.size.width, rRect.size.height)*dimPerPix;
	height = min(rRect.size.width, rRect.size.height)*dimPerPix;

	//moon = FindMoon(cnt, head, foot);
	double moonx = 0;
	get_nut_bending_tendency(cnt, headside, moonx);
	moon = moonx;

	double moonThresh = height<boundHeight && length<boundLength ? moonThresh1 : moonThresh2;
	if (moon<moonThresh)
	{
		failureCode = 3;
		return false;
	}
	if (length>overMaxLength || height<overMinHeight)
	{
		failureCode = 5;
		return false;
	}
	failureCode = 0;
	return true;
}



bool Cutter::Detect4(Bitmap^ bgrImage1, Bitmap^ bgrImage2)
{

	// the first process when the whole body of nut has came into sight.
	Mat1b src = BgrBmp2GrayMat(bgrImage2);

	std::vector<Point2i> cnt;

	if (!receive_image_find_nut_contour(src, 1, (int)(grayThresh), cnt, src))
	{
		failureCode = 1;
		return false;
	}

	// 2. finding head with 3 methods, and triming the branches.
	int headsides[] = { 0,0,0 };
	
	find_head_by_cluster_method(src, headsides[0]);
	trim_branch_and_find_head(src, branchwidththresh, branchlengththresh, headsides[1]); // if branch exists, it has been removed, thus the contour should be updated.

	std::vector<std::vector<Point2i>> cntz;
	find_valid_nut_contours(src, cntz);
	cnt = cntz[0]; // cnt has been updated.
	int floweroffset = 0;
	find_head_by_flower_neck(cnt, flowerheightthresh, mblankscale, mwidthscale, headsides[2], floweroffset);

	headside = headsides[0]; // headside is a public variable.
	if (headsides[1] != 0) { headside = headsides[1]; }
	if (headsides[2] != 0) { headside = headsides[2]; }

	// 3. finding head and foot points
	Point2i ct1, ct2, head, foot;
	find_the_left_right_centers(cnt, ct1, ct2);
	head = (headside < 0 ? ct1 : ct2);
	foot = (headside < 0 ? ct2 : ct1);

	// 4. finding rectangle
	RotatedRect rRect = find_nut_rectangle(cnt, 20); // this command should wait until cnt has been updated.

	RectPnts = RectForShow(rRect); // show rectangle

	//=============================================
	Head = ToPoint(head);
	Center = ToPoint(rRect.center);
	Contour = ContourForShow(cnt);
	//=============================================

	length = max(rRect.size.width, rRect.size.height)*dimPerPix;
	height = min(rRect.size.width, rRect.size.height)*dimPerPix;
	if (length>overMaxLength || height<overMinHeight)
	{
		failureCode = 5;
		return false;
	}

	if ((height<boundHeight && length<boundLength ? OneCut(cnt, head, foot, CutDatas) : TwoCut(cnt, head, foot, CutDatas)))
	{
		failureCode = 0;
		return true;
	}
	else
	{
		failureCode = 7;
		return false;
	}
}