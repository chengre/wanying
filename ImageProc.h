#pragma once

#pragma comment(lib, "C:/OpenCV2.1.0/lib/cv210.lib")
#pragma comment(lib, "C:/OpenCV2.1.0/lib/cxcore210.lib")
#pragma comment(lib, "C:/OpenCV2.1.0/lib/highgui210.lib")

#pragma comment(lib, "D2D1")  
#pragma comment(lib, "DWrite")

#define _USE_MATH_DEFINES
#define _CRT_NON_CONFORMING_SWPRINTFS

#include "C:/OpenCV2.1.0/include/cv.h"
#include "C:/OpenCV2.1.0/include/cxcore.h"
#include "C:/OpenCV2.1.0/include/highgui.h"

#include <D2D1.h>
#include <DWrite.h>   
#include <math.h>
#include <tchar.h>

namespace ImageProc
{
	using namespace cv;
	using namespace D2D1;
	using namespace System;
	using namespace System::IO;
	using namespace System::Windows::Forms;
	using namespace System::Drawing;
	using namespace System::Drawing::Imaging;
	using namespace System::Collections::Generic;
	using namespace System::Runtime::InteropServices;

	template<typename Type>  
	void SafeRelease(Type& p)  
	{  
		if(p)  
		{  
			p->Release();  
			p = 0;  
		}  
	}

	class pRender
	{
	public:
		ID2D1Factory* factory;
		ID2D1HwndRenderTarget* renderTarget;
		ID2D1SolidColorBrush* solidBrush;
		ID2D1Bitmap* bitmap;
		D2D1_HWND_RENDER_TARGET_PROPERTIES hwndRenderTargetProperties;
		D2D1_BITMAP_PROPERTIES bitmapProperties;
		D2D1_SIZE_U bitmapSize;
		D2D1_RECT_F renderRect;
	};

	public ref class ColorSpyEventArgs : EventArgs
	{
	public:
		int X;
		int Y;
		byte B;
		byte G;
		byte R;
		byte Gray;
		bool InRange;
		void Locate(IntPtr bgrPtr, int width, int height)
		{
			if (cv::Rect(1, 1, width, height).contains(Point2i(X, Y)))
			{
				Vec3b bgr = Mat(height, width, CV_8UC(3), (byte*)bgrPtr.ToPointer()).at<Vec3b>(height-Y-1, X);
				B = bgr[0];
				G = bgr[1];
				R = bgr[2];
				Gray = (R*38 + G*75 + B*15) >> 7;
			}
		}
		System::String^ ToStr()
		{
			return !InRange ? System::String::Empty :
				"U:" + X.ToString() + " " +
				"V:" + Y.ToString() + " " +
				"红:" + R.ToString() + " " +
				"绿:" + G.ToString() + " " +
				"蓝:" + B.ToString() + " " +
				"灰:" + Gray.ToString();
		}
	};

	public ref class Render : Control
	{
	public:
		EventHandler<ColorSpyEventArgs^>^ ColorSpyChanged;
		pRender* base;
		byte *bgra, *flip;
		float s;
		Render()
		{
			base = new pRender();
			D2D1CreateFactory(D2D1_FACTORY_TYPE_MULTI_THREADED, &base->factory);
			base->hwndRenderTargetProperties = HwndRenderTargetProperties((HWND)this->Handle.ToPointer());
			base->bitmapProperties = D2D1::BitmapProperties(D2D1::PixelFormat(DXGI_FORMAT_B8G8R8A8_UNORM, D2D1_ALPHA_MODE_IGNORE));
			bgra = new byte[640*180*4];
			flip = new byte[640*180*4];
			base->bitmapSize = SizeU(640,180);
			Mat1b(1, 640*180*4, bgra) *= 0;
			InitRenderTarget();
		}
		~Render()
		{
			SafeRelease(base->renderTarget);
			SafeRelease(base->solidBrush);
			SafeRelease(base->bitmap);
			SafeRelease(base->factory);
			delete base;
			delete []bgra; 
			delete []flip; 
		}
		bool InitRenderTarget()
		{
			base->hwndRenderTargetProperties.pixelSize = SizeU(this->ClientSize.Width, this->ClientSize.Height);
			HRESULT hr = base->factory->CreateHwndRenderTarget(RenderTargetProperties(), base->hwndRenderTargetProperties, &base->renderTarget);
			if (SUCCEEDED(hr))
			{
				base->renderTarget->CreateSolidColorBrush(ColorF(1.0f,0,0), &base->solidBrush);
				base->renderTarget->CreateBitmap(base->bitmapSize, base->bitmapProperties, &base->bitmap);
				base->bitmap->CopyFromMemory(&RectU(0, 0, base->bitmapSize.width, base->bitmapSize.height), bgra, 4*base->bitmapSize.width);
			}
			return SUCCEEDED(hr);
		}
		void SetParent(Control^ parent)
		{
			this->Location = System::Drawing::Point::Empty;
			this->Size = parent->Size;
			this->Anchor = (AnchorStyles)15;
			parent->Controls->Add(this);
		}
		void SetImage(IntPtr bgrPtr, int w, int h, bool needFlip)
		{
			int iw = base->bitmapSize.width;
			int ih = base->bitmapSize.height;
			if (iw != w || ih != h)
			{
				base->bitmapSize = SizeU(w, h);
				SafeRelease(base->bitmap);
				base->renderTarget->CreateBitmap(base->bitmapSize, base->bitmapProperties, &base->bitmap);
				delete []bgra;
				delete []flip;
				bgra = new byte[w*h*4];
				flip = new byte[w*h*4];
			}
			if (needFlip)
			{
				cvtColor(Mat(h, w, CV_8UC3, bgrPtr.ToPointer()), Mat(h, w, CV_8UC4, flip), CV_BGR2BGRA);
				cv::flip(Mat(h, w, CV_8UC4, flip), Mat(h, w, CV_8UC4, bgra), 0);
			}
			else
			{
				cvtColor(Mat(h, w, CV_8UC3, bgrPtr.ToPointer()), Mat(h, w, CV_8UC4, bgra), CV_BGR2BGRA);
			}
			base->bitmap->CopyFromMemory(&RectU(0, 0, w, h), bgra, 4*w);
		}
		void SetImage(Bitmap^ bgrImg)
		{
			BitmapData^ bmpData = bgrImg->LockBits(System::Drawing::Rectangle(0, 0, bgrImg->Width, bgrImg->Height), ImageLockMode::ReadWrite, bgrImg->PixelFormat);
			SetImage(bmpData->Scan0, bgrImg->Width, bgrImg->Height, false);
			bgrImg->UnlockBits(bmpData);
		}
		virtual void OnPaintBackground(PaintEventArgs^ e)override
		{

		}
		virtual void OnSizeChanged(EventArgs^ e)override
		{
			Draw();
		}
		virtual void OnPaint(PaintEventArgs^ e)override
		{
			Draw();
		}
		void Draw()
		{
			if(!base->renderTarget && !InitRenderTarget())  
			{
				return;
			}

			int cw = this->ClientSize.Width;
			int ch = this->ClientSize.Height;
			D2D1_SIZE_U ps = base->renderTarget->GetPixelSize();
			if (ps.width!=cw || ps.height!=ch)
			{
				base->renderTarget->Resize(SizeU(cw, ch));
			}

			float iw = float(base->bitmapSize.width);
			float ih = float(base->bitmapSize.height);
			if (float(cw)/ch<iw/ih)
			{
				base->renderRect = RectF(0, (ch - cw*ih/iw)/2, float(cw), (ch + cw*ih/iw)/2);
				s = cw/iw;
			}
			else
			{
				base->renderRect = RectF((cw - ch*iw/ih)/2, 0, (cw + ch*iw/ih)/2, float(ch));
				s = ch/ih;
			}

			base->renderTarget->BeginDraw();
			base->renderTarget->Clear(ColorF(ColorF::DarkSeaGreen));
			base->renderTarget->DrawBitmap(base->bitmap, base->renderRect);

			if (D2DERR_RECREATE_TARGET == base->renderTarget->EndDraw())
			{
				SafeRelease(base->renderTarget);
				SafeRelease(base->solidBrush);
				SafeRelease(base->bitmap);
			}
		}
		virtual void OnMouseMove(MouseEventArgs^ e)override
		{
			__super::OnMouseMove(e);
			if (ColorSpyChanged != nullptr)
			{
				ColorSpyEventArgs^ le = gcnew ColorSpyEventArgs();
				int x = cvRound((e->X - base->renderRect.left)/s);
				int y = cvRound((e->Y - base->renderRect.top)/s);
				int iw = base->bitmapSize.width;
				int ih = base->bitmapSize.height;
				le->InRange = cv::Rect(0, 0, iw, ih).contains(Point2i(x, y));
				if (le->InRange)
				{
					le->X = x;
					le->Y = y;
					le->B = bgra[(y*iw+x)*4 + 0];
					le->G = bgra[(y*iw+x)*4 + 1];
					le->R = bgra[(y*iw+x)*4 + 2];
					le->Gray = (le->R*38 + le->G*75 + le->B*15) >> 7;
				}
				ColorSpyChanged(this, le);
			}
		}
		virtual void OnMouseLeave(EventArgs^ e)override
		{
			__super::OnMouseLeave(e);
			if (ColorSpyChanged != nullptr)
			{
				ColorSpyEventArgs^ le = gcnew ColorSpyEventArgs();
				le->InRange = false;
				ColorSpyChanged(this, le);
			}
		}
	};

	public ref class CutData
	{
	public:
		int Flag;
		Color ColorForShow;
		double Length;
		double Height;
		double LengthScore;
		double HeightScore;
		System::Drawing::PointF CutPnt1;
		System::Drawing::PointF CutPnt2;
		double Ang;
		double Dist;
		///////////////////////////////////////////
		template<typename T> 
		CutData(double l, double h, double ls, double hs, Point_<T> &cp1, Point_<T> &cp2, int flag)
		{
			Length = l;
			Height = h;
			LengthScore = ls;
			HeightScore = hs;
			CutPnt1 = System::Drawing::PointF(float(cp1.x), float(cp1.y));
			CutPnt2 = System::Drawing::PointF(float(cp2.x), float(cp2.y));
			Flag = flag;
			if (flag==0)
			{
				ColorForShow = Color::White;
			}
			else if (flag==1)
			{
				ColorForShow = Color::Yellow;
			}
			else if (flag==2)
			{
				ColorForShow = Color::Magenta;
			}
		}
		CutData()
		{
			LengthScore = 0;
			HeightScore = 0;
		}
		double SumScore()
		{
			return LengthScore+HeightScore;
		}
		System::String^ ToStr()
		{
			System::String^ str = "角度：" + Ang.ToString("0.00") + " 距离：" + Dist.ToString("0.00") + "\r\n切长：" + Length.ToString("0.00") + " 背高：" + Height.ToString("0.00");
			if (Flag>0)
			{
				str+= "\r\n切长分：" + LengthScore.ToString("0.0") + " 背高分：" + HeightScore.ToString("0.0") + "\r\n";
			}
			return str;
		}
	};

	public ref class Cutter
	{
	public:
		/////////////////////输入
		double grayThresh;
		double minLength;
		double minHeight;

		double moonRetPct;
		double moonThresh1;
		double moonThresh2;

		double dimPerPix;//

		double matchDeadX1;//
		double matchDeadX2;//
		double matchThresh;

		double overMaxLength;
		double overMinHeight;

		double boundLength;
		double boundHeight;

		double oneCutHeadRet;
		double oneCutFootRet;
		double oneCutMinH;

		double twoCutHeadRet;
		double twoCutHeadOfsThresh;
		double twoCutHeadOfs1;
		double twoCutHeadOfs2;
		double twoCutFootRet;
		double twoCutFootOfs;
		double twoCutFootStep;

		double twoCutMinLength;
		double twoCutMaxLength;
		double twoCutMinHeight;
		double twoCutMaxHeight;
		double twoCutMode;

		double midRangePct;

		array<double>^ lengthList1;
		array<double>^ heightList1;
		array<double>^ lengthList2;
		array<double>^ heightList2;
		array<double>^ lengthScoreList1;
		array<double>^ heightScoreList1;
		array<double>^ lengthScoreList2;
		array<double>^ heightScoreList2;

		double JigTopY;// 夹具上坐标
		double JigBtmY;// 夹具下坐标
		double JigX;// 夹具中心x
		double JigY;// 夹具中心y
		double CutX;// 切刀中心x
		double CutY;// 切刀中心y

		double FlagCrcR;
		double FlagSft1;
		double FlagSft2;

		int branchwidththresh;
		int branchlengththresh;
		int flowerwidththresh;
		int flowerheightthresh;
		double mblankscale;
		double mwidthscale;
		int fStepWd;
		int fThresh;

		int headside; // headside = -1: head is left; headside = +1: head is right; 
		//int sideoffset; // the side offset around the most left/right points to obtain two shoulder points for each, the two center points from left to right can be derived to get the nut direction.
		//int rectoffset; // the side offset for the new 


		/////////////////////输出
		array<System::Drawing::Point>^ Contour;
		System::Drawing::Point Head;
		array<System::Drawing::PointF>^ RectPnts;
		System::Drawing::Point Center;
		array<CutData^>^ CutDatas;

		System::Drawing::PointF HeadOfsPnt1;
		System::Drawing::PointF HeadLimPnt1;
		System::Drawing::PointF HeadCtrPnt1;
		System::Drawing::PointF FootOfsPnt1;
		System::Drawing::PointF FootCtrPnt1;

		System::Drawing::PointF HeadOfsPnt2;
		System::Drawing::PointF HeadLimPnt2;
		System::Drawing::PointF HeadCtrPnt2;
		System::Drawing::PointF FootOfsPnt2;
		System::Drawing::PointF FootCtrPnt2;

		System::Drawing::PointF FlagCrcPnt1;
		System::Drawing::PointF FlagCrcPnt2;
		System::Drawing::PointF FlagCutPnt1;
		System::Drawing::PointF FlagCutPnt2;

		double length;
		double height;
		double topY;
		double btmY;
		double moon;
		double match;
		int failureCode;
		/////////////////////

		// declaration of functions.

		// PART I: function with contour as input.

		void find_the_left_right_points(const std::vector<Point2i>& cnt, int idx[4]);
		void find_the_left_right_centers(const std::vector<Point2i>& cnt, Point2i& ct1, Point2i& ct2);
		void find_the_left_right_centers(const std::vector<Point2i>& cnt, int indentleft, int indentright, Point2i& ct1, Point2i& ct2);
		RotatedRect find_nut_rectangle(const std::vector<Point2i>& cnt, int dwidth);

		void trim_contour_with_offset(std::vector<Point2i>& cnt, int dwidth);
		std::vector<Point2i> convert_contour2f_to_contour2i(std::vector<Point2f> contour2f);
		std::vector<Point2f> convert_contour2i_to_contour2f(std::vector<Point2i> contour2i);
		void rotate_contour_by_angle(std::vector<Point2f>&, RotatedRect rRect, float angle);
		//double norm_point(Point2i pt);

		void sort_contours_by_area(std::vector <std::vector<Point2i>>& cnts);
		void sort_contours_by_length(std::vector <std::vector<Point2i>>& cnts);
		double distance_from_point_to_line(Point2i pt, Point2i pt1, Point2i pt2);
		void get_nut_bending_tendency(std::vector<Point2i>& cnt, int headside, double& indicator);

		bool get_length_of_splited_contour(std::vector<Point2i>& cnt, Point2i pti, Point2i ptj, int headside, int floweroffset, int flowerindent, int idx[2], double& xlength);
		bool get_heights_of_splited_contour(std::vector<Point2i>& cnt, Point2i pt1, Point2i pt2, int headside, int floweroffset, double xheights[2]);

		bool remove_branches(Mat1b& grayMat, Mat1b& ctrMat, std::vector<Point2i>& contour, int branchwidththresh);

		bool FindHead(Mat1b, std::vector<Point2i>&, Point2i&, Point2i&);
		void FindFoot(std::vector<Point2i> &, Point2i &, Point2i &);

		// PART II: functions with image as input.

		bool receive_image_find_nut_contour(Mat1b& src, int nutside, int graythresh, std::vector<Point2i>& cnt, Mat1b& srx);
		void convert_colored_to_binary_image(Mat1b& src, Mat1b& srx, int graythresh);
		bool find_valid_nut_contours(Mat1b& srx, std::vector<std::vector<Point2i>>& cnts);

		void find_head_by_cluster_method(Mat1b& srx, int& headside);
		void trim_branch_and_find_head(Mat1b& srx, int branchwidththresh, int branchlengththresh, int& headside);
		void find_flower_and_find_head(Mat1b& srx, int flowerwidththresh, int flowerheightthresh, double mblankscale, double mwidthscale, int& headside, int& floweroffset);
		void find_head_by_flower_neck(std::vector<Point2i>& cnt, int flowerheightthresh, double mblankscale, double mwidthscale, int& headside, int& floweroffset);

		bool Detect2(Bitmap^ bgrImage);
		bool Detect4(Bitmap^ bgrImage1, Bitmap^ bgrImage2);
		// end of declaration.

		bool Detect1(Bitmap^ bgrImage)
		{
			// check if the whole body of the nut has came into the sight of camera.
			length = 0;
			height = 0;
			topY = 0;
			btmY = 0;
			moon = 0;
			match = 0;
			vector<Point2i> contour;
			if (FindContour(BgrBmp2GrayMat(bgrImage), true, contour))
			{
				failureCode = 0;
				return true;
			}
			else
			{
				failureCode = 1;
				return false;
			}
		}


		/*
		bool Detect2(Bitmap^ bgrImage)
		{
			// the first process when the whole body of nut has came into sight.
			Mat1b grayMat = BgrBmp2GrayMat(bgrImage), binMat;
			vector<Point2i> contour;
			if (!FindContour(grayMat, true, contour, binMat))
			{
				failureCode = 1;
				return false;
			}

			cv::Rect rect = boundingRect(Mat(contour));
			topY = rect.tl().y;
			btmY = rect.br().y;

			Point2i fp1, fp2;
			Point2i center = FindCenter(contour, fp1, fp2);
			Center = ToPoint(center);

			if (topY<=JigTopY || btmY>=JigBtmY)
			{
				failureCode = 4;
				return false;
			}

			// find head and foot with trimed contour
			std::vector<Point2i> contourx = contour;
			trim_contour_with_offset(contourx, 12);
			Point2i head, foot;

			if (!FindHead(grayMat, contourx, head, foot))
			{
				failureCode = 2;
				return false;
			}
			Head = ToPoint(head);
			/////////////////////////////////////////////
			//FixContour(binMat, contour, head);

			Contour = ContourForShow(contourx);
			/////////////////////////////////////////////
			//Point2i foot;

			// RotatedRect rRect = FindRect(contour, head, foot); // called by yang

			// finding the foot point. by chen
			//FindFoot(contour, head, foot); // finding the foot point by choosing the farest point to the head point.

			// finding the rectangle of the nut with offset. by chen
			//Point2f keypoints[2];
			//find_key_points_with_offset(contour, keypoints, 20);
			//RotatedRect rRect = FindRect(contour, keypoints, 0);

			Point2i ct1, ct2;
			find_the_left_right_centers(contour, ct1, ct2);
			Point2i ct = ct2 - ct1;
			double dl = norm(ct);
			double vech[2] = { ct.x / dl,ct.y / dl };
			RotatedRect rRect = find_nut_rectangle(contour, vech);
			// end of finding.


			RectPnts = RectForShow(rRect);
			length = max(rRect.size.width, rRect.size.height)*dimPerPix;
			height = min(rRect.size.width, rRect.size.height)*dimPerPix;
			moon = FindMoon(contour, head, foot);

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
		*/



		bool Detect3(Bitmap^ bgrImage)
		{
			Mat1b grayMat = BgrBmp2GrayMat(bgrImage), binMat;
			vector<Point2i> contour;
			if (!FindContour(grayMat, false, contour, binMat))
			{
				failureCode = 1;
				return false;
			}

			cv::Rect rect = boundingRect(Mat(contour));
			topY = rect.tl().y;
			btmY = rect.br().y;

			Point2i fp1, fp2;
			Point2i center = FindCenter(contour, fp1, fp2);
			Center = ToPoint(center);

			if (topY<=JigTopY || btmY>=JigBtmY)
			{
				failureCode = 4;
				return false;
			}

			if (matchDeadX1-rect.tl().x<3 || rect.br().x-matchDeadX2<3)
			{
				failureCode = 6;
				return false;
			}

			// find head and foot with trimed contour
			std::vector<Point2i> contourx = contour;
			trim_contour_with_offset(contourx, 12);
			Point2i head, foot;

			if (!FindHead(grayMat, contourx, head, foot))
			{
				failureCode = 2;
				return false;
			}
			Head = ToPoint(head);
			/////////////////////////////////////////////
			//FixContour(binMat, contour, head);

			Contour = ContourForShow(contourx);
			/////////////////////////////////////////////
			//Point2i foot;
			//RotatedRect rRect = FindRect(contour, head, foot);

			// finding the foot point. by chen
			//FindFoot(contour, head, foot); // finding the foot point by choosing the farest point to the head point.

			// finding the rectangle of the nut with offset. by chen
			//Point2f keypoints[2];
			//find_key_points_with_offset(contour, keypoints, 20);
			//RotatedRect rRect = FindRect(contour, keypoints, 10);
			Point2i ct1, ct2;
			find_the_left_right_centers(contour, ct1, ct2);
			Point2i ct = ct2 - ct1;
			double dl = norm(ct);
			double vech[2] = { ct.x / dl,ct.y / dl };
			RotatedRect rRect = find_nut_rectangle(contour, 20);
			// end of finding.

			RectPnts = RectForShow(rRect);
			length = max(rRect.size.width, rRect.size.height)*dimPerPix;
			height = min(rRect.size.width, rRect.size.height)*dimPerPix;
			if (length>overMaxLength || height<overMinHeight)
			{
				failureCode = 5;
				return false;
			}
			if (height<boundHeight && length<boundLength ? OneCut(contour, head, foot, CutDatas) : TwoCut(contour, head, foot, CutDatas))
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


		/*
		bool Detect4(Bitmap^ bgrImage1, Bitmap^ bgrImage2)
		{
			vector<Point2i> contour1, contour2, contour;
			Mat1b m1 = BgrBmp2GrayMat(bgrImage1), m2 = BgrBmp2GrayMat(bgrImage2);
			if (!FindContour(m1, false, contour1) || !FindContour(m2, false, contour2))
			{
				failureCode = 1;
				return false;
			}
			// trim contours with width offset, 6.16
			//trim_contour_with_offset(contour1, 10);
			//trim_contour_with_offset(contour2, 10);
			// trim end

			cv::Rect rect1 = boundingRect(Mat(contour1));
			cv::Rect rect2 = boundingRect(Mat(contour2));
			double l1 = max(rect1.width, rect1.height);
			double l2 = max(rect2.width, rect2.height);
			if (l1/l2<0.75 || l2/l1<0.75)
			{
				failureCode = 1;
				return false;
			}

			Mat1b binMat;
			match = Match(m1.size(), contour1, contour2, contour, binMat);
			if (match<matchThresh)
			{
				failureCode = 6;
				return false;
			}

			Point2i fp1, fp2;
			Point2i center = FindCenter(contour, fp1, fp2);
			Center = ToPoint(center);

			// find head and foot with trimed contour
			std::vector<Point2i> contourx = contour;
			trim_contour_with_offset(contourx, 12);
			Point2i head, foot;
			//if (!FindHead(binMat, contour, boundingRect(Mat(contour)), fp1, fp2, head))
			if (!FindHead(binMat, contourx, head, foot))
			{
				failureCode = 2;
				return false;
			}
			Head = ToPoint(head);
			/////////////////////////////////////////////
			//FixContour(binMat, contour, head);

			Contour = ContourForShow(contourx);
			/////////////////////////////////////////////
			//Point2i foot;
			//RotatedRect rRect = FindRect(contour, head, foot);

			// finding the foot point. by chen
			//FindFoot(contour, head, foot); // finding the foot point by choosing the farest point to the head point.

			// finding the rectangle of the nut with offset. by chen
			//Point2f keypoints[2];
			//find_key_points_with_offset(contour, keypoints, 20);
			//RotatedRect rRect = FindRect(contour, keypoints, 30);
			Point2i ct1, ct2;
			find_the_left_right_centers(contour, ct1, ct2);
			Point2i ct = ct2 - ct1;
			double dl = norm(ct);
			double vech[2] = { ct.x / dl,ct.y / dl };
			RotatedRect rRect = find_nut_rectangle(contour, 20);
			// end of finding.

			RectPnts = RectForShow(rRect);
			length = max(rRect.size.width, rRect.size.height)*dimPerPix;
			height = min(rRect.size.width, rRect.size.height)*dimPerPix;
			if (length>overMaxLength || height<overMinHeight)
			{
				failureCode = 5;
				return false;
			}

			if ((height<boundHeight && length<boundLength ? OneCut(contour, head, foot, CutDatas) : TwoCut(contour, head, foot, CutDatas)) && Check(m2))
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
		*/


		/////////////////////////////////////////////////////////////////
		Mat1b BgrBmp2GrayMat(Bitmap^ bgrImage)
		{
			Mat1b grayMat;
			BitmapData^ bmpData = bgrImage->LockBits(System::Drawing::Rectangle(0,0,bgrImage->Width,bgrImage->Height), ImageLockMode::ReadWrite, PixelFormat::Format24bppRgb);
			cvtColor(Mat(bgrImage->Height, bgrImage->Width, CV_8UC(3), (byte*)bmpData->Scan0.ToPointer(), bmpData->Stride), grayMat, CV_BGR2GRAY);
			bgrImage->UnlockBits(bmpData);
			return grayMat;
		}
		bool FindContour(Mat1b& grayMat, bool isLeft, vector<Point2i> &contour, Mat1b& ctrMat)
		{
			// convert grayscale to binary mat, and find the contours.
			Mat1b binMat;
			cv::threshold(grayMat, binMat, grayThresh, 255, THRESH_BINARY+THRESH_BINARY_INV);
			// the above type THRESH_BINARY_INV should not be omited.

			Mat se = getStructuringElement(MORPH_RECT, cv::Size(3, 3));
			erode(binMat, binMat, se);
			dilate(binMat, binMat, se);

			std::vector<std::vector<cv::Point>> contours;
			findContours(binMat, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);

			// find the nut contour on left or right side.
			double maxA = 0;
			int posi = -1;
			for (int i=0; i<(int)contours.size(); i++)
			{
				cv::Rect rect = boundingRect(Mat(contours[i]));
				if (rect.x>2 && rect.x+rect.width<=binMat.cols-2)
				{
					RotatedRect rRect = minAreaRect(Mat(contours[i]));
					float h = min(rRect.size.width, rRect.size.height);
					float l = max(rRect.size.width, rRect.size.height);
					float mx = rRect.center.x;
					if (h>minHeight/dimPerPix && l<overMaxLength/dimPerPix && l>minLength/dimPerPix)
					{
						if ((isLeft && mx<binMat.cols/2.0) || (!isLeft && mx>binMat.cols/2.0))
						{
							if (h*l>=maxA)
							{
								maxA = h*l;
								posi = i;
							}
						}
					}
				}
			}

			if (maxA==0)
			{
				return false;
			}


			// generate a black image with a white nut in, and return this image and contour.
			Mat1b fstMat = Mat1b::zeros(grayMat.size());
			drawContours(fstMat, contours, posi, Scalar(1), CV_FILLED);
			////////////////////////////////////////////////////
			/*
			ctrMat = fstMat;
			contour = contours[posi];
			*/

			//////////////////////////////////////////////// trim branches by chen.
			if(!remove_branches(fstMat, ctrMat, contour, branchwidththresh))
			{
				ctrMat = fstMat;
				contour = contours[posi];
			}

			//////////////////////////////////////////////////////////////////去花枝;
			// 			Mat1b tmpMat = fstMat.clone(); // nut without branches.
			// 			erode (tmpMat, tmpMat, getStructuringElement(MORPH_RECT, cv::Size(branchwidththresh, branchwidththresh)));
			// 			dilate(tmpMat, tmpMat, getStructuringElement(MORPH_RECT, cv::Size(branchwidththresh, branchwidththresh)));
			// 
			// 			Mat1b subMat = Mat1b::zeros(grayMat.size()); // branches.
			// 			for (int i=0;i<grayMat.rows;i++)
			// 			{
			// 				for (int j=0;j<grayMat.cols;j++)
			// 				{
			// 					if (fstMat(i,j)>0 && tmpMat(i,j)==0)
			// 					{
			// 						subMat(i,j) = 255;
			// 					}
			// 				}
			// 			}
			// 
			// 			findContours(subMat, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
			// 			for (int i=0; i<(int)contours.size(); i++)
			// 			{
			// 				RotatedRect rRec = minAreaRect(Mat(contours[i]));
			// 				if (max(rRec.size.width,rRec.size.height)<=branchwidththresh)
			// 				{
			// 					drawContours(tmpMat, contours, i, Scalar(255), CV_FILLED); // pots smaller than branches.
			// 				}
			// 			}
			// 
			// 			Mat1b tagMat = Mat1b::zeros(grayMat.size()); //nut and small pots, without branches.
			// 			for (int i=0;i<grayMat.rows;i++)
			// 			{
			// 				for (int j=0;j<grayMat.cols;j++)
			// 				{
			// 					if (fstMat(i,j)>0 && tmpMat(i,j)>0)
			// 					{
			// 						tagMat(i,j) = 255;
			// 					}
			// 				}
			// 			}
			// 
			// 			findContours(tagMat, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
			// 			double maxB = 0;
			// 			int posj = -1;
			// 			for (int j=0; j<(int)contours.size(); j++)
			// 			{
			// 				RotatedRect rRec = minAreaRect(Mat(contours[j]));
			// 				if (rRec.size.area()>=maxB)
			// 				{
			// 					maxB = rRec.size.area();
			// 					posj = j;
			// 				}
			// 			}
			// 			contour = contours[posj];
			// 			ctrMat = Mat1b::zeros(grayMat.size()); // only nut, without small pots and branches.
			// 			drawContours(ctrMat, contours, posj, Scalar(255), CV_FILLED);
			return true;
		}
		bool FindContour(Mat1b& grayMat, bool isLeft, vector<Point2i> &contour)
		{
			Mat1b ctrMat;
			return FindContour(grayMat, isLeft, contour, ctrMat);
		}
		void FixContour(Mat1b& ctrMat, vector<Point2i> &contour, Point2i &head)
		{
			RotatedRect rect = minAreaRect(Mat(contour));
			Point2f center = rect.center;
			float w = max(rect.size.width, rect.size.height);
			float h = min(rect.size.width, rect.size.height);

			Point2f pts1[4], pts2[4];
			rect.points(pts1);
			if (norm(pts1[0]-pts1[1])<norm(pts1[1]-pts1[2]))
			{
				pts2[0] = center + Point2f(-w/2, h/2);
				pts2[1] = center + Point2f(-w/2,-h/2);
				pts2[2] = center + Point2f( w/2,-h/2);
				pts2[3] = center + Point2f( w/2, h/2);
			}
			else
			{				
				pts2[0] = center + Point2f(-w/2,-h/2);
				pts2[1] = center + Point2f( w/2,-h/2);
				pts2[2] = center + Point2f( w/2, h/2);
				pts2[3] = center + Point2f(-w/2, h/2);
			}
			Mat1f trs1 = getAffineTransform(pts1,pts2);
			Mat1b tmpMat;
			cv::warpAffine(ctrMat, tmpMat, trs1, ctrMat.size());
			int tx = cvRound(trs1(0,0)*head.x + trs1(0,1)*head.y + trs1(0,2));
			int cx = cvRound(center.x);

			int fStep = 2*(fStepWd/2)+1;
			int fCount = (int)floor(0.5*w)/fStep;
			int fDrect = tx>cx ? -1 : 1;
			vector<double> aVec;
			for (int i=0;i<fCount;i++)
			{
				int mx = cx+fDrect*i*fStep;
				double sa = 0;
				for (int j=mx-fStep/2;j<=mx+fStep/2;j++)
				{
					sa += sum(tmpMat.col(j))[0];
				}
				aVec.push_back(sa/fStep);
			}

			int posi = 0;
			double maxA = 0;
			for (int i=0;i<(int)aVec.size()-1;i++)
			{
				for (int j=i+1;j<(int)aVec.size();j++)
				{
					if (aVec[j]-aVec[i]>=maxA)
					{
						posi = cx+fDrect*i*fStep;
						maxA = aVec[j]-aVec[i];
					}
				}
			}
			if (maxA<fThresh)
			{
				return;
			}

			Point2i p1, p2;
			for (int i=0;i<tmpMat.rows;i++)
			{
				if(tmpMat(i,posi)>0)
				{
					p1 = Point2i(posi,i);
					break;
				}
			}
			for (int i=tmpMat.rows-1;i>=0;i--)
			{
				if(tmpMat(i,posi)>0)
				{
					p2 = Point2i(posi,i);
					break;
				}
			}

			Point2i q1((cx+posi)/2,0), q2((cx+posi)/2,0);
			for (int i=cx-fStep/2;i<=cx+fStep/2;i++)
			{
				for (int j=0;j<tmpMat.rows;j++)
				{
					if (tmpMat(j,i)>0)
					{
						q1.y+=j;
						break;
					}
				}
				for (int j=tmpMat.rows-1;j>=0;j--)
				{
					if (tmpMat(j,i)>0)
					{
						q2.y+=j;
						break;
					}
				}
			}
			q1.y/=fStep;
			q2.y/=fStep;

			double params1[4], params2[4];
			GetLineParams(p1, q1, params1);
			GetLineParams(p2, q2, params2);

			int px = cvRound(cx+fDrect*w/2);
			Point2d r1(px, (-params1[0]*px-params1[2])/params1[1]);
			Point2d r2(px, (-params2[0]*px-params2[2])/params2[1]);

			vector<Point2i> ctr;
			ctr.push_back(p1);
			ctr.push_back(Point2i(r1));
			ctr.push_back(Point2i(r2));
			ctr.push_back(p2);
			vector<vector<Point2i>> ctrs; ctrs.push_back(ctr);
			Mat1b tagMat = Mat1b::zeros(tmpMat.size());
			drawContours(tagMat, ctrs, 0, Scalar(1), CV_FILLED);
			for (int i=min(posi,px);i<=max(posi,px);i++)
			{
				for (int j=0;j<tmpMat.rows;j++)
				{
					tmpMat(j,i) = tmpMat(j,i)*tagMat(j,i)>0 ? 1 : 0;
				}
			}

			Mat1b resMat;
			cv::warpAffine(tmpMat, resMat, getAffineTransform(pts2,pts1), tmpMat.size());
			vector<vector<cv::Point>> contours;
			findContours(resMat, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
			double maxB = 0;
			int posj = -1;
			for (int j=0; j<(int)contours.size(); j++)
			{
				RotatedRect rRec = minAreaRect(Mat(contours[j]));
				if (rRec.size.area()>=maxB)
				{
					maxB = rRec.size.area();
					posj = j;
				}
			}
			contour = contours[posj];
		}
		array<System::Drawing::Point>^ ContourForShow(vector<Point2i> &contour)
		{
			array<System::Drawing::Point>^ c = gcnew array<System::Drawing::Point>(contour.size());
			for (int i=0;i<(int)contour.size();i++)
			{
				c[i] = ToPoint(contour[i]);
			}
			return c;
		}
		bool FindHead(Mat1b& rData, vector<Point2i> &contour, cv::Rect &rect, Point2i &fp1, Point2i &fp2, Point2i &head)
		{

			int lx = rect.x, rx = lx + rect.width;
			int uy = rect.y, dy = uy + rect.height;

			double d = norm(fp1-fp2)/15;
			double area1 = 0, area2 = 0;
			for (int i=uy; i<dy; i++)
			{
				for (int j=lx; j<rx; j++)
				{
					if (rData(i,j)>0)
					{
						Point2d proj = Pnt2LineProj(fp1, fp2, Point2d(j,i));
						if (norm(proj-Point2d(fp1))<d)
						{
							area1++;
						} 
						if (norm(proj-Point2d(fp2))<d)
						{
							area2++;
						}
					}
				}
			}
			cv::Point2d hp = area1<area2 ? fp1 : fp2;

			Mat1b m0 = rData.clone();
			Mat1b m1 = Mat1b::zeros(rData.size());
			dilate(m0, m0, getStructuringElement(MORPH_RECT, cv::Size(15, 15)));
			erode (m0, m0, getStructuringElement(MORPH_RECT, cv::Size(15, 15)));
			Thin(Mat1b(m0,rect), Mat1b(m1,rect), 50);

			double params1[4], params2[4];
			GetLineParams(fp1, fp2, params1);
			GetVertParams(params1, (fp1+fp2)*0.5, params2);
			double mind = DBL_MAX, maxd = 0;
			Point2d kp1, kp2;
			for (int i=uy; i<dy; i++)
			{
				for (int j=lx; j<rx; j++)
				{
					if (m1(i,j)>0)
					{
						double d = Pnt2LineDist(params2, Point2d(j,i));
						if (d<=mind)
						{
							mind = d;
							kp1 = Point2d(j,i);
						}
						if (d>=maxd)
						{
							maxd = d;
							kp2 = Point2d(j,i);
						}
					}
				}
			}

			Point2d kp = norm(kp1-hp)<norm(kp2-hp) ? kp1 : kp2;
			double dist = norm(kp-hp);
			Point2d rp;
			for (int i=uy; i<dy; i++)
			{
				for (int j=lx; j<rx; j++)
				{
					if (m1(i,j)>0 && abs(norm(kp-Point2d(j,i))-dist)<=1)
					{
						rp = Point2d(j,i);
						break;
					}
				}
			}

			return Intersection(SplitContour(contour,params2,hp), kp, rp, head);
		}


		Point2i FindCenter(vector<Point2i> &contour, Point2i &p1, Point2i &p2)
		{
			int n = contour.size();
			double maxd = 0;
			int pi, pj;
			for (int i=0;i<n-1;i++)
			{
				for (int j=i+1;j<n;j++) 
				{
					double d = norm(contour[i]-contour[j]);
					if (d>=maxd)
					{
						maxd = d;
						pi = i;
						pj = j;
					}
				}
			}
			p1 = contour[pi];
			p2 = contour[pj];
			/////////////////////////////////////////////
			Point2d mp = (Point2d(p1)+Point2d(p2))*0.5;
			double params1[4], params2[4];
			GetLineParams(p1, p2, params1);
			GetVertParams(params1, mp, params2);
			/////////////////////////////////////////////
			int r = cvRound(maxd/4);
			Mat1d dMat = Mat1d::zeros(2*r+1, 2);
			Mat1i iMat(2*r+1, 2);
			for (int i=0;i<n;i++)
			{
				double d = Pnt2LineDist(params2, contour[i]);
				if (abs(d)<=r)
				{
					int pos = cvRound(d+r);
					d = Pnt2LineDist(params1, contour[i]);
					dMat(pos, d>0?0:1) = d;
					iMat(pos, d>0?0:1) = i;
				}
			}
			/////////////////////////////////////////////
			int gap = 4;
			double maxg = 0;
			int posg = r;
			for (int i=gap; i<dMat.rows-gap; i++)
			{
				int uc = 0, dc = 0;
				double us = 0, ds = 0;
				for (int j=-gap; j<=gap; j++)
				{
					if (dMat(i+j,0)>0)
					{
						us += dMat(i+j,0);
						uc ++;
					}
					if (dMat(i+j,1)<0)
					{
						ds += dMat(i+j,1);
						dc ++;
					}
				}
				if (uc*dc>0)
				{
					double g = us/uc - ds/dc;
					if (g>=maxg)
					{
						maxg = g;
						posg = i;
					}
				}
			}
			/////////////////////////////////////////////
			Point2i up = Point2i(), dp = Point2i();
			double uc = 0, dc = 0;
			for (int j=posg-gap; j<=posg+gap; j++)
			{
				if (dMat(j,0)>0)
				{
					up += contour[iMat(j,0)];
					uc ++;
				}
				if (dMat(j,1)<0)
				{
					dp += contour[iMat(j,1)];
					dc ++;
				}
			}

			if (uc*dc==0)
			{
				return Point2i(mp);
			}
			else
			{
				up *= (1/uc); 
				dp *= (1/dc);
				Point2d tmp = (up+dp)*0.5;
				return abs(tmp.x-mp.x)<abs(p2.x-p1.x)*0.5*midRangePct ? Point2i(tmp) : Point2i((tmp+mp)*0.5);
			}
		}
		array<System::Drawing::PointF>^ RectForShow(RotatedRect &rect)
		{
			Point2f pnts[4];
			rect.points(pnts);
			array<System::Drawing::PointF>^ ps = gcnew array<System::Drawing::PointF>(5);
			for (int i=0;i<4;i++)
			{
				ps[i] = ToPointF(pnts[i]);
			}
			ps[4] = ps[0];
			return ps;
		}
		double FindMoon(vector<Point2i> &contour, Point2i &head, Point2i &foot)
		{
			double retd = norm(head-foot)*moonRetPct;
			Point2d n = head - foot;
			Point2d n1 = Point2d(head) - (retd/norm(n))*n;
			Point2d n2 = Point2d(foot) + (retd/norm(n))*n;

			vector<Point2i> split1, split2;
			SplitContour(contour, head, foot, split1, split2);

			Point2i lp1, lp2, rp1, rp2;
			double params1[4], params2[4];
			GetLineParams(head, foot, params1);

			GetVertParams(params1, n1, params2);
			Intersection(split1, params2, lp1);
			Intersection(split2, params2, lp2);

			GetVertParams(params1, n2, params2);
			Intersection(split1, params2, rp1);
			Intersection(split2, params2, rp2);

			Point2d lmp = (lp1 + lp2)*0.5;
			Point2d rmp = (rp1 + rp2)*0.5;
			Point2d mmp = (lmp + rmp)*0.5;

			Point2i tp, dp;
			GetLineParams(lmp, rmp, params1);
			GetVertParams(params1, mmp, params2);
			Intersection(split1, params2, tp);
			Intersection(split2, params2, dp);

			double mvar = norm(Point2d(tp)-mmp)/norm(Point2d(dp)-mmp);
			return min(mvar, 1/mvar);
		}
		double Match(cv::Size &imgSize, vector<Point2i> &c1, vector<Point2i> &c2, vector<Point2i> &dst, Mat1b &dstMat)
		{
			int lx = cvRound(matchDeadX1);
			int rx = cvRound(matchDeadX2);

			vector<Point2i> tmp;
			for (vector<Point2i>::iterator i=c2.begin(); i!=c2.end(); i++)
			{
				if (i->x<lx || i->x>rx)
				{
					tmp.push_back(*i);
				}
			}

			int m = c1.size(), n = tmp.size();
			double ang = 0, tx = 0, ty = 0, epsilon = 0.0000001;
			Mat1d A(2*n,3), B(2*n,1), X(3,1);
			for (int i=0;i<20;i++)
			{
				double cosa = cos(ang);
				double sina = sin(ang);
				for (int j=0;j<n;j++)
				{
					double mind = DBL_MAX;
					for (int k=0;k<m;k++)
					{
						double dx = c1[k].x - ( cosa*tmp[j].x+sina*tmp[j].y+tx);
						double dy = c1[k].y - (-sina*tmp[j].x+cosa*tmp[j].y+ty);
						double d = dx*dx+dy*dy;
						if (d<=mind)
						{
							mind = d;
							B(j  ,0) = dx;
							B(j+n,0) = dy;
						}
					}
					A(j  ,0) =-sina*tmp[j].x+cosa*tmp[j].y;
					A(j  ,1) = 1;
					A(j  ,2) = 0;
					A(j+n,0) =-cosa*tmp[j].x-sina*tmp[j].y;
					A(j+n,1) = 0;
					A(j+n,2) = 1;
				}
				solve(A, B, X, DECOMP_SVD);
				ang+= X(0,0);
				tx += X(1,0);
				ty += X(2,0);
				if(abs(X(0,0))<epsilon)
				{
					break;
				}
			}

			double cosa = cos(-ang);
			double sina = sin(-ang);
			for (int i=0; i<m; i++)
			{
				double x = c1[i].x - tx;
				double y = c1[i].y - ty;
				c1[i] = Point2i(cvRound(cosa*x+sina*y), cvRound(-sina*x+cosa*y));
			}

			Mat1b m1 = Mat1b::zeros(imgSize);
			Mat1b m2 = Mat1b::zeros(imgSize);
			vector<vector<Point2i>> c1s(1); c1s[0] = c1;
			vector<vector<Point2i>> c2s(1); c2s[0] = c2;
			drawContours(m1, c1s, 0, Scalar(1), CV_FILLED);
			drawContours(m2, c2s, 0, Scalar(1), CV_FILLED);
			double ssl = sum(m1.col(lx))[0];
			double ssr = sum(m1.col(rx))[0];
			double sl = ssl==0 ? 1 : sum(m2.col(lx))[0]/ssl;
			double sr = ssr==0 ? 1 : sum(m2.col(rx))[0]/ssr;

			m1.colRange(lx,rx+1).copyTo(m2.colRange(lx,rx+1));
			vector<vector<cv::Point>> contours;
			findContours(m2.clone(), contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
			double maxA = 0;
			for (int i=0;i<(int)contours.size();i++)
			{
				RotatedRect rRect = minAreaRect(Mat(contours[i]));
				float h = min(rRect.size.width, rRect.size.height);
				float l = max(rRect.size.width, rRect.size.height);
				if (h*l>=maxA)
				{
					maxA = h*l;
					dst = contours[i];
				}
			}

			dstMat = Mat1b::zeros(imgSize);
			vector<vector<Point2i>> ccs(1); ccs[0] = dst;
			drawContours(dstMat, ccs, 0, Scalar(1), CV_FILLED);

			return min(sl, sr);
		}
		bool OneCut(vector<Point2i> &contour, Point2i &head, Point2i &foot, array<CutData^>^ %cutDatas)
		{
			double minh  = oneCutMinH/dimPerPix;
			double retd1 = oneCutHeadRet/dimPerPix;
			double retd2 = oneCutFootRet/dimPerPix;

			Point2d n = head - foot;
			Point2d n1 = Point2d(head) - (retd1/norm(n))*n;
			Point2d n2 = Point2d(foot) + (retd2/norm(n))*n;

			vector<Point2i> split1, split2;
			SplitContour(contour, head, foot, split1, split2);

			Point2i lp1, lp2, rp1, rp2;
			double params1[4], params2[4];
			GetLineParams(head, foot, params1);

			GetVertParams(params1, n1, params2);
			Intersection(split1, params2, lp1);
			Intersection(split2, params2, lp2);

			GetVertParams(params1, n2, params2);
			Intersection(split1, params2, rp1);
			Intersection(split2, params2, rp2);

			Point2d lmp = (lp1 + lp2)*0.5;
			Point2d rmp = (rp1 + rp2)*0.5;

			double h1 = MaxHeight(split1, lmp, rmp);
			double h2 = MaxHeight(split2, lmp, rmp);
			if (h1>=minh || h2>=minh)
			{
				GetVertParams(params1, (lmp+rmp)*0.5, params2);
				vector<Point2i> split3, split4;
				SplitContour(contour, params2, split3, split4);

				Point2d ns;
				double ls;
				if (h1<minh && h2>=minh)
				{
					ns = (rp2 - rp1)*0.5;
					ls = norm(ns);
				}
				else if (h1>=minh && h2<minh)
				{
					ns = (rp1 - rp2)*0.5;
					ls = norm(ns);
				}
				else
				{
					ns = Point2d();
					ls = 1;
				}
				for (int i=1; i<=cvRound(ls); i++)
				{
					Point2d rmps = rmp + (i/ls)*ns;
					vector<Point2i> split5, split6;
					SplitContour(contour, lmp, rmps, split5, split6);
					double h3 = MaxHeight(split5, lmp, rmps);
					double h4 = MaxHeight(split6, lmp, rmps);
					if (h3>=minh && h4>=minh)
					{
						Point2i bp1, bp2;
						Intersection(split3, lmp, rmps, bp1);
						Intersection(split4, lmp, rmps, bp2);
						double l = norm(bp1-bp2);

						cutDatas = gcnew array<CutData^>(2);
						cutDatas[0] = gcnew CutData(l*dimPerPix, h3*dimPerPix, 0, 0, lmp, rmps, 0);
						cutDatas[1] = gcnew CutData(l*dimPerPix, h4*dimPerPix, 0, 0, lmp, rmps, 0);
						AngDist(cutDatas[0]);
						AngDist(cutDatas[1]);
						return true;
					}
				}
			}
			return false;
		}
		bool TwoCut(vector<Point2i> &contour, Point2i &head, Point2i &foot, array<CutData^>^ %cutDatas)
		{
			double hot = twoCutHeadOfsThresh/dimPerPix;
			double ho1 = twoCutHeadOfs1/dimPerPix;
			double ho2 = twoCutHeadOfs2/dimPerPix;

			double retd1 = twoCutHeadRet/dimPerPix;
			Point2d n = head - foot;
			Point2d n1 = Point2d(head) - (retd1/norm(n))*n;
			Point2d n2 = Point2d(foot) + twoCutFootRet*n;

			vector<Point2i> split1, split2;
			SplitContour(contour, head, foot, split1, split2);

			double params1[4], params2[4];
			GetLineParams(head, foot, params1);

			Point2i hx1, hx2;
			GetVertParams(params1, n1, params2);
			Intersection(split1, params2, hx1);
			Intersection(split2, params2, hx2);

			Point2i fx1, fx2;
			GetVertParams(params1, n2, params2);
			Intersection(split1, params2, fx1);
			Intersection(split2, params2, fx2);

			Point2d hm = Point2d(hx1+hx2)*0.5;
			Point2d fm = Point2d(fx1+fx2)*0.5;
			SplitContour(contour, hm, fm, split1, split2);

			double hw = norm(hx1-hx2);
			double ho = hw<hot ? ho1 : ho2;
			Point2d hop1 = hm + Point2d(hx1-hx2)*(ho/hw)*0.5;
			Point2d hop2 = hm - Point2d(hx1-hx2)*(ho/hw)*0.5;

			double fw = norm(fx1-fx2);
			Point2d fop1 = fm + Point2d(fx1-fx2)*twoCutFootOfs*0.5;
			Point2d fop2 = fm - Point2d(fx1-fx2)*twoCutFootOfs*0.5; 
			//////////////////////////////////////////////////////////////////////////
			HeadOfsPnt1 = ToPointF(hop1);
			HeadCtrPnt1 = ToPointF(hx1);
			FootOfsPnt1 = ToPointF(fop1);
			FootCtrPnt1 = ToPointF(fx1);

			HeadOfsPnt2 = ToPointF(hop2);
			HeadCtrPnt2 = ToPointF(hx2);
			FootOfsPnt2 = ToPointF(fop2);
			FootCtrPnt2 = ToPointF(fx2);
			//////////////////////////////////////////////////////////////////////////
			vector<Point2i> c1 = SplitContour(split1, hop1, fop1, fx1);
			vector<Point2i> c2 = SplitContour(split2, hop2, fop2, fx2);
			vector<Point2i> ch = SplitContour(contour, hop1, hop2, head);
			CutData^ cd1 = CutByRot(MergeContour(ch, c1), hm, fm, hop1, lengthList1, lengthScoreList1, heightList1, heightScoreList1, 1);
			CutData^ cd2 = CutByRot(MergeContour(ch, c2), hm, fm, hop2, lengthList1, lengthScoreList1, heightList1, heightScoreList1, 1);

			CutData^ cd3 = gcnew CutData();
			CutData^ cd4 = gcnew CutData();
			if (twoCutFootStep==0)
			{
				cd3 = CutFree(c1, hm, fm, lengthList2, lengthScoreList2, heightList2, heightScoreList2, 2);
				cd4 = CutFree(c2, hm, fm, lengthList2, lengthScoreList2, heightList2, heightScoreList2, 2);
			}
			else if(twoCutFootStep>0)
			{
				Point2d stepn1 = Point2d(fx1) - fop1;
				int steps1 = cvRound(norm(stepn1)/twoCutFootStep);
				for (int i=0;i<steps1;i++)
				{
					Point2d rp = fop1 + (i*twoCutFootStep/norm(stepn1))*stepn1;
					CutData^ tmp = CutByRot(c1, hm, fm, rp, lengthList2, lengthScoreList2, heightList2, heightScoreList2, 2);
					if (tmp->SumScore()>=cd3->SumScore())
					{
						cd3 = tmp;
					}
				}
				Point2d stepn2 = Point2d(fx2) - fop2;
				int steps2 = cvRound(norm(stepn2)/twoCutFootStep);
				for (int i=0;i<steps2;i++)
				{
					Point2d rp = fop2 + (i*twoCutFootStep/norm(stepn2))*stepn2;
					CutData^ tmp = CutByRot(c2, hm, fm, rp, lengthList2, lengthScoreList2, heightList2, heightScoreList2, 2);
					if (tmp->SumScore()>=cd4->SumScore())
					{
						cd4 = tmp;
					}
				}
			}
			else
			{
				return false;
			}

			if (twoCutMode==1 && cd1->SumScore()>0 && cd2->SumScore()>0)
			{
				cutDatas = gcnew array<CutData^>(2);
				cutDatas[0] = cd1;
				cutDatas[1] = cd2;
				AngDist(cutDatas[0]);
				AngDist(cutDatas[1]);
				return true;
			}

			if (twoCutMode==2 || twoCutMode==3)
			{
				double s[4];
				s[0] = cd1->SumScore()>0 && cd2->SumScore()>0 ? cd1->SumScore()+cd2->SumScore() : 0;
				s[1] = cd2->SumScore()>0 && cd3->SumScore()>0 ? cd2->SumScore()+cd3->SumScore() : 0;
				s[2] = cd1->SumScore()>0 && cd4->SumScore()>0 ? cd1->SumScore()+cd4->SumScore() : 0;
				s[3] = cd3->SumScore()>0 && cd4->SumScore()>0 ? cd3->SumScore()+cd4->SumScore() : 0;
				if (twoCutMode==2)
				{
					s[3] = 0;
				}
				double maxs = 0;
				for (int i=0;i<4;i++)
				{
					if (s[i]>=maxs)
					{
						maxs = s[i];
					}
				}
				if (maxs>0)
				{
					cutDatas = gcnew array<CutData^>(2);
					if (maxs==s[0])
					{
						cutDatas[0] = cd1;
						cutDatas[1] = cd2;
					}
					if (maxs==s[1])
					{
						cutDatas[0] = cd2;
						cutDatas[1] = cd3;
					}
					if (maxs==s[2])
					{
						cutDatas[0] = cd1;
						cutDatas[1] = cd4;
					}
					if (maxs==s[3])
					{
						cutDatas[0] = cd3;
						cutDatas[1] = cd4;
					}
					AngDist(cutDatas[0]);
					AngDist(cutDatas[1]);
					return true;
				}
			}
			return false;
		}
		bool Check(Mat1b& grayMat)
		{
			Mat1b binMat, cutMat = grayMat(Range::all(), Range(cvRound(matchDeadX1), cvRound(matchDeadX2)));
			cv::threshold(cutMat, binMat, grayThresh, 255, THRESH_BINARY+THRESH_BINARY_INV);
			erode (binMat, binMat, getStructuringElement(MORPH_RECT, cv::Size(9, 9)));
			dilate(binMat, binMat, getStructuringElement(MORPH_RECT, cv::Size(9, 9)));

			vector<vector<cv::Point>> contours;
			findContours(binMat, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, Point2i(cvRound(matchDeadX1),0));

			cv::Point2f crcPnt1, crcPnt2, cutPnt1, cutPnt2;
			int count = 0;
			for (int i=0;i<(int)contours.size();i++)
			{
				cv::Rect rect = boundingRect(Mat(contours[i]));
				if (abs(rect.tl().x-matchDeadX1)>3 && abs(rect.br().x-matchDeadX2)>3)
				{
					RotatedRect rRct = fitEllipse(Mat(contours[i]));
					double rw = max(rRct.size.width, rRct.size.height);
					double rh = min(rRct.size.width, rRct.size.height);
					double mr = (rw+rh)/4;
					if (rw/rh<1.25 && mr>0.75*FlagCrcR && mr<1.25*FlagCrcR)
					{
						if (rRct.center.y<JigY)
						{
							crcPnt1 = rRct.center;
							cutPnt1 = crcPnt1 + cv::Point2f(0, (float)FlagSft1);
						}
						else
						{
							crcPnt2 = rRct.center;
							cutPnt2 = crcPnt2 - cv::Point2f(0, (float)FlagSft2);
						}
						count++;
					}
				}
			}

			if (count==2)
			{
				FlagCrcPnt1 = ToPointF(crcPnt1);
				FlagCrcPnt2 = ToPointF(crcPnt2);

				FlagCutPnt1 = ToPointF(cutPnt1);
				FlagCutPnt2 = ToPointF(cutPnt2);

				double params1[4], params2[4];
				GetLineParams(CutDatas[0]->CutPnt1, CutDatas[0]->CutPnt2, params1);
				GetLineParams(CutDatas[1]->CutPnt1, CutDatas[1]->CutPnt2, params2);

				return Pnt2LineDist(params1,cutPnt1)*Pnt2LineDist(params1,cutPnt2)<0 && Pnt2LineDist(params2,cutPnt1)*Pnt2LineDist(params2,cutPnt2)<0;
			}
			else
			{
				return false;
			}
		}
		////////////////////////////////////////////////////////////////
		template<typename T>
		System::Drawing::Point ToPoint(Point_<T> &p)
		{
			return System::Drawing::Point(cvRound(p.x), cvRound(p.y));
		}

		template<typename T>
		System::Drawing::PointF ToPointF(Point_<T> &p)
		{
			return System::Drawing::PointF(float(p.x), float(p.y));
		}

		template<typename T, typename T1, typename T2, typename T3> 
		CutData^ CutByRot(vector<Point_<T>> &contour, Point_<T1> &head, Point_<T2> &foot, Point_<T3> &rotPnt, array<double>^ lengthList, array<double>^ lengthScoreList, array<double>^ heightList, array<double>^ heightScoreList, int flag)
		{
			vector<Point_<T>> contour1, contour2;
			double params1[4], params2[4];
			GetLineParams(head, foot, params1);
			GetVertParams(params1, (Point2d(head)+Point2d(foot))*0.5, params2);
			SplitContour(contour, params2, contour1, contour2);

			double minl = twoCutMinLength/dimPerPix;
			double maxl = twoCutMaxLength/dimPerPix;
			double minh = twoCutMinHeight/dimPerPix;
			double maxh = twoCutMaxHeight/dimPerPix;
			double maxScore = 0, ml, mh, mls, mhs;
			int posi, posj;
			double params[4];
			for (int i=0;i<(int)contour1.size();i++)
			{
				GetLineParams(contour1[i], rotPnt, params);
				for (int j=0;j<(int)contour2.size();j++)
				{
					if(abs(Pnt2LineDist(params, contour2[j]))<1)
					{
						double l  = norm(contour1[i]-contour2[j]);
						if (l>minl && l<maxl)
						{
							double h  = MaxHeight(contour, contour1[i], contour2[j]);
							if (h>minh && h<maxh)
							{
								double ls = Score(lengthList, lengthScoreList, l*dimPerPix);
								double hs = Score(heightList, heightScoreList, h*dimPerPix);
								if (ls+hs>=maxScore)
								{
									maxScore = ls+hs;
									ml = l;
									mh = h;
									mls = ls;
									mhs = hs;
									posi = i;
									posj = j;
								}
							}
						}
					}
				}
			}
			return maxScore==0 ? gcnew CutData() : gcnew CutData(ml*dimPerPix, mh*dimPerPix, mls, mhs, contour1[posi], contour2[posj], flag);
		}

		template<typename T, typename T1, typename T2> 
		CutData^ CutFree(vector<Point_<T>> &contour, Point_<T1> &head, Point_<T2> &foot, array<double>^ lengthList, array<double>^ lengthScoreList, array<double>^ heightList, array<double>^ heightScoreList, int flag)
		{
			vector<Point_<T>> contour1, contour2;
			double params1[4], params2[4];
			GetLineParams(head, foot, params1);
			GetVertParams(params1, (Point2d(head)+Point2d(foot))*0.5, params2);
			SplitContour(contour, params2, contour1, contour2);

			double minl = twoCutMinLength/dimPerPix;
			double maxl = twoCutMaxLength/dimPerPix;
			double minh = twoCutMinHeight/dimPerPix;
			double maxh = twoCutMaxHeight/dimPerPix;
			double maxScore = 0, ml, mh, mls, mhs;
			int posi, posj;
			for (int i=0;i<(int)contour1.size();i++)
			{
				for (int j=0;j<(int)contour2.size();j++)
				{
					double l  = norm(contour1[i]-contour2[j]);
					if (l>minl && l<maxl)
					{
						double h  = MaxHeight(contour, contour1[i], contour2[j]);
						if (h>minh && h<maxh)
						{
							double ls = Score(lengthList, lengthScoreList, l*dimPerPix);
							double hs = Score(heightList, heightScoreList, h*dimPerPix);
							if (ls+hs>=maxScore)
							{
								maxScore = ls+hs;
								ml = l;
								mh = h;
								mls = ls;
								mhs = hs;
								posi = i;
								posj = j;
							}
						}
					}
				}
			}
			return maxScore==0 ? gcnew CutData() : gcnew CutData(ml*dimPerPix, mh*dimPerPix, mls, mhs, contour1[posi], contour2[posj], flag);
		}

		double Score(array<double>^ dimList, array<double>^ scoreList, double dim)
		{
			double score = 0;
			for (int i=0;i<dimList->Length-1;i++)
			{
				if (dim>=dimList[i] && dim<dimList[i+1])
				{
					score = scoreList[i];
				}
			}
			return score;
		}

		void Thin(Mat &srcImg, Mat &dstImg, int iterations)
		{
			IplImage* src = (IplImage*)&IplImage(srcImg);
			IplImage* dst = (IplImage*)&IplImage(dstImg);

			CvSize size = cvGetSize(src);
			cvCopy(src, dst);
			int n = 0,i = 0,j = 0;
			for(n=0; n<iterations; n++)
			{
				IplImage* t_image = cvCloneImage(dst);
				for(i=0; i<size.height;  i++)
				{
					for(j=0; j<size.width; j++)
					{
						if(CV_IMAGE_ELEM(t_image,byte,i,j)==1)
						{
							int ap=0;
							int p2 = (i==0)?0:CV_IMAGE_ELEM(t_image,byte, i-1, j);
							int p3 = (i==0 || j==size.width-1)?0:CV_IMAGE_ELEM(t_image,byte, i-1, j+1);
							if (p2==0 && p3==1)
							{
								ap++;
							}
							int p4 = (j==size.width-1)?0:CV_IMAGE_ELEM(t_image,byte,i,j+1);
							if(p3==0 && p4==1)
							{
								ap++;
							}
							int p5 = (i==size.height-1 || j==size.width-1)?0:CV_IMAGE_ELEM(t_image,byte,i+1,j+1);
							if(p4==0 && p5==1)
							{
								ap++;
							}
							int p6 = (i==size.height-1)?0:CV_IMAGE_ELEM(t_image,byte,i+1,j);
							if(p5==0 && p6==1)
							{
								ap++;
							}
							int p7 = (i==size.height-1 || j==0)?0:CV_IMAGE_ELEM(t_image,byte,i+1,j-1);
							if(p6==0 && p7==1)
							{
								ap++;
							}
							int p8 = (j==0)?0:CV_IMAGE_ELEM(t_image,byte,i,j-1);
							if(p7==0 && p8==1)
							{
								ap++;
							}
							int p9 = (i==0 || j==0)?0:CV_IMAGE_ELEM(t_image,byte,i-1,j-1);
							if(p8==0 && p9==1)
							{
								ap++;
							}
							if(p9==0 && p2==1)
							{
								ap++;
							}
							if((p2+p3+p4+p5+p6+p7+p8+p9)>1 && (p2+p3+p4+p5+p6+p7+p8+p9)<7)
							{
								if(ap==1)
								{
									if(!(p2 && p4 && p6))
									{
										if(!(p4 && p6 && p8)) 
										{
											CV_IMAGE_ELEM(dst,byte,i,j)=0;
										}
									}
								}
							}

						}
					}
				}
				cvReleaseImage(&t_image);
				t_image = cvCloneImage(dst);
				for(i=0; i<size.height;  i++)
				{
					for(int j=0; j<size.width; j++)
					{
						if(CV_IMAGE_ELEM(t_image,byte,i,j)==1)
						{
							int ap=0;
							int p2 = (i==0)?0:CV_IMAGE_ELEM(t_image,byte, i-1, j);
							int p3 = (i==0 || j==size.width-1)?0:CV_IMAGE_ELEM(t_image,byte, i-1, j+1);
							if (p2==0 && p3==1)
							{
								ap++;
							}
							int p4 = (j==size.width-1)?0:CV_IMAGE_ELEM(t_image,byte,i,j+1);
							if(p3==0 && p4==1)
							{
								ap++;
							}
							int p5 = (i==size.height-1 || j==size.width-1)?0:CV_IMAGE_ELEM(t_image,byte,i+1,j+1);
							if(p4==0 && p5==1)
							{
								ap++;
							}
							int p6 = (i==size.height-1)?0:CV_IMAGE_ELEM(t_image,byte,i+1,j);
							if(p5==0 && p6==1)
							{
								ap++;
							}
							int p7 = (i==size.height-1 || j==0)?0:CV_IMAGE_ELEM(t_image,byte,i+1,j-1);
							if(p6==0 && p7==1)
							{
								ap++;
							}
							int p8 = (j==0)?0:CV_IMAGE_ELEM(t_image,byte,i,j-1);
							if(p7==0 && p8==1)
							{
								ap++;
							}
							int p9 = (i==0 || j==0)?0:CV_IMAGE_ELEM(t_image,byte,i-1,j-1);
							if(p8==0 && p9==1)
							{
								ap++;
							}
							if(p9==0 && p2==1)
							{
								ap++;
							}
							if((p2+p3+p4+p5+p6+p7+p8+p9)>1 && (p2+p3+p4+p5+p6+p7+p8+p9)<7)
							{
								if(ap==1)
								{
									if(p2*p4*p8==0)
									{
										if(p2*p6*p8==0)
										{
											CV_IMAGE_ELEM(dst, byte,i,j)=0;
										}
									}
								}
							}                    
						}

					}

				}            
				cvReleaseImage(&t_image);
			}
		}

		void TransParam(Vec4f &cvline, double params[])
		{
			params[0] =-cvline[1];
			params[1] = cvline[0];
			params[2] =-params[0]*cvline[2] - params[1]*cvline[3];
			params[3] = sqrt(params[0]*params[0] + params[1]*params[1]);
		}

		template<typename T> 
		vector<Point_<T>> MergeContour(vector<Point_<T>> &split1, vector<Point_<T>> &split2)
		{
			vector<Point_<T>> contour;
			for (int i=0;i<(int)split1.size();i++)
			{
				contour.push_back(split1[i]);
			}
			for (int i=0;i<(int)split2.size();i++)
			{
				contour.push_back(split2[i]);
			}
			return contour;
		}

		template<typename T1, typename T2, typename T3, typename T4> 
		vector<Point_<T1>> SplitContour(vector<Point_<T1>> &contour, Point_<T2> &p1, Point_<T3> &p2, Point_<T4> &kPnt)
		{
			double params[4];
			GetLineParams(p1, p2, params);
			return SplitContour(contour, params, kPnt);
		}

		template<typename T1, typename T2> 
		vector<Point_<T1>> SplitContour(vector<Point_<T1>> &contour, double params[], Point_<T2> &kPnt)
		{
			vector<Point_<T1>> split;
			double flag = Pnt2LineDist(params, kPnt);
			for (vector<Point_<T1>>::iterator i=contour.begin(); i!=contour.end(); i++)
			{
				if (Pnt2LineDist(params, *i)*flag>0)
				{
					split.push_back(*i);
				}
			}
			return split;
		}

		template<typename T1, typename T2, typename T3> 
		void SplitContour(vector<Point_<T1>> &contour, Point_<T2> &p1, Point_<T3> &p2, vector<Point_<T1>> &split1, vector<Point_<T1>> &split2)
		{
			double params[4];
			GetLineParams(p1, p2, params);
			SplitContour(contour, params, split1, split2);
		}

		template<typename T> 
		void SplitContour(vector<Point_<T>> &contour, double params[], vector<Point_<T>> &split1, vector<Point_<T>> &split2)
		{
			split1.clear();
			split2.clear();
			for (vector<Point_<T>>::iterator i=contour.begin(); i!=contour.end(); i++)
			{
				double d = Pnt2LineDist(params, *i);
				if (d>0)
				{
					split1.push_back(*i);
				}
				if (d<0)
				{
					split2.push_back(*i);
				}
			}
		}

		template<typename T1, typename T2, typename T3> 
		double MaxHeightPnt(vector<Point_<T1>> &contour, Point_<T2> &p1, Point_<T3> &p2, Point_<T1> &mp)
		{
			double md = 0;
			for (vector<Point_<T1>>::iterator i=contour.begin(); i!=contour.end(); i++)
			{
				Point2d proj = Pnt2LineProj(p1, p2, Point2d(*i));
				Point2d n1 = Point2d(p1) - proj;
				Point2d n2 = Point2d(p2) - proj;
				if (n1.dot(n2)<0)
				{
					double dx = proj.x - i->x;
					double dy = proj.y - i->y;
					double d = dx*dx + dy*dy;
					if (d>=md)
					{
						md = d;
						mp = *i;
					}
				}
			}
			return sqrt(md);
		}

		template<typename T1, typename T2, typename T3> 
		double MaxHeight(vector<Point_<T1>> &contour, Point_<T2> &p1, Point_<T3> &p2)
		{
			double md = 0;
			for (vector<Point_<T1>>::iterator i=contour.begin(); i!=contour.end(); i++)
			{
				Point2d proj = Pnt2LineProj(p1, p2, Point2d(*i));
				Point2d n1 = Point2d(p1) - proj;
				Point2d n2 = Point2d(p2) - proj;
				if (n1.dot(n2)<0)
				{
					double dx = proj.x - i->x;
					double dy = proj.y - i->y;
					double d = dx*dx + dy*dy;
					if (d>=md)
					{
						md = d;
					}
				}
			}
			return sqrt(md);
		}

		template<typename T1, typename T2, typename T3> 
		bool Intersection(vector<Point_<T3>> &contour, Point_<T1> &p1, Point_<T2> &p2, Point_<T3> &x)
		{
			double params[4];
			GetLineParams(p1, p2, params);
			return Intersection(contour, params, x);
		}

		template<typename T> 
		bool Intersection(vector<Point_<T>> &contour, double params[], Point_<T> &x)
		{
			double md = DBL_MAX;
			for (vector<Point_<T>>::iterator i=contour.begin(); i!=contour.end(); i++)
			{
				double d = abs(Pnt2LineDist(params, *i));
				if (d<=md)
				{
					md = d;
					x = *i;
				}
			}
			return md>1 ? false : true;
		}

		template<typename T> 
		void GetVertParams(double iParams[], Point_<T> &pnt, double oParams[])
		{
			oParams[0] = iParams[1];
			oParams[1] =-iParams[0];
			oParams[2] =-oParams[0]*pnt.x-oParams[1]*pnt.y;
			oParams[3] = iParams[3];
		}

		template<typename T1, typename T2>
		void GetLineParams(Point_<T1> &p1, Point_<T2> &p2, double params[])
		{
			double x1 = p1.x;
			double y1 = p1.y;
			double x2 = p2.x;
			double y2 = p2.y;
			params[0] = y1 - y2;
			params[1] = x2 - x1;
			params[2] = -x1*params[0]-y1*params[1];
			params[3] = sqrt(params[0]*params[0] + params[1]*params[1]);
		}

		void GetLineParams(System::Drawing::PointF p1, System::Drawing::PointF p2, double params[])
		{
			GetLineParams(Point2f(p1.X, p1.Y), Point2f(p2.X, p2.Y), params);
		}

		template<typename T> 
		double Pnt2LineDist(double params[], Point_<T> &pnt)
		{
			double d = params[0]*pnt.x + params[1]*pnt.y + params[2];
			return d/params[3];
		}

		template<typename T1, typename T2, typename T3> 
		Point_<T3> Pnt2LineProj(Point_<T1> &P1, Point_<T2> &P2, Point_<T3> &P3)
		{
			double a1 = P2.x - P1.x;
			double b1 = P2.y - P1.y;
			double a1a1 = a1*a1;
			double b1b1 = b1*b1;
			double denominator = a1a1 + b1b1;
			if (denominator == 0) 
			{
				return P3;
			}
			double x1y2 = P1.x*P2.y;
			double x2y1 = P2.x*P1.y;
			double a1b1 = a1*b1;
			double moleculey = b1b1*P3.y + a1b1*P3.x - a1*x1y2 + a1*x2y1;
			double moleculex = a1a1*P3.x + a1b1*P3.y - b1*x2y1 + b1*x1y2;
			Point2d p = Point2d(moleculex/denominator, moleculey/denominator);
			return Point_<T3>(p);
		}

		void AngDist(CutData^ %cutData)
		{
			double dy = cutData->CutPnt1.Y - cutData->CutPnt2.Y ;
			double dx = cutData->CutPnt2.X - cutData->CutPnt1.X;
			cutData->Ang = atan(dy/dx)*180/M_PI;
			double x1 = (cutData->CutPnt1.X - JigX)*dimPerPix;
			double y1 =-(cutData->CutPnt1.Y - JigY)*dimPerPix;
			double x2 = (cutData->CutPnt2.X - JigX)*dimPerPix;
			double y2 =-(cutData->CutPnt2.Y - JigY)*dimPerPix;
			double a = y1 - y2;
			double b = x2 - x1;
			double c =-x1*a - y1*b;
			cutData->Dist = CutY + (a*CutX+c)/b;
		}
	};

	public ref class Calibration
	{
	public:
		static void Detect(IntPtr imgPtr, int imgW, int imgH, int cornerRow, int cornerCol)
		{
			Mat img(imgH, imgW, CV_8UC(3), (byte*)imgPtr.ToPointer(), (size_t)ceil(imgW*3/4.0)*4);
			vector<Point2f> corner;
			cv::Size cornerSize(cornerCol, cornerRow);
			bool found = false;
			int winSize = 11;
			if (findChessboardCorners(img, cornerSize, corner, CV_CALIB_CB_FAST_CHECK + CV_CALIB_CB_ADAPTIVE_THRESH + CV_CALIB_CB_FILTER_QUADS) || 
				findChessboardCorners(img, cornerSize, corner, CV_CALIB_CB_FAST_CHECK + CV_CALIB_CB_ADAPTIVE_THRESH + CALIB_CB_NORMALIZE_IMAGE))
			{
				found = true;
			}
			/////////////////////////////////////
			cv::line(img, Point2i(0, imgH/2), Point2i(imgW-1, imgH/2), Scalar(0,255,0), 3);
			cv::line(img, Point2i(imgW/2, 0), Point2i(imgW/2, imgH-1), Scalar(0,255,0), 3);
			/////////////////////////////////////
			if (found)
			{
				int k = 3;
				for (int i=0;i<cornerRow;i++)
				{
					for (int j=0;j<cornerCol;j++)
					{
						Point2i mp = Point2i(corner[i*cornerCol+j]);
						Point2i p1 = mp + Point2i(-k,0);
						Point2i p2 = mp + Point2i(+k,0);
						Point2i p3 = mp + Point2i(0,-k);
						Point2i p4 = mp + Point2i(0,+k);
						if (i==cornerRow/2 && j==cornerCol/2)
						{
							cv::line(img, p1, p2, Scalar(255,0,0), 3);
							cv::line(img, p3, p4, Scalar(255,0,0), 3);
						}
						else
						{
							cv::line(img, p1, p2, Scalar(0,0,255), 3);
							cv::line(img, p3, p4, Scalar(0,0,255), 3);
						}
					}
				}
			}
		}
		static bool Detect(IntPtr imgPtr, int imgW, int imgH, int cornerRow, int cornerCol, double dimCell, [Out] double% dimPerPix)
		{
			Mat img(imgH, imgW, CV_8UC(3), (byte*)imgPtr.ToPointer(), (size_t)ceil(imgW*3/4.0)*4);
			vector<Point2f> corner;
			cv::Size cornerSize(cornerCol, cornerRow);
			bool found = false;
			int winSize = 11;
			if (findChessboardCorners(img, cornerSize, corner, CV_CALIB_CB_FAST_CHECK + CV_CALIB_CB_ADAPTIVE_THRESH + CV_CALIB_CB_FILTER_QUADS) || 
				findChessboardCorners(img, cornerSize, corner, CV_CALIB_CB_FAST_CHECK + CV_CALIB_CB_ADAPTIVE_THRESH + CALIB_CB_NORMALIZE_IMAGE))
			{
				//cornerSubPix(img, corner, cv::Size(winSize,winSize), cv::Size(-1,-1), TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 20, DBL_EPSILON));
				double d = norm(corner[0]-corner[cornerCol-1]) + norm(corner[(cornerRow-1)*cornerCol]-corner[cornerRow*cornerCol-1]);
				dimPerPix = 2*dimCell*(cornerCol-1)/d;
				found = true;
			}
			/////////////////////////////////////
			cv::line(img, Point2i(0, imgH/2), Point2i(imgW-1, imgH/2), Scalar(0,255,0), 3);
			cv::line(img, Point2i(imgW/2, 0), Point2i(imgW/2, imgH-1), Scalar(0,255,0), 3);
			/////////////////////////////////////
			if (found)
			{
				int k = 3;
				for (int i=0;i<cornerRow;i++)
				{
					for (int j=0;j<cornerCol;j++)
					{
						Point2i mp = Point2i(corner[i*cornerCol+j]);
						Point2i p1 = mp + Point2i(-k,0);
						Point2i p2 = mp + Point2i(+k,0);
						Point2i p3 = mp + Point2i(0,-k);
						Point2i p4 = mp + Point2i(0,+k);
						if (i==cornerRow/2 && j==cornerCol/2)
						{
							cv::line(img, p1, p2, Scalar(255,0,0), 3);
							cv::line(img, p3, p4, Scalar(255,0,0), 3);
						}
						else
						{
							cv::line(img, p1, p2, Scalar(0,0,255), 3);
							cv::line(img, p3, p4, Scalar(0,0,255), 3);
						}
					}
				}
			}
			return found;
		}
	};
} 

