// CalibrateView.cpp : CCalibrateView 类的实现
//

#include "stdafx.h"
#include "Calibrate.h"

#include<stdlib.h>
#include <stdio.h>
#include "CalibrateDoc.h"
#include "CalibrateView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


//灰度阈值
void /*__stdcall*/ ThreshGray_callback(int Pos, void* usrdata)    //回调函数1
{
	Mat binarizedImage;
	Mat dst;
	binarizedImage = *(Mat*)(usrdata);
	Canny(binarizedImage,dst,Pos,Pos*2,3);
	imshow("SimpleBlobDetector",dst);
}

void /*__stdcall*/ ThreshEdge_callback(int Pos, void* usrdata)    //
{
	Mat binarizedImage;
	Mat dst;
	binarizedImage = *(Mat*)(usrdata);
	threshold(binarizedImage,dst,Pos,255,THRESH_BINARY);
	imshow("SimpleBlobDetector",dst);
}

// CCalibrateView

IMPLEMENT_DYNCREATE(CCalibrateView, CFormView)

BEGIN_MESSAGE_MAP(CCalibrateView, CFormView)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_BN_CLICKED(IDC_BTN_DCali, &CCalibrateView::OnBnClickedBtnDcali)
	ON_BN_CLICKED(IDC_BTN_GETP, &CCalibrateView::OnBnClickedBtnGetp)
END_MESSAGE_MAP()

// CCalibrateView 构造/析构

CCalibrateView::CCalibrateView()
	: CFormView(CCalibrateView::IDD)
{
	// TODO: 在此处添加构造代码
	TestNum = 14;							   //测试图像对数

	m_LMarkers = 30;						   //标定板上的相邻标记点间距离

	m_iImgCounts = 14;						   //标定图像个数
	m_iRowCount=8;							   //标定板行的标记点个数
	m_iColCount=8;							   //标定板列的标记点个数
	m_iPointCounts = m_iRowCount*m_iColCount;  //标定板上标记点的个数

	m_iLeftThreshold = 23;					   //左图像阈值
	m_iRightThreshold = 25;					   //右图像阈值
}

CCalibrateView::~CCalibrateView()
{
}

void CCalibrateView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_PROGRESS1, m_CProgress);
	DDX_Control(pDX, IDC_LIST_SHOWEROR, m_CShowError);
}

BOOL CCalibrateView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: 在此处通过修改
	//  CREATESTRUCT cs 来修改窗口类或样式

	return CFormView::PreCreateWindow(cs);
}

void CCalibrateView::OnInitialUpdate()
{
	CFormView::OnInitialUpdate();
	GetParentFrame()->RecalcLayout();
	ResizeParentToFit();


	//1. 初始化表格控件
	CRect rect;
	m_CShowError.GetWindowRect(&rect);
	m_CShowError.InsertColumn(0,_T("序号"),LVCFMT_LEFT,rect.Width()*1/4,-1);
	m_CShowError.InsertColumn(1,_T("标定板的编号"),LVCFMT_LEFT,rect.Width()*1/4,-1);
	m_CShowError.InsertColumn(2,_T("标定板的平均误差"),LVCFMT_LEFT,rect.Width()*1/4,-1);
	m_CShowError.InsertColumn(3,_T("总相对误差"),LVCFMT_LEFT,rect.Width()*1/4,-1);
	m_CShowError.SetExtendedStyle(m_CShowError.GetExtendedStyle()|LVS_EX_FULLROWSELECT | LVS_EX_GRIDLINES|LVS_EX_SUBITEMIMAGES );//此函数使得点击的某行全亮

	//2.初始化实验测量路径
	memcpy( m_aLeftReadDirectory,"G:\\Calibration\\ChessboardImg\\LeftCamera\\",MAX_PATH );
	memcpy( m_aLeftSaveDirectory,"G:\\Calibration\\Detecte\\Left\\",MAX_PATH );
	memcpy( m_aRightReadDirectory,"G:\\Calibration\\ChessboardImg\\RightCamera\\",MAX_PATH );
	memcpy( m_aRightSaveDirectory,"G:\\Calibration\\Detecte\\Right\\",MAX_PATH );
	memcpy( m_aTestReadDirectory,"G:\\Calibration\\Test\\Read\\",MAX_PATH );
	memcpy( m_aTestSaveDirectory,"G:\\Calibration\\Test\\Save\\",MAX_PATH );

	//3.初始化opencv的矩阵、图像，主要是分配空间
	InitCameraParams();

}

void CCalibrateView::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CCalibrateView::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CCalibrateView 诊断

#ifdef _DEBUG
void CCalibrateView::AssertValid() const
{
	CFormView::AssertValid();
}

void CCalibrateView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}

CCalibrateDoc* CCalibrateView::GetDocument() const // 非调试版本是内联的
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CCalibrateDoc)));
	return (CCalibrateDoc*)m_pDocument;
}
#endif //_DEBUG


// CCalibrateView 消息处理程序

/***************************************************************************
* 函数名称：   CCalibrationView::InitCameraParams
* 摘    要：   初始化相机内外参数
* 全局影响：   public 
* 参    数：   int m_iImgCounts
* 参    数：   int m_iPointCounts
* 参    数：   int m_iChessboardCount
* 返回值：     void
* 
* 修改记录： 7.28
*  [日期]     [作者/修改者]  [修改原因]
*2015/07/25       罗佳丽          添加
***************************************************************************/
void CCalibrateView::InitCameraParams()
{
	//初始化多块内参
	m_pcvBig_dimeObj_Left = cvCreateMat(CHESSNUM * m_iPointCounts * m_iImgCounts , 3 ,CV_32FC1);
	m_pcvBig_dimeObj_Right = cvCreateMat(CHESSNUM * m_iPointCounts * m_iImgCounts , 3 ,CV_32FC1);
	m_pcvImgPoint_Left = cvCreateMat(CHESSNUM * m_iPointCounts * m_iImgCounts, 2, CV_32FC1);
	m_pcvImgPoint_Right = cvCreateMat(CHESSNUM * m_iPointCounts * m_iImgCounts, 2, CV_32FC1);

	// 摄像机内参数
	m_pcvIntrinsicMat_Left = cvCreateMat(3, 3, CV_64FC1);
	m_pcvIntrinsicMat_Right = cvCreateMat(3, 3, CV_64FC1);
	m_pcvDistortionCoeffs_Left = cvCreateMat(5, 1, CV_64FC1);
	m_pcvDistortionCoeffs_Right = cvCreateMat(5, 1, CV_64FC1);

	// 摄像机外参数
	m_pcvRotationVct_Left = cvCreateMat(m_iImgCounts, 3, CV_64FC1);
	m_pcvRotationVct_Right = cvCreateMat(m_iImgCounts, 3, CV_64FC1);
	m_pcvTranslationVct_Left = cvCreateMat(m_iImgCounts, 3, CV_64FC1);
	m_pcvTranslationVct_Right = cvCreateMat(m_iImgCounts, 3, CV_64FC1);

	//每张照片中，标记点的个数
	m_pcvPointCounts = cvCreateMat(m_iImgCounts, 1, CV_32SC1); 

	m_pcvTestImg_Left = cvCreateMat( CHESSNUM * m_iPointCounts,2,CV_32FC1);
	m_pcvTestImg_Right = cvCreateMat( CHESSNUM * m_iPointCounts,2,CV_32FC1);

	m_pcvImgPoint_LeftReProjectiom = cvCreateMat(CHESSNUM * m_iImgCounts * m_iPointCounts,2, CV_32FC1);
	m_pcvImgPoint_RightReProjectiom = cvCreateMat(CHESSNUM * m_iImgCounts * m_iPointCounts,2, CV_32FC1);
}

int CCalibrateView::dbCmp( const void* _a, const void* _b )
{
	double a = *(const double*)_a;
	double b = *(const double*)_b;

	return (a > b) - (a < b);
}

/***************************************************************************
* 函数名称：   CCalibrationView::OnBnClickedBtnDcali
* 摘    要：   双目标定
* 全局影响：   public 
* 返回值：     void
* 
* 修改记录：
*  [日期]     [作者/修改者]  [修改原因]
*2016/08/09       熊磊         添加
***************************************************************************/
void CCalibrateView::OnBnClickedBtnDcali()
{
	//进度条
	m_CProgress.SetRange(0,2*m_iImgCounts+2);

	//处理左相机图片
	

	ProcessPicture(m_aLeftReadDirectory,m_aLeftSaveDirectory,m_iLeftThreshold,m_iImgCounts,m_pcvBig_dimeObj_Left,m_pcvImgPoint_Left,m_pcvPointCounts,m_imageSize_Left,
		m_aRightReadDirectory,m_aRightSaveDirectory,m_iRightThreshold,m_pcvBig_dimeObj_Right,m_pcvImgPoint_Right,m_pcvPointCounts,m_imageSize_Right);

		CvMat* FundamentalMat = cvCreateMat(3, 3, CV_64FC1);		// 基础矩阵,即learning OpenCv 425页的F矩阵
		CvMat* EigenValueMat = cvCreateMat(3, 3, CV_64FC1);			// 本征矩阵，即learning OpenCv 425页的E矩阵

		CvMat* RotatonStereo = cvCreateMat(3, 3, CV_64FC1);			// 右目相机坐标系转换到左目相机坐标系所需的旋转矩阵，即learning OpenCv 428页的R矩阵
		CvMat* TranslationStereo = cvCreateMat(3, 1, CV_64FC1);		// 右目相机坐标系转换到左目相机坐标系所需的平移矩阵，即learning OpenCv 428页的T矩阵
		CvMat* ProjectionMat_Left = cvCreateMat(3, 4, CV_64FC1);    // 左目相机投影方程，即3*4的矩阵Pl
		CvMat* ProjectionMat_Right = cvCreateMat(3, 4, CV_64FC1);	// 右目相机投影方程，即3*4的矩阵Pr
		CvMat* Rl = cvCreateMat(3, 3, CV_64FC1);					//左目相机行对准矩阵
		CvMat* Rr = cvCreateMat(3, 3, CV_64FC1);					//右目相机行对准矩阵
		m_pcvProjectionMat = cvCreateMat(4, 4, CV_64FC1);			// 重投影矩阵，即4*4的矩阵Q


		m_pcvMapX_Left = cvCreateMat( m_imageSize_Left.height,m_imageSize_Left.width, CV_32FC1);
		m_pcvMapY_Left = cvCreateMat( m_imageSize_Left.height,m_imageSize_Left.width, CV_32FC1);
		m_pcvMapX_Right = cvCreateMat( m_imageSize_Right.height,m_imageSize_Right.width, CV_32FC1);
		m_pcvMapY_Right = cvCreateMat( m_imageSize_Right.height,m_imageSize_Right.width, CV_32FC1); 

		cvStereoCalibrate1(m_pcvBig_dimeObj_Left,m_pcvBig_dimeObj_Right,
			m_pcvImgPoint_Left,m_pcvImgPoint_Right,m_pcvPointCounts,
			m_pcvIntrinsicMat_Left,m_pcvDistortionCoeffs_Left,
			m_pcvIntrinsicMat_Right,m_pcvDistortionCoeffs_Right,
			m_imageSize_Left,RotatonStereo,TranslationStereo,
			EigenValueMat,FundamentalMat,
			cvTermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS,100,1e-5),CV_CALIB_FIX_INTRINSIC);//CV_CALIB_USE_INTRINSIC_GUESS);///立体标定


		cvStereoRectify(m_pcvIntrinsicMat_Left,m_pcvIntrinsicMat_Right,m_pcvDistortionCoeffs_Left,m_pcvDistortionCoeffs_Right,
			m_imageSize_Left,RotatonStereo,TranslationStereo,Rl,Rr,ProjectionMat_Left,
			ProjectionMat_Right,m_pcvProjectionMat,CV_CALIB_ZERO_DISPARITY);
		cvInitUndistortRectifyMap(m_pcvIntrinsicMat_Left, m_pcvDistortionCoeffs_Left, Rl, 
			ProjectionMat_Left, m_pcvMapX_Left, m_pcvMapY_Left);
		cvInitUndistortRectifyMap(m_pcvIntrinsicMat_Right, m_pcvDistortionCoeffs_Right, Rr,
			ProjectionMat_Right, m_pcvMapX_Right, m_pcvMapY_Right);

		MessageBox(_T("双目标定完成"));
}


double CCalibrateView::cvStereoCalibrate1( const CvMat* _objectPoints1,const CvMat* _objectPoints2, 
						const CvMat* _imagePoints1,const CvMat* _imagePoints2, const CvMat* _npoints,
                        CvMat* _cameraMatrix1, CvMat* _distCoeffs1,
                        CvMat* _cameraMatrix2, CvMat* _distCoeffs2,
                        CvSize imageSize, CvMat* matR, CvMat* matT,
                        CvMat* matE, CvMat* matF,
                        CvTermCriteria termCrit,
                        int flags )
{
    const int NINTRINSIC = 12;
    Ptr<CvMat> npoints, err, J_LR, Je, Ji, imagePoints[2], objectPoints[2], RT0;
    CvLevMarq solver;
    double reprojErr = 0;

    double A[2][9], dk[2][8]={{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}}, rlr[9];
    CvMat K[2], Dist[2], om_LR, T_LR;
    CvMat R_LR = cvMat(3, 3, CV_64F, rlr);
    int i, k, p, ni = 0, ofs, nimages, pointsTotal, maxPoints = 0;
    int nparams;
    bool recomputeIntrinsics = false;
    double aspectRatio[2] = {0,0};

    CV_Assert( CV_IS_MAT(_imagePoints1) && CV_IS_MAT(_imagePoints2) &&
               CV_IS_MAT(_objectPoints1) && CV_IS_MAT(_npoints) &&
               CV_IS_MAT(matR) && CV_IS_MAT(matT) );

    CV_Assert( CV_ARE_TYPES_EQ(_imagePoints1, _imagePoints2) &&
               CV_ARE_DEPTHS_EQ(_imagePoints1, _objectPoints1) );

    CV_Assert( (_npoints->cols == 1 || _npoints->rows == 1) &&
               CV_MAT_TYPE(_npoints->type) == CV_32SC1 );

    nimages = _npoints->cols + _npoints->rows - 1;
    npoints = cvCreateMat( _npoints->rows, _npoints->cols, _npoints->type );
    cvCopy( _npoints, npoints );

    for( i = 0, pointsTotal = 0; i < nimages; i++ )
    {
        maxPoints = MAX(maxPoints, npoints->data.i[i]);
        pointsTotal += npoints->data.i[i];
    }

    for( k = 0; k < 2; k++ )
    {
        const CvMat* points = k == 0 ? _imagePoints1 : _imagePoints2;
		const CvMat* objpoints = k == 0 ? _objectPoints1 : _objectPoints2;
        const CvMat* cameraMatrix = k == 0 ? _cameraMatrix1 : _cameraMatrix2;
        const CvMat* distCoeffs = k == 0 ? _distCoeffs1 : _distCoeffs2;

        int cn = CV_MAT_CN(_imagePoints1->type);
        CV_Assert( (CV_MAT_DEPTH(_imagePoints1->type) == CV_32F ||
                CV_MAT_DEPTH(_imagePoints1->type) == CV_64F) &&
               ((_imagePoints1->rows == pointsTotal && _imagePoints1->cols*cn == 2) ||
                (_imagePoints1->rows == 1 && _imagePoints1->cols == pointsTotal && cn == 2)) );

        K[k] = cvMat(3,3,CV_64F,A[k]);
        Dist[k] = cvMat(1,8,CV_64F,dk[k]);

		objectPoints[k] = cvCreateMat(objpoints->rows,objpoints->cols,CV_64FC(CV_MAT_CN(objpoints->type)));
        cvConvert(objpoints,objectPoints[k]);
		cvReshape( objectPoints[k], objectPoints[k], 3, 1 );
        imagePoints[k] = cvCreateMat( points->rows, points->cols, CV_64FC(CV_MAT_CN(points->type)));
        cvConvert( points, imagePoints[k] );
        cvReshape( imagePoints[k], imagePoints[k], 2, 1 );

        if( flags & (CV_CALIB_FIX_INTRINSIC|CV_CALIB_USE_INTRINSIC_GUESS|
            CV_CALIB_FIX_ASPECT_RATIO|CV_CALIB_FIX_FOCAL_LENGTH) )
            cvConvert( cameraMatrix, &K[k] );

        if( flags & (CV_CALIB_FIX_INTRINSIC|CV_CALIB_USE_INTRINSIC_GUESS|
            CV_CALIB_FIX_K1|CV_CALIB_FIX_K2|CV_CALIB_FIX_K3|CV_CALIB_FIX_K4|CV_CALIB_FIX_K5|CV_CALIB_FIX_K6) )
        {
            CvMat tdist = cvMat( distCoeffs->rows, distCoeffs->cols,
                CV_MAKETYPE(CV_64F,CV_MAT_CN(distCoeffs->type)), Dist[k].data.db );
            cvConvert( distCoeffs, &tdist );
        }

        if( !(flags & (CV_CALIB_FIX_INTRINSIC|CV_CALIB_USE_INTRINSIC_GUESS)))
        {
            cvCalibrateCamera2( objectPoints[k], imagePoints[k],
                npoints, imageSize, &K[k], &Dist[k], 0, 0, flags );
        }
    }

    if( flags & CV_CALIB_SAME_FOCAL_LENGTH )
    {
        static const int avg_idx[] = { 0, 4, 2, 5, -1 };
        for( k = 0; avg_idx[k] >= 0; k++ )
            A[0][avg_idx[k]] = A[1][avg_idx[k]] = (A[0][avg_idx[k]] + A[1][avg_idx[k]])*0.5;
    }

    if( flags & CV_CALIB_FIX_ASPECT_RATIO )
    {
        for( k = 0; k < 2; k++ )
            aspectRatio[k] = A[k][0]/A[k][4];
    }

    recomputeIntrinsics = (flags & CV_CALIB_FIX_INTRINSIC) == 0;

    err = cvCreateMat( maxPoints*2, 1, CV_64F );
    Je = cvCreateMat( maxPoints*2, 6, CV_64F );
    J_LR = cvCreateMat( maxPoints*2, 6, CV_64F );
    Ji = cvCreateMat( maxPoints*2, NINTRINSIC, CV_64F );
    cvZero( Ji );

    // we optimize for the inter-camera R(3),t(3), then, optionally,
    // for intrinisic parameters of each camera ((fx,fy,cx,cy,k1,k2,p1,p2) ~ 8 parameters).
    nparams = 6*(nimages+1) + (recomputeIntrinsics ? NINTRINSIC*2 : 0);

    // storage for initial [om(R){i}|t{i}] (in order to compute the median for each component)
    RT0 = cvCreateMat( 6, nimages, CV_64F );

    solver.init( nparams, 0, termCrit );
    if( recomputeIntrinsics )
    {
        uchar* imask = solver.mask->data.ptr + nparams - NINTRINSIC*2;
        if( !(flags & CV_CALIB_RATIONAL_MODEL) )
            flags |= CV_CALIB_FIX_K4 | CV_CALIB_FIX_K5 | CV_CALIB_FIX_K6;
        if( flags & CV_CALIB_FIX_ASPECT_RATIO )
            imask[0] = imask[NINTRINSIC] = 0;
        if( flags & CV_CALIB_FIX_FOCAL_LENGTH )
            imask[0] = imask[1] = imask[NINTRINSIC] = imask[NINTRINSIC+1] = 0;
        if( flags & CV_CALIB_FIX_PRINCIPAL_POINT )
            imask[2] = imask[3] = imask[NINTRINSIC+2] = imask[NINTRINSIC+3] = 0;
        if( flags & CV_CALIB_ZERO_TANGENT_DIST )
            imask[6] = imask[7] = imask[NINTRINSIC+6] = imask[NINTRINSIC+7] = 0;
        if( flags & CV_CALIB_FIX_K1 )
            imask[4] = imask[NINTRINSIC+4] = 0;
        if( flags & CV_CALIB_FIX_K2 )
            imask[5] = imask[NINTRINSIC+5] = 0;
        if( flags & CV_CALIB_FIX_K3 )
            imask[8] = imask[NINTRINSIC+8] = 0;
        if( flags & CV_CALIB_FIX_K4 )
            imask[9] = imask[NINTRINSIC+9] = 0;
        if( flags & CV_CALIB_FIX_K5 )
            imask[10] = imask[NINTRINSIC+10] = 0;
        if( flags & CV_CALIB_FIX_K6 )
            imask[11] = imask[NINTRINSIC+11] = 0;
    }

    for( i = ofs = 0; i < nimages; ofs += ni, i++ )
    {
        ni = npoints->data.i[i];
        CvMat objpt_i[2];
        double _om[2][3], r[2][9], t[2][3];
        CvMat om[2], R[2], T[2], imgpt_i[2];

        
        for( k = 0; k < 2; k++ )
        {
			objpt_i[k] = cvMat(1, ni, CV_64FC3, objectPoints[k]->data.db + ofs*3);
            imgpt_i[k] = cvMat(1, ni, CV_64FC2, imagePoints[k]->data.db + ofs*2);
            om[k] = cvMat(3, 1, CV_64F, _om[k]);
            R[k] = cvMat(3, 3, CV_64F, r[k]);
            T[k] = cvMat(3, 1, CV_64F, t[k]);

            // FIXME: here we ignore activePoints[k] because of
            // the limited API of cvFindExtrnisicCameraParams2
            cvFindExtrinsicCameraParams2( &objpt_i[k], &imgpt_i[k], &K[k], &Dist[k], &om[k], &T[k] );
            cvRodrigues2( &om[k], &R[k] );
            if( k == 0 )
            {
                // save initial om_left and T_left
                solver.param->data.db[(i+1)*6] = _om[0][0];
                solver.param->data.db[(i+1)*6 + 1] = _om[0][1];
                solver.param->data.db[(i+1)*6 + 2] = _om[0][2];
                solver.param->data.db[(i+1)*6 + 3] = t[0][0];
                solver.param->data.db[(i+1)*6 + 4] = t[0][1];
                solver.param->data.db[(i+1)*6 + 5] = t[0][2];
            }
        }
        cvGEMM( &R[1], &R[0], 1, 0, 0, &R[0], CV_GEMM_B_T );
        cvGEMM( &R[0], &T[0], -1, &T[1], 1, &T[1] );
        cvRodrigues2( &R[0], &T[0] );
        RT0->data.db[i] = t[0][0];
        RT0->data.db[i + nimages] = t[0][1];
        RT0->data.db[i + nimages*2] = t[0][2];
        RT0->data.db[i + nimages*3] = t[1][0];
        RT0->data.db[i + nimages*4] = t[1][1];
        RT0->data.db[i + nimages*5] = t[1][2];
    }

    // find the medians and save the first 6 parameters
    for( i = 0; i < 6; i++ )
    {
        qsort( RT0->data.db + i*nimages, nimages, CV_ELEM_SIZE(RT0->type), dbCmp );
        solver.param->data.db[i] = nimages % 2 != 0 ? RT0->data.db[i*nimages + nimages/2] :
            (RT0->data.db[i*nimages + nimages/2 - 1] + RT0->data.db[i*nimages + nimages/2])*0.5;
    }

    if( recomputeIntrinsics )
        for( k = 0; k < 2; k++ )
        {
            double* iparam = solver.param->data.db + (nimages+1)*6 + k*NINTRINSIC;
            if( flags & CV_CALIB_ZERO_TANGENT_DIST )
                dk[k][2] = dk[k][3] = 0;
            iparam[0] = A[k][0]; iparam[1] = A[k][4]; iparam[2] = A[k][2]; iparam[3] = A[k][5];
            iparam[4] = dk[k][0]; iparam[5] = dk[k][1]; iparam[6] = dk[k][2];
            iparam[7] = dk[k][3]; iparam[8] = dk[k][4]; iparam[9] = dk[k][5];
            iparam[10] = dk[k][6]; iparam[11] = dk[k][7];
        }

    om_LR = cvMat(3, 1, CV_64F, solver.param->data.db);
    T_LR = cvMat(3, 1, CV_64F, solver.param->data.db + 3);

    for(;;)
    {
        const CvMat* param = 0;
        CvMat tmpimagePoints;
        CvMat *JtJ = 0, *JtErr = 0;
        double *_errNorm = 0;
        double _omR[3], _tR[3];
        double _dr3dr1[9], _dr3dr2[9], /*_dt3dr1[9],*/ _dt3dr2[9], _dt3dt1[9], _dt3dt2[9];
        CvMat dr3dr1 = cvMat(3, 3, CV_64F, _dr3dr1);
        CvMat dr3dr2 = cvMat(3, 3, CV_64F, _dr3dr2);
        //CvMat dt3dr1 = cvMat(3, 3, CV_64F, _dt3dr1);
        CvMat dt3dr2 = cvMat(3, 3, CV_64F, _dt3dr2);
        CvMat dt3dt1 = cvMat(3, 3, CV_64F, _dt3dt1);
        CvMat dt3dt2 = cvMat(3, 3, CV_64F, _dt3dt2);
        CvMat om[2], T[2], imgpt_i[2],objpt_i[2];
        CvMat dpdrot_hdr, dpdt_hdr, dpdf_hdr, dpdc_hdr, dpdk_hdr;
        CvMat *dpdrot = &dpdrot_hdr, *dpdt = &dpdt_hdr, *dpdf = 0, *dpdc = 0, *dpdk = 0;

        if( !solver.updateAlt( param, JtJ, JtErr, _errNorm ))
            break;
        reprojErr = 0;

        cvRodrigues2( &om_LR, &R_LR );
        om[1] = cvMat(3,1,CV_64F,_omR);
        T[1] = cvMat(3,1,CV_64F,_tR);

        if( recomputeIntrinsics )
        {
            double* iparam = solver.param->data.db + (nimages+1)*6;
            double* ipparam = solver.prevParam->data.db + (nimages+1)*6;
            dpdf = &dpdf_hdr;
            dpdc = &dpdc_hdr;
            dpdk = &dpdk_hdr;
            if( flags & CV_CALIB_SAME_FOCAL_LENGTH )
            {
                iparam[NINTRINSIC] = iparam[0];
                iparam[NINTRINSIC+1] = iparam[1];
                ipparam[NINTRINSIC] = ipparam[0];
                ipparam[NINTRINSIC+1] = ipparam[1];
            }
            if( flags & CV_CALIB_FIX_ASPECT_RATIO )
            {
                iparam[0] = iparam[1]*aspectRatio[0];
                iparam[NINTRINSIC] = iparam[NINTRINSIC+1]*aspectRatio[1];
                ipparam[0] = ipparam[1]*aspectRatio[0];
                ipparam[NINTRINSIC] = ipparam[NINTRINSIC+1]*aspectRatio[1];
            }
            for( k = 0; k < 2; k++ )
            {
                A[k][0] = iparam[k*NINTRINSIC+0];
                A[k][4] = iparam[k*NINTRINSIC+1];
                A[k][2] = iparam[k*NINTRINSIC+2];
                A[k][5] = iparam[k*NINTRINSIC+3];
                dk[k][0] = iparam[k*NINTRINSIC+4];
                dk[k][1] = iparam[k*NINTRINSIC+5];
                dk[k][2] = iparam[k*NINTRINSIC+6];
                dk[k][3] = iparam[k*NINTRINSIC+7];
                dk[k][4] = iparam[k*NINTRINSIC+8];
                dk[k][5] = iparam[k*NINTRINSIC+9];
                dk[k][6] = iparam[k*NINTRINSIC+10];
                dk[k][7] = iparam[k*NINTRINSIC+11];
            }
        }

        for( i = ofs = 0; i < nimages; ofs += ni, i++ )
        {
            ni = npoints->data.i[i];
            CvMat  _part;

            om[0] = cvMat(3,1,CV_64F,solver.param->data.db+(i+1)*6);
            T[0] = cvMat(3,1,CV_64F,solver.param->data.db+(i+1)*6+3);

            if( JtJ || JtErr )
                cvComposeRT( &om[0], &T[0], &om_LR, &T_LR, &om[1], &T[1], &dr3dr1, 0,
                             &dr3dr2, 0, 0, &dt3dt1, &dt3dr2, &dt3dt2 );
            else
                cvComposeRT( &om[0], &T[0], &om_LR, &T_LR, &om[1], &T[1] );

            
            err->rows = Je->rows = J_LR->rows = Ji->rows = ni*2;
            cvReshape( err, &tmpimagePoints, 2, 1 );

            cvGetCols( Ji, &dpdf_hdr, 0, 2 );
            cvGetCols( Ji, &dpdc_hdr, 2, 4 );
            cvGetCols( Ji, &dpdk_hdr, 4, NINTRINSIC );
            cvGetCols( Je, &dpdrot_hdr, 0, 3 );
            cvGetCols( Je, &dpdt_hdr, 3, 6 );

            for( k = 0; k < 2; k++ )
            {

                double maxErr, l2err;
                imgpt_i[k] = cvMat(1, ni, CV_64FC2, imagePoints[k]->data.db + ofs*2);
				objpt_i[k] = cvMat(1, ni, CV_64FC3, objectPoints[k]->data.db + ofs*3);

                if( JtJ || JtErr )
                    cvProjectPoints2( &objpt_i[k], &om[k], &T[k], &K[k], &Dist[k],
                            &tmpimagePoints, dpdrot, dpdt, dpdf, dpdc, dpdk,
                            (flags & CV_CALIB_FIX_ASPECT_RATIO) ? aspectRatio[k] : 0);
                else
                    cvProjectPoints2( &objpt_i[k], &om[k], &T[k], &K[k], &Dist[k], &tmpimagePoints );
                cvSub( &tmpimagePoints, &imgpt_i[k], &tmpimagePoints );

                l2err = cvNorm( &tmpimagePoints, 0, CV_L2 );
                maxErr = cvNorm( &tmpimagePoints, 0, CV_C );

                if( JtJ || JtErr )
                {
                    int iofs = (nimages+1)*6 + k*NINTRINSIC, eofs = (i+1)*6;
                    assert( JtJ && JtErr );

                    if( k == 1 )
                    {
                        // d(err_{x|y}R) ~ de3
                        // convert de3/{dr3,dt3} => de3{dr1,dt1} & de3{dr2,dt2}
                        for( p = 0; p < ni*2; p++ )
                        {
                            CvMat de3dr3 = cvMat( 1, 3, CV_64F, Je->data.ptr + Je->step*p );
                            CvMat de3dt3 = cvMat( 1, 3, CV_64F, de3dr3.data.db + 3 );
                            CvMat de3dr2 = cvMat( 1, 3, CV_64F, J_LR->data.ptr + J_LR->step*p );
                            CvMat de3dt2 = cvMat( 1, 3, CV_64F, de3dr2.data.db + 3 );
                            double _de3dr1[3], _de3dt1[3];
                            CvMat de3dr1 = cvMat( 1, 3, CV_64F, _de3dr1 );
                            CvMat de3dt1 = cvMat( 1, 3, CV_64F, _de3dt1 );

                            cvMatMul( &de3dr3, &dr3dr1, &de3dr1 );
                            cvMatMul( &de3dt3, &dt3dt1, &de3dt1 );

                            cvMatMul( &de3dr3, &dr3dr2, &de3dr2 );
                            cvMatMulAdd( &de3dt3, &dt3dr2, &de3dr2, &de3dr2 );

                            cvMatMul( &de3dt3, &dt3dt2, &de3dt2 );

                            cvCopy( &de3dr1, &de3dr3 );
                            cvCopy( &de3dt1, &de3dt3 );
                        }

                        cvGetSubRect( JtJ, &_part, cvRect(0, 0, 6, 6) );
                        cvGEMM( J_LR, J_LR, 1, &_part, 1, &_part, CV_GEMM_A_T );

                        cvGetSubRect( JtJ, &_part, cvRect(eofs, 0, 6, 6) );
                        cvGEMM( J_LR, Je, 1, 0, 0, &_part, CV_GEMM_A_T );

                        cvGetRows( JtErr, &_part, 0, 6 );
                        cvGEMM( J_LR, err, 1, &_part, 1, &_part, CV_GEMM_A_T );
                    }

                    cvGetSubRect( JtJ, &_part, cvRect(eofs, eofs, 6, 6) );
                    cvGEMM( Je, Je, 1, &_part, 1, &_part, CV_GEMM_A_T );

                    cvGetRows( JtErr, &_part, eofs, eofs + 6 );
                    cvGEMM( Je, err, 1, &_part, 1, &_part, CV_GEMM_A_T );

                    if( recomputeIntrinsics )
                    {
                        cvGetSubRect( JtJ, &_part, cvRect(iofs, iofs, NINTRINSIC, NINTRINSIC) );
                        cvGEMM( Ji, Ji, 1, &_part, 1, &_part, CV_GEMM_A_T );
                        cvGetSubRect( JtJ, &_part, cvRect(iofs, eofs, NINTRINSIC, 6) );
                        cvGEMM( Je, Ji, 1, &_part, 1, &_part, CV_GEMM_A_T );
                        if( k == 1 )
                        {
                            cvGetSubRect( JtJ, &_part, cvRect(iofs, 0, NINTRINSIC, 6) );
                            cvGEMM( J_LR, Ji, 1, &_part, 1, &_part, CV_GEMM_A_T );
                        }
                        cvGetRows( JtErr, &_part, iofs, iofs + NINTRINSIC );
                        cvGEMM( Ji, err, 1, &_part, 1, &_part, CV_GEMM_A_T );
                    }
                }

                reprojErr += l2err*l2err;
            }
        }
        if(_errNorm)
            *_errNorm = reprojErr;
    }

    cvRodrigues2( &om_LR, &R_LR );
    if( matR->rows == 1 || matR->cols == 1 )
        cvConvert( &om_LR, matR );
    else
        cvConvert( &R_LR, matR );
    cvConvert( &T_LR, matT );

    if( recomputeIntrinsics )
    {
        cvConvert( &K[0], _cameraMatrix1 );
        cvConvert( &K[1], _cameraMatrix2 );

        for( k = 0; k < 2; k++ )
        {
            CvMat* distCoeffs = k == 0 ? _distCoeffs1 : _distCoeffs2;
            CvMat tdist = cvMat( distCoeffs->rows, distCoeffs->cols,
                CV_MAKETYPE(CV_64F,CV_MAT_CN(distCoeffs->type)), Dist[k].data.db );
            cvConvert( &tdist, distCoeffs );
        }
    }

    if( matE || matF )
    {
        double* t = T_LR.data.db;
        double tx[] =
        {
            0, -t[2], t[1],
            t[2], 0, -t[0],
            -t[1], t[0], 0
        };
        CvMat Tx = cvMat(3, 3, CV_64F, tx);
        double e[9], f[9];
        CvMat E = cvMat(3, 3, CV_64F, e);
        CvMat F = cvMat(3, 3, CV_64F, f);
        cvMatMul( &Tx, &R_LR, &E );
        if( matE )
            cvConvert( &E, matE );
        if( matF )
        {
            double ik[9];
            CvMat iK = cvMat(3, 3, CV_64F, ik);
            cvInvert(&K[1], &iK);
            cvGEMM( &iK, &E, 1, 0, 0, &E, CV_GEMM_A_T );
            cvInvert(&K[0], &iK);
            cvMatMul(&E, &iK, &F);
            cvConvertScale( &F, matF, fabs(f[8]) > 0 ? 1./f[8] : 1 );
        }
    }

    return std::sqrt(reprojErr/(pointsTotal*2));
}

/***************************************************************************
* 函数名称：   CCalibrationView::ProcessPicture
* 摘    要：   对图片处理，得到单目标定需要的参数///前两个参数为输入，后面为输出
* 全局影响：   public 
* 参    数：   char picaddress[]
* 参    数：   int pictcount
* 参    数：   CvMat * objpoints
* 参    数：   CvMat * imgpoints
* 参    数：   CvMat * everypicpointcnt
* 参    数：   CvSize imgsize
* 返回值：     BOOL
* 
* 修改记录：
*  [日期]     [作者/修改者]  [修改原因]
*2016/03/13       韩杨杨         添加
***************************************************************************/
BOOL CCalibrateView::ProcessPicture(char picreadaddress_L[],char picsaveaddr_L[],int ithreshold_L, int pictcount,CvMat * objpoints_L, CvMat * imgpoints_L,CvMat * everypicpointcnt_L,CvSize & imgsize_L,char picreadaddress_R[],char picsaveaddr_R[],int ithreshold_R, CvMat * objpoints_R, CvMat * imgpoints_R,CvMat * everypicpointcnt_R,CvSize & imgsize_R )
{
	//0.定义函数内部的变量
	CvMat * _ObjPoints_L,*_ObjPoints_R;
	_ObjPoints_L = cvCreateMat(objpoints_L->rows,objpoints_L->cols,CV_32FC1);
	_ObjPoints_R = cvCreateMat(objpoints_R->rows,objpoints_R->cols,CV_32FC1);
	
	CvMat * _ImgPoints_L,* _ImgPoints_R;
    _ImgPoints_L = cvCreateMat(imgpoints_L->rows,imgpoints_L->cols,CV_32FC1);
	_ImgPoints_R = cvCreateMat(imgpoints_R->rows,imgpoints_R->cols,CV_32FC1);

	CvMat * _PointsCnt_L, *_PointsCnt_R;
	_PointsCnt_L = cvCreateMat(everypicpointcnt_L->rows,everypicpointcnt_L->cols,CV_32SC1);
	_PointsCnt_R = cvCreateMat(everypicpointcnt_R->rows,everypicpointcnt_R->cols,CV_32SC1);
	
	CvSize _ImgSize_L,_ImgSize_R;


	vector<KeyPoint> key_points_L,key_points_R;  
	vector<ST_CIRCLE_POINT>Circle_Ponits_L,Circle_Ponits_R;

	//存放排序好的点
	VEC_POINTS Computer_points_L(CHESSNUM),Computer_points_R(CHESSNUM);

	// 棋盘格原点坐标
	PST_UpLeft_POINTS pstUpLeftPoints_L,pstUpLeftPoints_R;

	//棋盘格标记点
	PST_ORGIN_POINTS pstOrginPoints_L,pstOrginPoints_R;


	//形成每个小靶标的像素坐标
	CvMat *pcvImg_L[CHESSNUM];
	CvMat *pcvImg_R[CHESSNUM];


	//定义文件读取、保存路径
	char FileNameTail[10];



	CvMat* m_pcvH[CHESSNUM][CAMERANUM];                     //记录各个小靶标的透视矩阵
	CvMat* m_pcvRV[CHESSNUM][CAMERANUM];					 //记录各个小靶标转换到靶标1的转换矩阵


	//每幅图的点数，主要是用于调用函数cvCalibrateCamera2_H
	CvMat * iPointCont1_L;
	CvMat * iPointCont1_R;
	iPointCont1_L  = cvCreateMat(1,1,CV_32SC1);
	iPointCont1_R  = cvCreateMat(1,1,CV_32SC1);

	//保存每幅图片的标记点个数（64*2）
	int * piPointNum_L;
	int * piPointNum_R;
	piPointNum_L = new int[pictcount];
	piPointNum_R = new int[pictcount];


	//形成每个小靶标的物理坐标
	float * obj_L = NULL;
	float * obj_R = NULL;

	obj_R = new float[m_iPointCounts*3];
	obj_L = new float[m_iPointCounts*3];

	CvMat * pcvObj_R;
	pcvObj_R = cvCreateMat(m_iPointCounts,1,CV_32FC3);
	(*pcvObj_R)=cvMat(m_iPointCounts,1,CV_32FC3,obj_R);
	obj_L = new float[m_iPointCounts*3];


	CvMat * pcvObj_L;
	pcvObj_L = cvCreateMat(m_iPointCounts,1,CV_32FC3);
	(*pcvObj_L)=cvMat(m_iPointCounts,1,CV_32FC3,obj_L);

	for (int m = 0;m < m_iRowCount;m++)
	{
		for (int n = 0;n < m_iColCount;n++)
		{
			obj_L[ m*m_iRowCount*3 + (n*3)+0]= m * m_LMarkers + 10;
			obj_L[ m*m_iRowCount*3 + (n*3)+1] = n * m_LMarkers + 10;
			obj_L[ m*m_iRowCount*3 + (n*3)+2]= 0;
		}
	}


	for (int m = 0;m < m_iRowCount;m++)
	{
		for (int n = 0;n < m_iColCount;n++)
		{
			obj_R[ m*m_iRowCount*3 + (n*3)+0]= m * m_LMarkers + 10;
			obj_R[ m*m_iRowCount*3 + (n*3)+1] = n * m_LMarkers + 10;
			obj_R[ m*m_iRowCount*3 + (n*3)+2]= 0;
		}
	}

	IplImage* img_R = NULL;
	
for (int n = 0;n < CAMERANUM;n++)
{
	for (int imgcnt=0;imgcnt<pictcount;imgcnt++)
	{
		//进度条
		m_CProgress.SetStep(1);
		m_CProgress.StepIt();

			if (n ==0)
			{
				sprintf(FileNameTail,"_LB");

				for (int i=0;i<CHESSNUM;i++)
				{
					m_pcvH[i][0] = cvCreateMat(3,3,CV_64FC1);
					m_pcvRV[i][0] = cvCreateMat(3,3,CV_64FC1);
					pcvImg_L[i] = cvCreateMat(m_iPointCounts,CHESSNUM,CV_32FC1);
				}


				sprintf(m_FileName,"%sChessboardImg1%s (%d).bmp",picreadaddress_L,FileNameTail,imgcnt + 1);	

				//彩色图片
				IplImage* img1=cvLoadImage(m_FileName,-1);

				IplImage* img_L= cvCreateImage(cvGetSize(img1),8,1);
				cvCvtColor(img1,img_L,CV_BGR2GRAY);

				cvReleaseImage(&img1);
				img1 = NULL;

				pTemp = cvCreateImage(cvGetSize(img_L),IPL_DEPTH_8U,1);
				uchar*ptr=(uchar*)img_L->imageData;


				int height=img_L->height;
				int width=img_L->width;
				int pixpbyte=img_L->nChannels;

				_ImgSize_L = cvSize(width,height);


				int iGrayVar = 0;
				if (ithreshold_L <= 50)
				{
					iGrayVar = ithreshold_L;
				}
				else if (ithreshold_L <= 100)
				{
					iGrayVar = ithreshold_L*2/3;
				}
				else  
				{
					iGrayVar = ithreshold_L/2;
				}

				int iMemThres = ithreshold_L-iGrayVar;
				if (iMemThres <= 0)
				{
					iMemThres = 0;
				}

				FindCircleCenter(img_L,iMemThres,iMemThres+iGrayVar,(iGrayVar/10),50,10000,0.05f/*凸度*/,0.6f/*惯性率*/,0.6f/*圆度*/,key_points_L,m_iPointCounts,CHESSNUM,0);
				cvCanny(img_L, pTemp, (iMemThres+iGrayVar/2),(iMemThres+iGrayVar/2)*2, 3 );


				sprintf(m_FileName,"");
				sprintf(m_FileName,"%scvCannyImg%d%s.bmp",picsaveaddr_L,imgcnt+1,FileNameTail);
				cvSaveImage(m_FileName,pTemp);

				CalCircleCentrePoint(pTemp,key_points_L,m_iPointCounts,CHESSNUM);
				cvReleaseImage(&pTemp);

				//将圆心画在圆上
				for (int i=0;i<key_points_L.size();i++)
				{
					int x=(int) key_points_L.at(i).pt.x;
					int y=(int) key_points_L.at(i).pt.y;
					*(ptr+y*width+x)=(uchar)255;
				}

				
				sprintf(m_FileName,"");
				sprintf(m_FileName,"%sDetectedImg%d%s.bmp",picsaveaddr_L,imgcnt+1,FileNameTail);
				cvSaveImage(m_FileName,img_L);

				cout<<"记录所有找到的标记点坐标，并存为文件"<<endl;
				//记录下第一块板上的所有坐标
				FILE *pf=NULL;
				sprintf(m_FileName,"");
				sprintf(m_FileName,"%sgridcorner[%d].txt",picsaveaddr_L,imgcnt+1);//其中0：左相机；1:右相机
				pf=fopen(m_FileName,"w");
				int countNumber=0;

				for (int i=0;i<key_points_L.size();i++)
				{
					ST_CIRCLE_POINT CirclePoints;
					CirclePoints.x=key_points_L.at(i).pt.x;
					CirclePoints.y=key_points_L.at(i).pt.y;
					CirclePoints.size=key_points_L.at(i).size;
					CirclePoints.distance=0;
					countNumber++;
					Circle_Ponits_L.push_back(CirclePoints);
					fprintf(pf,"%.3f %.3f %d",key_points_L.at(i).pt.x,key_points_L.at(i).pt.y,countNumber);
					fprintf(pf,"\r\n");
				}

				fclose(pf);

				pf=NULL;


				//获得标定板的标记点坐标
				GetOrignPoints(Circle_Ponits_L,CHESSNUM,pstOrginPoints_L);

				//获得标定板的原点坐标
				GetUpLeftPoints(Circle_Ponits_L,pstOrginPoints_L,pstUpLeftPoints_L,CHESSNUM);

				//圆心排序
				SortPoints(Computer_points_L,Circle_Ponits_L,pstUpLeftPoints_L,CHESSNUM,m_iPointCounts);

				//将顺序画在图上
				for (int k=0;k<CHESSNUM;k++)
				{
					if ((int)Computer_points_L[k].size()!=m_iPointCounts)
					{
						CString str;
						str.Format(_T("第%d张标定图片不能找到角点"),imgcnt+1);
						MessageBox(str);

					}
					for (int j=0;j<(int)Computer_points_L[k].size()-1;j++)
					{

						float x= Computer_points_L[k].at(j).x;
						float y= Computer_points_L[k].at(j).y;

						float x1= Computer_points_L[k].at(j+1).x;
						float y1= Computer_points_L[k].at(j+1).y;

						CvPoint pt1;
						CvPoint pt2;
						pt1.x=x;
						pt1.y=y;
						pt2.x=x1;
						pt2.y=y1;

						cvLine(img_L,pt1,pt2,cvScalar(255,0,0),1,8,0);

					}
				}


				sprintf_s(m_FileName,"");
				sprintf_s(m_FileName,"%sDetectedSortedImg%d%s.bmp",picsaveaddr_L,imgcnt+1,FileNameTail);
				cvSaveImage(m_FileName,img_L);

				//将排列好的标记点赋值给矩阵，

				float *pfdata=NULL;
				float *pfdata1 = NULL;
				//pfdata=m_pcvImgPoints_Left->data.fl;
				int istep=0;
				int istep1 =0;

				int cont =0;
				for (int chessnumber=0;chessnumber<CHESSNUM;chessnumber++)
				{   

					//保存排序好的标记点的像素坐标
					sprintf_s(m_FileName,"");
					sprintf_s(m_FileName,"%sImgPoints_[%d]_%d.txt",picsaveaddr_L,imgcnt+1,chessnumber);///////
					pf=fopen(m_FileName,"w");

					pfdata=_ImgPoints_L->data.fl;
					istep=_ImgPoints_L->step/sizeof(float);

					pfdata1 = pcvImg_L[chessnumber]->data.fl;       //用来求得每个图片上的小靶标的透视矩阵
					istep1 = pcvImg_L[chessnumber]->step/sizeof(float);

					for (int i=0;i< Computer_points_L[chessnumber].size();i++)
					{
						pfdata1[i*istep1] = Computer_points_L[chessnumber].at(i).x;
						pfdata1[i*istep1+1] = Computer_points_L[chessnumber].at(i).y;


						//给矩阵赋值
						cont = imgcnt*CHESSNUM*m_iPointCounts*istep + chessnumber*m_iPointCounts*istep + i*istep;
						pfdata[cont] = Computer_points_L[chessnumber].at(i).x;
						pfdata[cont + 1] = Computer_points_L[chessnumber].at(i).y;



						fprintf(pf,"%.3f %.3f %d ",pfdata[cont],pfdata[cont + 1],i);
						fprintf(pf,"\r\n");

					}
					pfdata=NULL;
					pfdata1 = NULL;

					fclose(pf);
					pf=NULL;
					sprintf(m_FileName,"");

				}

				piPointNum_L[imgcnt]=m_iPointCounts*CHESSNUM;
				cvSetData(iPointCont1_L,&m_iPointCounts,sizeof(int));

				Circle_Ponits_L.clear();
				for (int chessnumber=0;chessnumber<CHESSNUM;chessnumber++)
				{
					Computer_points_L[chessnumber].clear();

				}

				//0对每幅图片的多块标定板分别获得各自的透视矩阵H[i]
				for (int i=0;i<CHESSNUM;i++)
				{

					cvCalibrateCamera2_H(pcvObj_L,pcvImg_L[i],iPointCont1_L,_ImgSize_L,m_pcvH[i][0],i,0);

				}

				//1 求得各个矩阵到H1的转换矩阵
				CvMat* matH0_Invert_L = cvCreateMat(3, 3, CV_64FC1);
				cvZero(matH0_Invert_L);
				cvInvert(m_pcvH[0][0], matH0_Invert_L,CV_LU);//改动

				float * pData = NULL;    //接收m_pcvObjectPoints[j]
				int step = 0;
				float * pData1 = NULL;   //m_pcvThree_dimeObj
				int step11 = 0;
				double * pD_RV = NULL;   //接收m_pcvRV[j]
				int step_RV = 0;

				double* pD_R = NULL;   //临时输出H0的逆
				int step_R = 0;

				pD_R= matH0_Invert_L->data.db;
				step_R = matH0_Invert_L->step/sizeof(double);


				sprintf(m_FileName,"%sH0_INV[%d]",picsaveaddr_L,imgcnt+1);   //输出每幅图片中靶标1的H1的逆
				pf =fopen(m_FileName,"w");

				for (int j = 0;j<3;j++)
				{
					fprintf(pf,"%.3f %.3f %.3f " , pD_R[j*step_R+0],pD_R[j*step_R+1],pD_R[j*step_R+2]);
					fprintf(pf,"\r\n");
				} 
				fclose(pf);
				pf=NULL;
				pD_R =NULL;
				sprintf(m_FileName,"");

				float point[3] = {0.,0.,0.};

				//1 将第二块靶标上的点的坐标转换到第一块靶标上
				for (int j =/*0*/1;j< /*1*/CHESSNUM;j++)   //改动
				{
					cvMatMul(matH0_Invert_L,m_pcvH[j][0], m_pcvRV[j][0]);

					pD_RV = m_pcvRV[j][0]->data.db;
					step_RV = m_pcvRV[j][0]->step/sizeof(double);

					FILE * pf=NULL;
					sprintf(m_FileName,"%sH1-2[%d]",picsaveaddr_L,imgcnt+1);     //输出每幅图片中R = H1的逆*H2
					pf = fopen(m_FileName,"w");

					for (int n = 0;n < 3 ;n++)
					{

						fprintf(pf,"%.3f %.3f %.3f ", pD_RV[n*step_RV + 0],pD_RV[n*step_RV + 1],pD_RV[n*step_RV + 2]);
						fprintf(pf,"\r\n");

					}
					fclose(pf);
					pf=NULL;
					sprintf(m_FileName,"");

					pData = pcvObj_L->data.fl;
					step = pcvObj_L->step/sizeof(float);

					pData1 = _ObjPoints_L->data.fl;
					step11 = _ObjPoints_L->step/sizeof(float); 

					for (int n = 0;n < m_iPointCounts ;n++)
					{
						point[0] = pD_RV[0]*pData[n*step+0] + pD_RV[1]*pData[n*step+1] + pD_RV[2];
						point[1] = pD_RV[step_RV + 0]*pData[n*step+0] + pD_RV[step_RV + 1]*pData[n*step+1] + pD_RV[step_RV + 2];
						point[2] = pD_RV[2*step_RV + 0]*pData[n*step+0] + pD_RV[2*step_RV + 1]*pData[n*step+1] + pD_RV[2*step_RV + 2];

						cont = imgcnt*CHESSNUM*m_iPointCounts*step11 + n*step11 + j*m_iPointCounts*step11 ;
						pData1[cont + 0] = point[0]/point[2];   //待改变
						pData1[cont + 1] = point[1]/point[2];
						pData1[cont + 2] = 0 ;

					}

					pData = NULL;
					pData1 = NULL;
					pD_RV = NULL;

				}

				//............2 以坐标点128个，标定摄像机参数................
				//2.1 将靶标2的标记点的转换后的物理坐标，加到统一的坐标系（三维坐标）中（放在最前面）
				pData = pcvObj_L->data.fl;
				step = pcvObj_L->step/sizeof(float);

				pData1 = _ObjPoints_L->data.fl;
				step11 = _ObjPoints_L->step/sizeof(float);

				for (int n=0;n<m_iPointCounts;n++)
				{
					cont = imgcnt*CHESSNUM*m_iPointCounts*step11 + n*step11 /*+ m_iPointCounts*step11*/; //改动   //放在每张图片物理坐标的第一个  
					pData1[cont + 0] = pData[n*step+0];
					pData1[cont + 1] = pData[n*step+1];
					pData1[cont + 2] = 0;
				}
				pData = NULL;

				sprintf(m_FileName,"%sthree_dimen_obj02_new[%d].txt",picsaveaddr_L,imgcnt+1);   //将每幅图片的新的物理坐标输出
				pf =fopen(m_FileName,"w");

				for (int n = 0;n < m_iPointCounts*CHESSNUM ;n++)
				{
					cont = imgcnt*CHESSNUM*m_iPointCounts*step11 ;
					fprintf(pf,"%.3f %.3f %.3f    %d ", pData1[cont + n*step11 + 0],pData1[cont + n*step11 + 1],pData1[cont + n*step11 + 2],n);
					fprintf(pf,"\r\n");

				}
				fclose(pf);
				pf=NULL;
				pData1 = NULL;
				sprintf(m_FileName,"");

			}
			else if (n == 1)
			{
				sprintf(FileNameTail,"_RB");

				for (int i=0;i<CHESSNUM;i++)
				{
					m_pcvH[i][1] = cvCreateMat(3,3,CV_64FC1);
					m_pcvRV[i][1] = cvCreateMat(3,3,CV_64FC1);
					pcvImg_R[i] = cvCreateMat(m_iPointCounts,CHESSNUM,CV_32FC1);
				}


				sprintf(m_FileName,"%sChessboardImg1%s (%d).bmp",picreadaddress_R,FileNameTail,imgcnt + 1);	

				//彩色图片
				IplImage* img1=cvLoadImage(m_FileName,-1);

				img_R= cvCreateImage(cvGetSize(img1),8,1);
				cvCvtColor(img1,img_R,CV_BGR2GRAY);

				cvReleaseImage(&img1);
				img1 = NULL;

				pTemp = cvCreateImage(cvGetSize(img_R),IPL_DEPTH_8U,1);
				uchar*ptr=(uchar*)img_R->imageData;


				int height=img_R->height;
				int width=img_R->width;
				int pixpbyte=img_R->nChannels;

				_ImgSize_R = cvSize(width,height);



				int iGrayVar = 0;
				if (ithreshold_R <= 50)
				{
					iGrayVar = ithreshold_R;
				}
				else if (ithreshold_R <= 100)
				{
					iGrayVar = ithreshold_R*2/3;
				}
				else  
				{
					iGrayVar = ithreshold_R/2;
				}

				int iMemThres = ithreshold_R-iGrayVar;
				if (iMemThres <= 0)
				{
					iMemThres = 0;
				}

				FindCircleCenter(img_R,iMemThres,iMemThres+iGrayVar,(iGrayVar/10),50,10000,0.05f/*凸度*/,0.6f/*惯性率*/,0.6f/*圆度*/,key_points_R,m_iPointCounts,CHESSNUM,1);
				cvCanny(img_R, pTemp, (iMemThres+iGrayVar/2),(iMemThres+iGrayVar/2)*2, 3 );

				sprintf(m_FileName,"");
				sprintf(m_FileName,"%scvCannyImg%d%s.bmp",picsaveaddr_R,imgcnt+1,FileNameTail);
				cvSaveImage(m_FileName,pTemp);

				CalCircleCentrePoint(pTemp,key_points_R,m_iPointCounts,CHESSNUM);
				cvReleaseImage(&pTemp);

				//将圆心画在圆上
				for (int i=0;i<key_points_R.size();i++)
				{
					int x=(int) key_points_R.at(i).pt.x;
					int y=(int) key_points_R.at(i).pt.y;
					*(ptr+y*width+x)=(uchar)255;

				}


				sprintf(m_FileName,"");
				sprintf(m_FileName,"%sDetectedImg%d%s.bmp",picsaveaddr_R,imgcnt+1,FileNameTail);
				cvSaveImage(m_FileName,img_R);

				cout<<"记录所有找到的标记点坐标，并存为文件"<<endl;
				//记录下第一块板上的所有坐标
				FILE * pf=NULL;
				sprintf(m_FileName,"");
				sprintf(m_FileName,"%sgridcorner[%d].txt",picsaveaddr_R,imgcnt+1);//其中0：左相机；1:右相机
				pf=fopen(m_FileName,"w");
				int countNumber=0;

				for (int i=0;i<key_points_R.size();i++)
				{

					ST_CIRCLE_POINT CirclePoints;
					CirclePoints.x=key_points_R.at(i).pt.x;
					CirclePoints.y=key_points_R.at(i).pt.y;
					CirclePoints.size=key_points_R.at(i).size;
					CirclePoints.distance=0;
					countNumber++;
					Circle_Ponits_R.push_back(CirclePoints);
					fprintf(pf,"%.3f %.3f %d",key_points_R.at(i).pt.x,key_points_R.at(i).pt.y,countNumber);
					fprintf(pf,"\r\n");

				}

				fclose(pf);

				pf=NULL;


				//获得标定板的标记点坐标
				//获得标定板的标记点坐标
				GetOrignPoints(Circle_Ponits_R,CHESSNUM,pstOrginPoints_R);

				//获得标定板的原点坐标
				GetUpLeftPoints(Circle_Ponits_R,pstOrginPoints_R,pstUpLeftPoints_R,CHESSNUM);

				SortPoints(Computer_points_R,Circle_Ponits_R,pstUpLeftPoints_R,CHESSNUM,m_iPointCounts);
				//将顺序画在图上
				for (int k=0;k<CHESSNUM;k++)
				{
					if ((int)Computer_points_R[k].size()!=m_iPointCounts)
					{
						CString str;
						str.Format(_T("第%d张标定图片不能找到角点"),imgcnt+1);
						MessageBox(str);
					}
					for (int j=0;j<(int)Computer_points_R[k].size()-1;j++)
					{

						float x= Computer_points_R[k].at(j).x;
						float y= Computer_points_R[k].at(j).y;

						float x1= Computer_points_R[k].at(j+1).x;
						float y1= Computer_points_R[k].at(j+1).y;

						CvPoint pt1;
						CvPoint pt2;
						pt1.x=x;
						pt1.y=y;
						pt2.x=x1;
						pt2.y=y1;

						cvLine(img_R,pt1,pt2,cvScalar(255,0,0),1,8,0);

					}
				}

				sprintf_s(m_FileName,"");
				sprintf_s(m_FileName,"%sDetectedSortedImg%d%s.bmp",picsaveaddr_R,imgcnt+1,FileNameTail);
				cvSaveImage(m_FileName,img_R);

				//将排列好的标记点赋值给矩阵，

				float *pfdata=NULL;
				float *pfdata1 = NULL;

				int istep=0;
				int istep1 =0;

				int cont =0;
				for (int chessnumber=0;chessnumber<CHESSNUM;chessnumber++)
				{   

					//保存排序好的标记点的像素坐标
					sprintf_s(m_FileName,"");
					sprintf_s(m_FileName,"%sImgPoints_[%d]_%d.txt",picsaveaddr_R,imgcnt+1,chessnumber);///////
					pf=fopen(m_FileName,"w");

					pfdata=_ImgPoints_R->data.fl;
					istep=_ImgPoints_R->step/sizeof(float);

					pfdata1 = pcvImg_R[chessnumber]->data.fl;       //用来求得每个图片上的小靶标的透视矩阵
					istep1 = pcvImg_R[chessnumber]->step/sizeof(float);

					for (int i=0;i< Computer_points_R[chessnumber].size();i++)
					{
						pfdata1[i*istep1] = Computer_points_R[chessnumber].at(i).x;
						pfdata1[i*istep1+1] = Computer_points_R[chessnumber].at(i).y;

						//给矩阵赋值
						cont = imgcnt*CHESSNUM*m_iPointCounts*istep + chessnumber*m_iPointCounts*istep + i*istep;
						pfdata[cont] = Computer_points_R[chessnumber].at(i).x;
						pfdata[cont + 1] = Computer_points_R[chessnumber].at(i).y;

						fprintf(pf,"%.3f %.3f %d ",pfdata[cont],pfdata[cont + 1],i);
						fprintf(pf,"\r\n");

					}
					pfdata=NULL;
					pfdata1 = NULL;

					fclose(pf);
					pf=NULL;
					sprintf(m_FileName,"");

				}

				piPointNum_R[imgcnt]=m_iPointCounts*CHESSNUM;
				cvSetData(iPointCont1_R,&m_iPointCounts,sizeof(int));

				Circle_Ponits_R.clear();
				for (int chessnumber=0;chessnumber<CHESSNUM;chessnumber++)
				{
					Computer_points_R[chessnumber].clear();

				}

				//0对每幅图片的多块标定板分别获得各自的透视矩阵H[i]
				for (int i=0;i<CHESSNUM;i++)
				{

					cvCalibrateCamera2_H(pcvObj_R,pcvImg_R[i],iPointCont1_R,_ImgSize_R,m_pcvH[i][1],i,0);

				}

				//1 求得各个矩阵到H1的转换矩阵
				CvMat* matH0_Invert_R = cvCreateMat(3, 3, CV_64FC1);
				cvZero(matH0_Invert_R);
				cvInvert(m_pcvH[0][1], matH0_Invert_R,CV_LU);//改动



				float * pData = NULL;    //接收m_pcvObjectPoints[j]
				int step = 0;
				float * pData1 = NULL;   //m_pcvThree_dimeObj
				int step11 = 0;
				double * pD_RV = NULL;   //接收m_pcvRV[j]
				int step_RV = 0;

				double* pD_R = NULL;   //临时输出H0的逆
				int step_R = 0;

				pD_R= matH0_Invert_R->data.db;
				step_R = matH0_Invert_R->step/sizeof(double);


				sprintf(m_FileName,"%sH0_INV[%d]",picsaveaddr_R,imgcnt+1);   //输出每幅图片中靶标1的H1的逆
				pf =fopen(m_FileName,"w");

				for (int j = 0;j<3;j++)
				{
					fprintf(pf,"%.3f %.3f %.3f " , pD_R[j*step_R+0],pD_R[j*step_R+1],pD_R[j*step_R+2]);
					fprintf(pf,"\r\n");
				} 
				fclose(pf);
				pf=NULL;
				pD_R =NULL;
				sprintf(m_FileName,"");

				float point[3] = {0.,0.,0.};

				//1 将第二块靶标上的点的坐标转换到第一块靶标上
				for (int j =/*0*/1;j< /*1*/CHESSNUM;j++)   //改动
				{
					cvMatMul(matH0_Invert_R,m_pcvH[j][1], m_pcvRV[j][1]);

					pD_RV = m_pcvRV[j][1]->data.db;
					step_RV = m_pcvRV[j][1]->step/sizeof(double);

					FILE * pf=NULL;
					sprintf(m_FileName,"%sH1-2[%d]",picsaveaddr_R,imgcnt+1);     //输出每幅图片中R = H1的逆*H2
					pf = fopen(m_FileName,"w");

					for (int n = 0;n < 3 ;n++)
					{

						fprintf(pf,"%.3f %.3f %.3f ", pD_RV[n*step_RV + 0],pD_RV[n*step_RV + 1],pD_RV[n*step_RV + 2]);
						fprintf(pf,"\r\n");

					}
					fclose(pf);
					pf=NULL;
					sprintf(m_FileName,"");

					pData = pcvObj_R->data.fl;
					step = pcvObj_R->step/sizeof(float);

					pData1 = _ObjPoints_R->data.fl;
					step11 = _ObjPoints_R->step/sizeof(float); 

					for (int n = 0;n < m_iPointCounts ;n++)
					{
						point[0] = pD_RV[0]*pData[n*step+0] + pD_RV[1]*pData[n*step+1] + pD_RV[2];
						point[1] = pD_RV[step_RV + 0]*pData[n*step+0] + pD_RV[step_RV + 1]*pData[n*step+1] + pD_RV[step_RV + 2];
						point[2] = pD_RV[2*step_RV + 0]*pData[n*step+0] + pD_RV[2*step_RV + 1]*pData[n*step+1] + pD_RV[2*step_RV + 2];

						cont = imgcnt*CHESSNUM*m_iPointCounts*step11 + n*step11 + j*m_iPointCounts*step11 ;
						pData1[cont + 0] = point[0]/point[2];   //待改变
						pData1[cont + 1] = point[1]/point[2];
						pData1[cont + 2] = 0 ;

					}

					pData = NULL;
					pData1 = NULL;
					pD_RV = NULL;

				}

				//............2 以坐标点128个，标定摄像机参数................
				//2.1 将靶标2的标记点的转换后的物理坐标，加到统一的坐标系（三维坐标）中（放在最前面）
				pData = pcvObj_R->data.fl;
				step = pcvObj_R->step/sizeof(float);

				pData1 = _ObjPoints_R->data.fl;
				step11 = _ObjPoints_R->step/sizeof(float);

				for (int n=0;n<m_iPointCounts;n++)
				{
					cont = imgcnt*CHESSNUM*m_iPointCounts*step11 + n*step11 /*+ m_iPointCounts*step11*/; //改动   //放在每张图片物理坐标的第一个  
					pData1[cont + 0] = pData[n*step+0];
					pData1[cont + 1] = pData[n*step+1];
					pData1[cont + 2] = 0;
				}
				pData = NULL;

				sprintf(m_FileName,"%sthree_dimen_obj02_new[%d].txt",picsaveaddr_R,imgcnt+1);   //将每幅图片的新的物理坐标输出
				pf =fopen(m_FileName,"w");

				for (int n = 0;n < m_iPointCounts*CHESSNUM ;n++)
				{
					cont = imgcnt*CHESSNUM*m_iPointCounts*step11 ;
					fprintf(pf,"%.3f %.3f %.3f    %d ", pData1[cont + n*step11 + 0],pData1[cont + n*step11 + 1],pData1[cont + n*step11 + 2],n);
					fprintf(pf,"\r\n");

				}
				fclose(pf);
				pf=NULL;
				pData1 = NULL;
				sprintf(m_FileName,"");


			}

		}
	}

	//释放资源

	destroyWindow("SimpleBlobDetector_L");

	cvSetData(_PointsCnt_L,piPointNum_L,sizeof(int));

	//释放空间

	cvReleaseMat(&pcvObj_L);
	cvReleaseMat(&iPointCont1_L);
	obj_L = NULL;
	for (int i = 0;i<CHESSNUM;i++)
	{
		cvReleaseMat(&pcvImg_L[i]);
		cvReleaseMat(&m_pcvH[i][0]);
		cvReleaseMat(&m_pcvRV[i][0]);
	}

	//将结果传给输出变量
	imgsize_L = _ImgSize_L;
	cvConvert(_ImgPoints_L,imgpoints_L);
	cvConvert(_ObjPoints_L,objpoints_L);
	cvConvert(_PointsCnt_L,everypicpointcnt_L);

	cvReleaseMat(&_ImgPoints_L);
	cvReleaseMat(&_ObjPoints_L);
	cvReleaseMat(&_PointsCnt_L);


	//释放资源

	destroyWindow("SimpleBlobDetector_R");

	cvSetData(_PointsCnt_R,piPointNum_R,sizeof(int));

	//释放空间

	cvReleaseMat(&pcvObj_R);
	cvReleaseMat(&iPointCont1_R);
	obj_R = NULL;
	for (int i = 0;i<CHESSNUM;i++)
	{
		cvReleaseMat(&pcvImg_R[i]);
	}

	for (int i=0;i<CHESSNUM;i++)
	{
		cvReleaseMat(&m_pcvH[i][1]);
		cvReleaseMat(&m_pcvRV[i][1]);
	}
	//将结果传给输出变量
	imgsize_R = _ImgSize_R;
	cvConvert(_ImgPoints_R,imgpoints_R);
	cvConvert(_ObjPoints_R,objpoints_R);
	cvConvert(_PointsCnt_R,everypicpointcnt_R);

	cvReleaseMat(&_ImgPoints_R);
	cvReleaseMat(&_ObjPoints_R);
	cvReleaseMat(&_PointsCnt_R);


	return TRUE;
}


/***************************************************************************
* 函数名称：   CCalibrationView::FindCircleCenter
* 摘    要：   需找标记点点圆心
* 全局影响：   public 
* 参    数：   IplImage * img
* 参    数：   int minThreshold
* 参    数：   int MaxThreshold
* 参    数：   int shresholdStep
* 参    数：   int minArea
* 参    数：   int MaxArea
* 参    数：   float minConvexity
* 参    数：   float minInertiaRatio
* 参    数：   float minCircularity
* 参    数：   vector<KeyPoint> & Vct_ponits
* 参    数：   int ipointscounts
* 参    数：   int iChesscount
* 返回值：     BOOL
* 
* 修改记录：
*  [日期]     [作者/修改者]  [修改原因]
*2015/07/25       罗佳丽          添加
***************************************************************************/
BOOL CCalibrateView::FindCircleCenter(IplImage*img, int minThreshold,int MaxThreshold,int shresholdStep, int minArea,int MaxArea,float minConvexity,float minInertiaRatio,float minCircularity,vector<KeyPoint>&Vct_ponits,int ipointscounts,int iChesscount,int flag)
{
	CvSize m_imageSize;
	m_imageSize=cvSize(img->width,img->height); 
	Mat mtx(img);

	Mat binaryImage(img);
	SimpleBlobDetector::Params params; 
	params.filterByCircularity=true;
	params.minThreshold = minThreshold;     //最小的阈值
	params.maxThreshold = MaxThreshold;    //最大的阈值
	params.thresholdStep = shresholdStep;    //多重阈值分割的阈值步长
	params.minArea =minArea;              //检测圆的面积最小值
	params.maxArea = MaxArea;              //检测圆的面积最大值
	params.minConvexity = minConvexity;   //斑点的凸度
	params.minInertiaRatio = minInertiaRatio;//斑点的最小惯性率
	params.blobColor=0;                    //只提取黑色斑点
	params.minCircularity=minCircularity;  //斑点的最小圆度
	SimpleBlobDetector detector(params);  
	detector.detect(mtx,Vct_ponits);  
	if (Vct_ponits.size()<ipointscounts*iChesscount+2)
	{
		cout<<"图片中未找到所有标记点"<<endl;
		return FALSE;
	}
	else
	{

		cout<<"成功找到所有标记点"<<endl;
		drawKeypoints( mtx, Vct_ponits, mtx, Scalar(0,0,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS ); 
		IplImage imgTemp =IplImage(mtx);

		//显示图片修改

		if (flag == 0)
		{
			namedWindow("SimpleBlobDetector_L",CV_WINDOW_NORMAL); 
			imshow("SimpleBlobDetector_L",mtx);
			cvWaitKey(0);
			destroyWindow("SimpleBlobDetector_L");
		}
		else if(flag == 1)
		{
			namedWindow("SimpleBlobDetector_R",CV_WINDOW_NORMAL); 
			imshow("SimpleBlobDetector_R",mtx);
			cvWaitKey(0);
			destroyWindow("SimpleBlobDetector_R");
		}

		cout<<"标记圆精确定位"<<endl;

		return  TRUE;
	}
}


/***************************************************************************
* 函数名称：   CCalibrateView::CalCircleCentrePoint
* 摘    要：   寻找 圆心
* 全局影响：   public 
* 参    数：   IplImage *  & pImage
* 参    数：   vector<KeyPoint> & Vct_ponits
* 参    数：   int ipointscounts
* 参    数：   int iChesscount
* 返回值：     void
* 
* 修改记录：
*  [日期]     [作者/修改者]  [修改原因]
*2015/07/25       罗佳丽          添加
***************************************************************************/
void CCalibrateView::CalCircleCentrePoint(IplImage* &pImage,vector<KeyPoint>&Vct_ponits,int ipointscounts,int iChesscount)
{
	for (int iCirclenCounts = 0; iCirclenCounts < Vct_ponits.size(); iCirclenCounts++)
	{
		/// 粗略的圆心坐标
		int iCentreX =(int)Vct_ponits.at(iCirclenCounts).pt.x;
		int iCentreY =(int)Vct_ponits.at(iCirclenCounts).pt.y;

		///扫描到的边缘点集
		int m_iPos[72][2];

		memset(m_iPos, 0, sizeof(int) * 72 * 2);

		/// 以0°（水平向右）为起点，沿逆时针方向每隔5°在圆周上取一点，并向圆心方向搜索，若某一点的像素值为255时，则为边缘点。
		for (int iAngle = 0; iAngle < 360; iAngle += 5)
		{
			for (int iR = 1; iR < 25; iR++)
			{
				int iX = (int)(iCentreX + (iR * cos((double)iAngle * PI / 180) + (((iAngle >= 0 && iAngle <= 90) || (iAngle >= 270 && iAngle < 360)) ? 0.5 : -0.5)));		// 待判断点的X坐标，式中的0.5为修正值
				int iY = (int)(iCentreY - (iR * sin((double)iAngle * PI / 180) + ((iAngle >= 0 && iAngle < 180) ? 0.5 : -0.5)));		// 待判断点的Y坐标，式中的0.5为修正值

				/// 判断是否超出图像的范围
				if (iX < 0 || iX >= pImage->width || iY < 0 || iY >= pImage->height)
				{
					m_iPos[iAngle / 5][0] = iX - iCentreX;					// 扫描到的边缘点集（iPos）是以（iCentreX，iCentreY）为坐标中心的坐标
					m_iPos[iAngle / 5][1] = iY - iCentreY;
					break;
				}

				uchar ucVal = (uchar)*(pImage->imageData + iY * pImage->widthStep + iX * pImage->nChannels + 0);
				if (255 == ucVal)
				{
					m_iPos[iAngle / 5][0] = iX - iCentreX;					// 扫描到的边缘点集（iPos）是以（iCentreX，iCentreY）为坐标中心的坐标
					m_iPos[iAngle / 5][1] = iY - iCentreY;
					break;
				}
			}
		}

		/// 采用交换排序法对扫描到的边缘点集按其到圆心的距离由近到远进行排序
		for (int i = 0; i < 72; i++)
		{
			double iDistance = sqrt((double)(m_iPos[i][0]* m_iPos[i][0] + m_iPos[i][1]/* -iCentreY*/ * m_iPos[i][1]));
			for (int j = i + 1; j < 72; j++)
			{
				double iTemp = sqrt((double)(m_iPos[j][0]* m_iPos[j][0] + m_iPos[j][1]/* -iCentreY*/ * m_iPos[j][1]));
				if (iTemp - iDistance < 0)
				{
					int x = m_iPos[i][0];
					int y = m_iPos[i][1];
					m_iPos[i][0] = m_iPos[j][0];
					m_iPos[i][1] = m_iPos[j][1];
					m_iPos[j][0] = x;
					m_iPos[j][1] = y;

					iDistance = sqrt((double)(m_iPos[i][0] * m_iPos[i][0] + m_iPos[i][1]/* -iCentreY*/ * m_iPos[i][1]));
				}
			}
		}

		/// 取位于iPos中间的36个点，用于拟合边缘圆。
		int iFitPos[36][2];				// 拟合点数组
		memset(iFitPos, 0, sizeof(int) * 36 * 2);
		for (int i = 0; i < 36; i++)
		{
			iFitPos[i][0] = m_iPos[24 + i][0];
			iFitPos[i][1] = m_iPos[24 + i][1];
		}

		/// 以下是最小二乘拟合求解圆心坐标及半径

		/// 矩阵A
		CvMat* matA = cvCreateMat(3, 3, CV_64FC1);
		cvZero(matA);
		for (int i = 0; i < 36; i++)
		{
			*(matA->data.db + 0) += iFitPos[i][0] * iFitPos[i][0];
			*(matA->data.db + 1) += iFitPos[i][0] * iFitPos[i][1];
			*(matA->data.db + 2) += iFitPos[i][0];
			*(matA->data.db + 3) += iFitPos[i][0] * iFitPos[i][1];
			*(matA->data.db + 4) += iFitPos[i][1] * iFitPos[i][1];
			*(matA->data.db + 5) += iFitPos[i][1];
			*(matA->data.db + 6) += iFitPos[i][0];
			*(matA->data.db + 7) += iFitPos[i][1];
		}
		*(matA->data.db + 8) = 36;

		/// 矩阵B
		CvMat* matB = cvCreateMat(3, 1, CV_64FC1);
		cvZero(matB);
		for (int i = 0; i < 36; i++)
		{
			*(matB->data.db + 0) -= iFitPos[i][0] * iFitPos[i][0] * iFitPos[i][0] + iFitPos[i][1] * iFitPos[i][1] * iFitPos[i][0];
			*(matB->data.db + 1) -= iFitPos[i][0] * iFitPos[i][0] * iFitPos[i][1] + iFitPos[i][1] * iFitPos[i][1] * iFitPos[i][1];
			*(matB->data.db + 2) -= iFitPos[i][0] * iFitPos[i][0] + iFitPos[i][1] * iFitPos[i][1];
		}

		/// 求A的逆矩阵
		CvMat* matA_Invert = cvCreateMat(3, 3, CV_64FC1);
		cvZero(matA_Invert);
		cvInvert(matA, matA_Invert);

		/// 求解圆心坐标及半径
		CvMat* matResult = cvCreateMat(3, 1, CV_64FC1);
		cvMatMul(matA_Invert, matB, matResult);

		double dbResultX = -matResult->data.db[0] / 2;
		double dbResultY = -matResult->data.db[1] / 2;
		double dbResultR = sqrt(dbResultX * dbResultX + dbResultY * dbResultY - matResult->data.db[2]);

		Vct_ponits.at(iCirclenCounts).pt.x = (float)(dbResultX+iCentreX);
		Vct_ponits.at(iCirclenCounts).pt.y=dbResultY+iCentreY;

		if (dbResultR!=0)
		{
			Vct_ponits.at(iCirclenCounts).size=dbResultR*dbResultR;
		}

		else
		{
			Vct_ponits.at(iCirclenCounts).size=Vct_ponits.at(iCirclenCounts).size*Vct_ponits.at(iCirclenCounts).size;
		}
	}
}


/***************************************************************************
* 函数名称：   CCalibrationView::GetChessborardMark
* 摘    要：   获得棋盘标记点坐标
* 全局影响：   public 
* 参    数：   vector<ST_CIRCLE_POINT> & Vct_points
* 参    数：   int m_iChessboardCount
* 参    数：   PST_ORGIN_POINTS & pst_orginPoints
* 返回值：     void
* 
* 修改记录：
*  [日期]     [作者/修改者]  [修改原因]
*2015/07/25       罗佳丽          添加
***************************************************************************/
void CCalibrateView::GetOrignPoints(vector<ST_CIRCLE_POINT>&Vct_points,int chessboardcount,PST_ORGIN_POINTS &pst_orginPoints)
{
	//分配空间 
	pst_orginPoints = new ST_ORGIN_POINTS[chessboardcount];
	int temp;
	for (int i=0;i<chessboardcount;i++)
	{

		float sizeMax=0.0; 
		temp=0;
		int iCount = (int)Vct_points.size();
		for(int j = 0; j < iCount; j++ )
		{
			//按照圆面积的大小 棋盘的标记点      
			if (Vct_points.at(j).size>sizeMax)
			{
				temp=j;
				sizeMax=Vct_points.at(j).size;

			}


		}

		//===============================================================添加保护
		if ( iCount != 0 )
		{
			pst_orginPoints[i].x = Vct_points.at(temp).x;
			pst_orginPoints[i].y = Vct_points.at(temp).y;
			pst_orginPoints[i].size = Vct_points.at(temp).size;
			Vct_points.erase( Vct_points.begin()+temp );
		}

	}
}


/***************************************************************************
* 函数名称：   CCalibrationView::GetUpLeftPoints
* 摘    要：   获得靶标左上原点坐标
* 全局影响：   public 
* 参    数：   vector<ST_CIRCLE_POINT> & Vct_points
* 参    数：   PST_ORGIN_POINTS & pst_OringPonit
* 参    数：   int m_iChessboardCount
* 返回值：     bool
* 
* 修改记录：
*  [日期]     [作者/修改者]  [修改原因]
*2015/07/25       罗佳丽          添加
***************************************************************************/
BOOL CCalibrateView::GetUpLeftPoints(vector<ST_CIRCLE_POINT>&Vct_points,PST_ORGIN_POINTS& pst_OringPonit,PST_UpLeft_POINTS& pst_UpLeftPoints,int chessboardcount)
{
	if (pst_OringPonit==NULL || Vct_points.size()==0)
	{
		return false;
	}

	pst_UpLeftPoints = new ST_UpLeft_POINTS[chessboardcount];

	//暂存变量，用来保存位置
	int  temp;
	for (int i=0;i<chessboardcount;i++)
	{
		temp=0;
		//
		double xcoordinate=pst_OringPonit[i].x;
		double ycoordinate=pst_OringPonit[i].y;


		//计算距离
		for (int j=0;j< Vct_points.size();j++)
		{

			double dx=fabs(Vct_points.at(j).x-xcoordinate);
			double dy=fabs(Vct_points.at(j).y-ycoordinate);
			Vct_points.at(j).distance= sqrt(dx*dx+dy*dy);
		}

		//计算最小距离
		double distance;
		distance=Vct_points.at(0).distance;

		for (int k=0;k< Vct_points.size();k++)
		{
			if (Vct_points.at(k).distance<distance)
			{
				distance=Vct_points.at(k).distance;
				//最小值在向量中的位置
				temp=k;
			}
		}

		pst_UpLeftPoints[i].x=Vct_points.at(temp).x;
		pst_UpLeftPoints[i].y=Vct_points.at(temp).y;
		pst_UpLeftPoints[i].size=Vct_points.at(temp).size;
	}
	return true;
}


/***************************************************************************
* 函数名称：   CCalibrationView::SortPoints
* 摘    要：   对一副图像的所有特征点进行排序
* 全局影响：   public 
* 参    数：   VEC_POINTS & Vec_points
* 参    数：   vector<ST_CIRCLE_POINT> & Circle_Ponits
* 参    数：   PST_ORGIN_POINTS Pst_orgin_points
* 参    数：   int Chessnum
* 参    数：   int Pointsnum
* 返回值：     void
* 
* 修改记录：
*  [日期]     [作者/修改者]  [修改原因]
*2015/07/25       罗佳丽          添加
***************************************************************************/
void CCalibrateView::SortPoints(VEC_POINTS&Vec_points,vector<ST_CIRCLE_POINT>&Circle_Ponits ,PST_UpLeft_POINTS Pst_UpLeft_points,int Chessnum,int Pointsnum)
{
	//图像原点
	ST_UpLeft_POINTS stUpperLeft;
	//距离原点最近的两点
	ST_CIRCLE_POINT NearstPoints[2];
	//两快标定板特征点排序
	for (int ChessNum = 0; ChessNum < Chessnum; ChessNum++)
	{


		stUpperLeft.x=Pst_UpLeft_points[ChessNum].x;
		stUpperLeft.y=Pst_UpLeft_points[ChessNum].y;
		stUpperLeft.size=Pst_UpLeft_points[ChessNum].size;


		//计算所有点距离原点的位置，距1号的原点
		CalcPonitDistance(stUpperLeft,Circle_Ponits);

		//计算最小距离
		double MinValue = GetMinDistance(Circle_Ponits);

		int Count=0;
		for (int i=0;i< Circle_Ponits.size();i++)
		{

			if (Circle_Ponits.at(i).distance<=MinValue+4 && Circle_Ponits.at(i).distance>0.0)
			{

				NearstPoints[Count].x=Circle_Ponits.at(i).x;
				NearstPoints[Count].y=Circle_Ponits.at(i).y;
				NearstPoints[Count].size=Circle_Ponits.at(i).size;

				NearstPoints[Count].distance=Circle_Ponits.at(i).distance;
				Count++;
			}
			if(Count==2)
				break;
		}

		//判断哪个是大点
		float dx1,dy1,dx,dy;

		if (NearstPoints[0].size > NearstPoints[1].size)
		{
			//竖直方向X轴变化率
			dx1=NearstPoints[0].x-stUpperLeft.x;


			//竖直方向Y轴变化率
			dy1=NearstPoints[0].y-stUpperLeft.y;
			//水平方向X轴变换率
			dx=NearstPoints[1].x-stUpperLeft.x;


			//水平方向Y轴变化率
			dy=NearstPoints[1].y-stUpperLeft.y;
		}
		else
		{
			//竖直方向X轴变化率
			dx=NearstPoints[0].x-stUpperLeft.x;


			//竖直方向Y轴变化率
			dy=NearstPoints[0].y-stUpperLeft.y;
			//水平方向X轴变换率
			dx1=NearstPoints[1].x-stUpperLeft.x;


			//水平方向Y轴变化率
			dy1=NearstPoints[1].y-stUpperLeft.y;
		}


		//原点坐标

		float x = Pst_UpLeft_points[ChessNum].x;
		float y = Pst_UpLeft_points[ChessNum].y;

		//点都存在 vector<ST_CIRCLE_POINT>Circle_Ponits里面 

		ST_CIRCLE_POINT ComputerPonits;
		ComputerPonits.x=x;
		ComputerPonits.y=y;


		//暂存
		float tempx=0;
		float tempy=0;

		tempx=x;
		tempy=y;
		//计算值与理论值匹配的点数目
		int findnumber=0;
		for (int i=0;i<m_iRowCount;i++)
		{
			for (int j=0;j<m_iRowCount;j++)
			{ 

				//寻找测量点中与理论点中最接近的点，并赋值
				int iChessSize = (int)Circle_Ponits.size();
				for (int k = 0;k < iChessSize; k++ )
				{
					if (fabs((Circle_Ponits.at(k).x-tempx))<=5 && fabs((Circle_Ponits.at(k).y-tempy))<=5)
					{
						ComputerPonits.dex = Circle_Ponits.at(k).x-tempx;
						ComputerPonits.dey = Circle_Ponits.at(k).y-tempy;
						tempx=Circle_Ponits.at(k).x;
						tempy=Circle_Ponits.at(k).y;


						ComputerPonits.x = tempx;
						ComputerPonits.y = tempy;

						Vec_points[ChessNum].push_back(ComputerPonits);

						tempx=tempx+dx;
						tempy=tempy+dy;

						findnumber++;

						//到下一行则变化量 

 						if (findnumber>0 && findnumber%m_iRowCount==0 )
						{
							tempx=Vec_points[ChessNum].at(findnumber-m_iRowCount).x+dx1;//竖直方向X
							tempy=Vec_points[ChessNum].at(findnumber-m_iRowCount).y+dy1;//竖直方向Y
						}

						Circle_Ponits.erase(Circle_Ponits.begin()+k);

						break;
					}
				}  
			}
		}
	}
}



/***************************************************************************
* 函数名称：   CCalibrationView::CalcPonitDistance
* 摘    要：   计算所有圆点阵列所有的点
* 全局影响：   public 
* 参    数：   ST_CIRCLE_POINT UpperLeftPonits
* 参    数：   vector<ST_CIRCLE_POINT> & Vct_Points
* 返回值：     BOOL
* 
* 修改记录：
*  [日期]     [作者/修改者]  [修改原因]
*2015/07/25       罗佳丽          添加
***************************************************************************/
BOOL CCalibrateView::CalcPonitDistance(ST_UpLeft_POINTS UpLeftPonits,vector<ST_CIRCLE_POINT>&Vct_Points )
{
	//double distance;
	double xcoordinate=UpLeftPonits.x;
	double ycoordinate=UpLeftPonits.y;
	if (Vct_Points.size()==0)
	{
		return false;
	}
	else
	{
		for (int i=0;i<Vct_Points.size();i++)
		{

			if (Vct_Points.at(i).x!=xcoordinate && Vct_Points.at(i).y!= ycoordinate)
			{

				double dx=fabs(Vct_Points.at(i).x-xcoordinate);
				double dy=fabs(Vct_Points.at(i).y-ycoordinate);

				Vct_Points.at(i).distance= sqrt(dx*dx+dy*dy);
			}

			else 
			{
				Vct_Points.at(i).distance=0.0;
			}

		}
	}

	return true;
}


/***************************************************************************
	* 函数名称：   CCalibrationView::GetMinDistance
	* 摘    要：   得到距原点最近的两个点的距离
	* 全局影响：   public 
	* 参    数：   vector<ST_CIRCLE_POINT> & Vct_ponits
	* 返回值：     double
	* 
	* 修改记录：
	*  [日期]     [作者/修改者]  [修改原因]
*2015/07/25       罗佳丽          添加
	***************************************************************************/
double CCalibrateView::GetMinDistance(vector<ST_CIRCLE_POINT>&Vct_ponits)
{
	double MinDistValue;
	if (Vct_ponits.size()==0)
	{
		return 0;
	}     
	else
	{
		MinDistValue=Vct_ponits.at(0).distance;
		for (int i=0;i < Vct_ponits.size();i++)
		{
			if (Vct_ponits.at(i).distance<MinDistValue && Vct_ponits.at(i).distance>0)
			{
				MinDistValue=Vct_ponits.at(i).distance;
			}

		}
	}
	return MinDistValue;
}

void CCalibrateView::cvCalibrateCamera2_H( const CvMat* objectPoints,
								  const CvMat* imagePoints, const CvMat* npoints,
								  CvSize imageSize, CvMat* H,int cnt,
								  int flags )
{
	const int NINTRINSIC = 12;
	Ptr<CvMat> matM, _m, _Ji, _Je, _err,h;
	CvLevMarq solver;
	double reprojErr = 0;

	double A[9], k[8] = {0,0,0,0,0,0,0,0};
	CvMat matA = cvMat(3, 3, CV_64F, A);
	int i , nimages, maxPoints = 0, ni = 0, total = 0, nparams, npstep;
	double aspectRatio = 0.;

	// 0. check the parameters & allocate buffers
	if( !CV_IS_MAT(objectPoints) || !CV_IS_MAT(imagePoints) ||
		!CV_IS_MAT(npoints))
		CV_Error( CV_StsBadArg, "One of required vector arguments is not a valid matrix" );

	if( imageSize.width <= 0 || imageSize.height <= 0 )
		CV_Error( CV_StsOutOfRange, "image width and height must be positive" );

	if( CV_MAT_TYPE(npoints->type) != CV_32SC1 ||
		(npoints->rows != 1 && npoints->cols != 1) )
		CV_Error( CV_StsUnsupportedFormat,
		"the array of point counters must be 1-dimensional integer vector" );

	nimages = npoints->rows*npoints->cols;
	npstep = npoints->rows == 1 ? 1 : npoints->step/CV_ELEM_SIZE(npoints->type);

	if( (CV_MAT_TYPE(H->type) != CV_32FC1 &&
		CV_MAT_TYPE(H->type) != CV_64FC1) ||
		H->rows != 3 || H->cols != 3 )
		CV_Error( CV_StsBadArg,
		"H parameters must be 3x3 floating-point matrix" );

	for( i =0; i < nimages; i++ )
	{
		ni = npoints->data.i[i*npstep];
		if( ni < 4 )
		{
			char buf[100];
			sprintf( buf, "The number of points in the view #%d is < 4", i );
			CV_Error( CV_StsOutOfRange, buf );
		}
		maxPoints = MAX( maxPoints, ni );
		total += ni;
	}

	matM = cvCreateMat( 1, total, CV_64FC3 );
	_m = cvCreateMat( 1, total, CV_64FC2 );
	h = cvCreateMat(3,3,CV_64FC1);

	cvConvertPointsHomogeneous( objectPoints, matM );
	cvConvertPointsHomogeneous( imagePoints, _m );

	nparams = NINTRINSIC + nimages*6;
	_Ji = cvCreateMat( maxPoints*2, NINTRINSIC, CV_64FC1 );
	_Je = cvCreateMat( maxPoints*2, 6, CV_64FC1 );
	_err = cvCreateMat( maxPoints*2, 1, CV_64FC1 );
	cvZero( _Ji );


	// 1. initialize intrinsic parameters & LM solver
	CvScalar mean, sdv;
	cvAvgSdv( matM, &mean, &sdv );
	if( fabs(mean.val[2]) > 1e-5 || fabs(sdv.val[2]) > 1e-5 )
		CV_Error( CV_StsBadArg,
		"For non-planar calibration rigs the initial intrinsic matrix must be specified" );
	/*for( i = 0; i < total; i++ )
		((CvPoint3D64f*)matM->data.db)[i].z = 0.;*/

	
	cvInitIntrinsicParams2D_H( matM, _m, npoints, imageSize, h, aspectRatio ,cnt);
	
	cvConvert(h,H);

}

void CCalibrateView::cvInitIntrinsicParams2D_H( const CvMat* objectPoints,
									 const CvMat* imagePoints, const CvMat* npoints,
									 CvSize imageSize, CvMat* cameraMatrix,
									 double aspectRatio,int cnt )
{
	Ptr<CvMat> matA, _b, _allH, _allK;

	int i, pos, nimages, total, ni = 0;
	double H[9];
	CvMat matH = cvMat( 3, 3, CV_64F, H );

	assert( CV_MAT_TYPE(npoints->type) == CV_32SC1 &&
		CV_IS_MAT_CONT(npoints->type) );
	nimages = npoints->rows + npoints->cols - 1;

	if( (CV_MAT_TYPE(objectPoints->type) != CV_32FC3 &&CV_MAT_TYPE(objectPoints->type) != CV_64FC3) ||
		(CV_MAT_TYPE(imagePoints->type) != CV_32FC2 &&CV_MAT_TYPE(imagePoints->type) != CV_64FC2) )
		CV_Error( CV_StsUnsupportedFormat, "Both object points and image points must be 2D" );

	if( objectPoints->rows != 1 || imagePoints->rows != 1 )
		CV_Error( CV_StsBadSize, "object points and image points must be a single-row matrices" );

	_allH = cvCreateMat( nimages, 9, CV_64F );

	total = cvRound(cvSum(npoints).val[0]);

	// extract vanishing points in order to obtain initial value for the focal length
	for( i = 0, pos = 0; i < nimages; i++, pos += ni )
	{
		ni = npoints->data.i[i];
		CvMat _m, matM;
		cvGetCols( objectPoints, &matM, pos, pos + ni );
		cvGetCols( imagePoints, &_m, pos, pos + ni );

		cvFindHomography( &matM, &_m, &matH );
		memcpy( _allH->data.db + i*9, H, sizeof(H) );

	}

	cvConvert( &matH, cameraMatrix );
}


void CCalibrateView::OnBnClickedBtnGetp()
{
	// TODO: 在此添加控件通知处理程序代码
	MessageBox(_T("图片获取成功!"),NULL,MB_OK);
}
