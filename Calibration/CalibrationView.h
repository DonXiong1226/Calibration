
// CalibrateView.h : CCalibrateView 类的接口
//

#pragma once

#include <iostream>
#include "define.h"
#include "resource.h"
#include "cv.h"
#include "highgui.h"
#include<vector>


using namespace cv;
#define CHESSNUM 2    //标定板块数
#define CAMERANUM 2   //双目标定时摄像机数

//控制台  头文件
#include<conio.h>
#include "afxwin.h"
#include "afxcmn.h"

class CCalibrateView : public CFormView
{
protected: // 仅从序列化创建
	CCalibrateView();
	DECLARE_DYNCREATE(CCalibrateView)

public:
	enum{ IDD = IDD_CALIBRATE_FORM };

// 特性
public:
	CCalibrateDoc* GetDocument() const;

// 操作
public:

// 重写
public:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持
	virtual void OnInitialUpdate(); // 构造后第一次调用

// 实现
public:
	virtual ~CCalibrateView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// 生成的消息映射函数
protected:
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()

public:

CvMat* m_pcvBig_dimeObj_Left;		// 左相机转化为大靶标的物理坐标矩阵
CvMat* m_pcvBig_dimeObj_Right;		// 右相机转化为大靶标的物理坐标矩阵
CvMat* m_pcvImgPoint_Left;			// 左目相机M幅图像中所有角点像素坐标矩阵
CvMat* m_pcvImgPoint_Right;			// 右目相机M幅图像中所有角点像素坐标矩阵

CvMat* m_pcvTestImg_Left;			//保存测试图像中标记点的像素坐标
CvMat* m_pcvTestImg_Right;			//保存测试图像中标记点的像素坐标

CvMat* m_pcvPointCounts;			// 所有用于标定相机的图像的角点个数组成的矩阵(左右相机一样)
CvSize m_imageSize_Left;			//左目相机M幅图像的大小
CvSize m_imageSize_Right;			//右目相机M幅图像的大小

// 摄像机内参数
CvMat* m_pcvIntrinsicMat_Left;		// 左摄像机参数矩阵
CvMat* m_pcvIntrinsicMat_Right;		// 右摄像机参数矩阵
CvMat* m_pcvDistortionCoeffs_Left;  // 左摄像机畸变系数
CvMat* m_pcvDistortionCoeffs_Right; //右摄像机畸变系数

// 摄像机外参数
CvMat* m_pcvRotationVct_Left;		// 棋盘绕左目摄像机坐标系坐标轴的旋转矩阵
CvMat* m_pcvRotationVct_Right;		// 棋盘绕右目摄像机坐标系坐标轴的旋转矩阵
CvMat* m_pcvTranslationVct_Left;	// 棋盘绕左目摄像机坐标系坐标轴的平移矩阵
CvMat* m_pcvTranslationVct_Right;	// 棋盘绕右目摄像机坐标系坐标轴的平移矩阵

// 校正矩阵
CvMat* m_pcvMapX_Left;				//左目相机x方向校正映射表，即mapx1
CvMat* m_pcvMapX_Right;				//右目相机x方向校正映射表，即mapy1
CvMat* m_pcvMapY_Left;				//左目相机y方向校正映射表，即mapx1
CvMat* m_pcvMapY_Right;				//右目相机y方向校正映射表，即mapy1


CvMat* m_pcvProjectionMat;

//重投影坐标/单目标定求重投影误差时用的
CvMat* m_pcvImgPoint_LeftReProjectiom;
CvMat* m_pcvImgPoint_RightReProjectiom;

char m_FileName[MAX_PATH];			//文件名

int m_LMarkers;						//靶标上任意相邻靶标间圆心距

//行列
int m_iRowCount;					//标定板上的行标记点个数
int m_iColCount;					//标定板上的列标记点个数
int m_iPointCounts;					//标定板上的标记点个数
int m_iImgCounts;					//标定图像个数

int TestNum;						//保存共处理多少张测试图片
int m_iLeftThreshold;				//左图像的阈值
int m_iRightThreshold;				//右图像阈值

IplImage* pTemp;					//临时变量

char m_aLeftReadDirectory[MAX_PATH]; //左图像读取文件地址
char m_aLeftSaveDirectory[MAX_PATH]; //左图像存储文件地址
char m_aRightReadDirectory[MAX_PATH];//右图像读取文件地址
char m_aRightSaveDirectory[MAX_PATH];//右图像存储文件地址
char m_aTestReadDirectory[MAX_PATH]; //测试图像读取文件地址
char m_aTestSaveDirectory[MAX_PATH]; //测试图像存储文件地址

//重投影误差值
float m_fError_Left;
float m_fError_Right;

//未用到的变量
float* pfObjectPoints;				 // 角点物理坐标指针

//3.1 初始化相机内外参数
void InitCameraParams();

//寻找圆心，最小二乘法
void CalCircleCentrePoint(IplImage* &pImage, vector<KeyPoint>&Vct_points,int ipointscount,int iChesscount);

//椭圆拟合
void CalElipseCenterPoint(IplImage* &pImage, vector<KeyPoint>&Vct_points,int ipointscount,int iChesscount);

//计算所有圆点阵列所有的点
BOOL CalcPonitDistance(ST_UpLeft_POINTS UpLeftPonits,vector<ST_CIRCLE_POINT>&Vct_Points );

//获得靶标最大圆心坐标
void GetOrignPoints(vector<ST_CIRCLE_POINT>&Vct_points,int chessboardcount,PST_ORGIN_POINTS &pst_orginPoints);

//得到距原点最近的两个点
double GetMinDistance(vector<ST_CIRCLE_POINT>&Vct_points);

//获得棋盘格原点
BOOL GetUpLeftPoints(vector<ST_CIRCLE_POINT>&Vct_points,PST_ORGIN_POINTS& pst_OringPonit,PST_UpLeft_POINTS& pst_UpLeftPoints,int chessboardcount);

//需找标记点点圆心
BOOL FindCircleCenter(IplImage* img,int minThreshold,int MaxThreshold,int shresholdStep, int minArea,int MaxArea, float minConvexity,float minInertiaRatio,float minCircularity,vector<KeyPoint>&Vct_Points,int iPointsCounts,int ichesscount,int flag);

//对一副图像的所有特征点进行排序
void SortPoints(VEC_POINTS&Vec_points,vector<ST_CIRCLE_POINT>&Circle_Ponits ,PST_UpLeft_POINTS Pst_UpLeft_points,int Chessnum,int Pointsnum);

//计算重投影坐标误差
BOOL CalcReprjection(CvMat*ObjectPoints,CvMat*IntrinsicMat,CvMat*ReprjectMat,CvMat*RotationVec,CvMat*TranslationVec,int iIPointsCounts,int iIiamgeCounts);

//计算重投影误差返回一个平均误差
float CalcReprjectionError(CvMat*Objectpoints,CvMat*ReprjectMat);


void  cvCalibrateCamera2_H( const CvMat* objectPoints,
	const CvMat* imagePoints, const CvMat* npoints,
	CvSize imageSize, CvMat* H,int cnt,
	int flags );

void cvInitIntrinsicParams2D_H( const CvMat* objectPoints, const CvMat* imagePoints, const CvMat* npoints,CvSize imageSize, CvMat* cameraMatrix,
	double aspectRatio,int cnt );


BOOL  ProcessPicture(char picreadaddress_L[],char picsaveaddr_L[],int ithreshold_L, int pictcount,CvMat * objpoints_L, CvMat * imgpoints_L,CvMat * everypicpointcnt_L,CvSize & imgsize_L,
	char picreadaddress_R[],char picsaveaddr_R[],int ithreshold_R, CvMat * objpoints_R, CvMat * imgpoints_R,CvMat * everypicpointcnt_R ,CvSize & imgsize_R);


double cvStereoCalibrate1( const CvMat* _objectPoints1,const CvMat* _objectPoints2, const CvMat* _imagePoints1,
	const CvMat* _imagePoints2, const CvMat* _npoints,
	CvMat* _cameraMatrix1, CvMat* _distCoeffs1,
	CvMat* _cameraMatrix2, CvMat* _distCoeffs2,
	CvSize imageSize, CvMat* matR, CvMat* matT,
	CvMat* matE, CvMat* matF,
	CvTermCriteria termCrit,
	int flags );

static int dbCmp( const void* _a, const void* _b );

void Get3DPos(double* pdbQ, double* pdb2D, double *pdbSolution);


CProgressCtrl m_CProgress;
CListCtrl m_CShowError;
afx_msg void OnBnClickedBtnDcali();
afx_msg void OnBnClickedBtnGetp();
};

#ifndef _DEBUG  // CalibrateView.cpp 中的调试版本
inline CCalibrateDoc* CCalibrateView::GetDocument() const
   { return reinterpret_cast<CCalibrateDoc*>(m_pDocument); }
#endif

