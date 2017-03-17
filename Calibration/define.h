#pragma once

#include "stdafx.h"
#include <vector>

using namespace std;

#ifndef CALIBRATION
#define CALIBRATION

#define  PI 3.141592

//圆心坐标
typedef struct CIRCLE_POINT
{
	double x;
	double y;
	double size;
	double distance;
	float  dex;
	float  dey;

}ST_CIRCLE_POINT,*PST_CIRLE_POINT;

//存储圆心的容器的容器
typedef vector<vector<CIRCLE_POINT>> VEC_POINTS;

//靶标左上原点
typedef struct UpLeft_POINTS
{
	double x;
	double y;
	double size;
}ST_UpLeft_POINTS,*PST_UpLeft_POINTS;

//靶标标记点
typedef struct ORGIN_POINTS
{
	double x;
	double y;
	double size;
}ST_ORGIN_POINTS,*PST_ORGIN_POINTS;

//棋盘格物理坐标
typedef struct OBJECT_POINTS
{
	double x;
	double y;
	double z;
}ST_OBJECT_POINTS,*PST_OBJECT_POINTS;


#endif
