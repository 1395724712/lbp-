#include<stdio.h>
#include<vector>
#include<iostream>
#include<opencv/highgui.h>
#include<math.h>
#include<algorithm>
using namespace std;
using namespace cv;
#pragma once

//未知量
const double Beta = 10;
const double Alpha = 10;
const int Iteration_nums = 2;
//存储点对信息,三维左边和二维坐标都整合到一起
struct Pointpair {
	int stereoleftx, stereolefty, stereoleftz;
	int stereorightx, stereorighty, stereorightz;
	int flatleftx, flatlefty;
	int flatrightx, flatrighty;
};
//存储标签信息
struct Label {
	//3*len
	vector<Pointpair> pointpairs;

	int len;
};

//以投影仪所投影的图片为依据，记录每个像素点所需要的信息
struct projector {
	//绝对相位信息
	const vector<double> abs_phase;
	//所有像素对应的点对信息 h*w
	const vector<Label> labs;
	//投影仪图像尺度
	//const unsigned int h, w;
};

struct msg_inf {
	//从论文可见，这里的msg需要独立记录
	double Left, Right, Down, Up;
};

//标签信息
enum TASK {initial,updata_msg};
enum DIRECTION { LEFT, RIGHT, UP, DOWN, DATA };

class LBPStereoMatch
{
public:
	LBPStereoMatch(int high,int width,int pic_high,int pic_hidth,projector pro);
	~LBPStereoMatch();
	//目标函数,功能：lbp部分，输出部分绝对相位图
	void PartPhaseUnwrapping();
	
	
	//展示相位
	Mat Showphase(DIRECTION dir);

private:
	const int h, w;
	const int pic_h, pic_w;//图像的尺寸
	const projector pro;//避免资源错误释放，这里不使用引用
	vector<int> bel;//最优标记
	Mat L_absphase;//左图绝对相位
	Mat R_absphase;//右图绝对相位
	//计算相位
	void Computephase();
	//初始化代价
	void Msg_Recompute(vector<vector<msg_inf> >&, TASK task, vector<vector<msg_inf> >&);
	//功能函数,从msg_recompute中分离复杂的计算部分
	void Get_Msg(vector<vector<msg_inf> >&, vector<vector<msg_inf> >msg, const DIRECTION, int x, int y, const TASK);
};


