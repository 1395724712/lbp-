#include<stdio.h>
#include<vector>
#include<iostream>
#include<opencv/highgui.h>
#include<math.h>
#include<algorithm>
using namespace std;
using namespace cv;
#pragma once

//δ֪��
const double Beta = 10;
const double Alpha = 10;
const int Iteration_nums = 2;
//�洢�����Ϣ,��ά��ߺͶ�ά���궼���ϵ�һ��
struct Pointpair {
	int stereoleftx, stereolefty, stereoleftz;
	int stereorightx, stereorighty, stereorightz;
	int flatleftx, flatlefty;
	int flatrightx, flatrighty;
};
//�洢��ǩ��Ϣ
struct Label {
	//3*len
	vector<Pointpair> pointpairs;

	int len;
};

//��ͶӰ����ͶӰ��ͼƬΪ���ݣ���¼ÿ�����ص�����Ҫ����Ϣ
struct projector {
	//������λ��Ϣ
	const vector<double> abs_phase;
	//�������ض�Ӧ�ĵ����Ϣ h*w
	const vector<Label> labs;
	//ͶӰ��ͼ��߶�
	//const unsigned int h, w;
};

struct msg_inf {
	//�����Ŀɼ��������msg��Ҫ������¼
	double Left, Right, Down, Up;
};

//��ǩ��Ϣ
enum TASK {initial,updata_msg};
enum DIRECTION { LEFT, RIGHT, UP, DOWN, DATA };

class LBPStereoMatch
{
public:
	LBPStereoMatch(int high,int width,int pic_high,int pic_hidth,projector pro);
	~LBPStereoMatch();
	//Ŀ�꺯��,���ܣ�lbp���֣�������־�����λͼ
	void PartPhaseUnwrapping();
	
	
	//չʾ��λ
	Mat Showphase(DIRECTION dir);

private:
	const int h, w;
	const int pic_h, pic_w;//ͼ��ĳߴ�
	const projector pro;//������Դ�����ͷţ����ﲻʹ������
	vector<int> bel;//���ű��
	Mat L_absphase;//��ͼ������λ
	Mat R_absphase;//��ͼ������λ
	//������λ
	void Computephase();
	//��ʼ������
	void Msg_Recompute(vector<vector<msg_inf> >&, TASK task, vector<vector<msg_inf> >&);
	//���ܺ���,��msg_recompute�з��븴�ӵļ��㲿��
	void Get_Msg(vector<vector<msg_inf> >&, vector<vector<msg_inf> >msg, const DIRECTION, int x, int y, const TASK);
};


