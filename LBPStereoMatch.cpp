#include "LBPStereoMatch.h"



//LBPStereoMatch::LBPStereoMatch()
//{
//}


LBPStereoMatch::LBPStereoMatch(int high, int width, int pic_high, int pic_width, projector pro)
	:h(high),
	w(width),
	pic_w(pic_width),
	pic_h(pic_high),
	pro(pro)
{
}

LBPStereoMatch::~LBPStereoMatch()
{
}

void LBPStereoMatch::Msg_Recompute(vector<vector<msg_inf> >&datacost,TASK task,vector<vector<msg_inf> >&msg)//by wh
{//麻烦,需要后续处理，建议并行实现和优化数据结构
	//这一部分的任务根据TASK进行，以datacost处理目标进行初始化，还是以更新msg为目标
	//在lbp初始化时计算论文中给出的msg方程的前半部分
	//更新msg


	//即：对于投影仪中一像素a（有n个label）,要求计算得到一个n维数据，每一维对应它从四邻域得到的数据值
	//by wh
	
	vector<vector<msg_inf> >temp_target(h*w);

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{//简化信息
			if (i != 0)	Get_Msg(temp_target, msg, UP, i, j, task);
			if (i != h-1)	Get_Msg(temp_target, msg, DOWN, i, j, task);
			if (j != 0)	Get_Msg(temp_target, msg, LEFT, i, j, task);
			if (j != w-1)	Get_Msg(temp_target, msg, RIGHT, i, j, task);


		}
	}

	//根据目标进行赋值
	if (task==initial)
	{
		datacost.assign(temp_target.begin(), temp_target.end());
	}
	else
	{
		msg.assign(temp_target.begin(), temp_target.end());
	}
}

//void LBPStereoMatch::Update_msg()
//{
//	//更新msg
//	double npqU, npqD, npqR, npqL;
//	//获得传递的信息
//	vector<vector<double> > msg_from_U(h*w), msg_from_D(h*w), msg_from_R(h*w), msg_from_L(h*w);
//
//	for (int i = 0; i < h; i++)
//	{
//		for (int j = 0; j < w; j++)
//		{
//			if (i == 0)
//			{
//				//for msg_from_U
//				//第一层没有数据
//				int len = msgD
//				msg_from_U[i*w + j] ;
//			}
//		}
//	}
//
//	//更新npqU
//
//}

void LBPStereoMatch::PartPhaseUnwrapping()
{
	//初始化代价，为msg确定大小
	vector<vector<msg_inf> > datacost;
	vector<vector<msg_inf> > placeholder;
	//注意这里的placeholder没有意义，仅仅是用于占位的
	Msg_Recompute(datacost,initial,placeholder);

	//初始化msg
	//这边我想用容器的所带来的方便性不及数组
	//但有一维只能在运行过程中得到值
	//by wh
	
	vector<vector<msg_inf> > msg(datacost);


	for (int i = 0; i < Iteration_nums; i++)
	{//迭代信息
		//更新msg
		Msg_Recompute(datacost, updata_msg, msg);
	}

	//获得所最小代价label
	//依据：msg之和最小的便是目标
	//以第一份python代码为依据
	//by wh
	vector<int> bel(h*w);
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			double min_temp = 0xffffffff;
			double temp;
			for (int n = 0; n < msg[i*w+h].size(); n++)
			{
				temp = msg[i*w + h][n].Down + msg[i*w + h][n].Left + msg[i*w + h][n].Right + msg[i*w + h][n].Up;
				if (temp<min_temp)
				{
					min_temp = temp;
					bel[i*w + h] = n;
				}
			}
		}
	}

	Computephase();

	//至此lbp部分完成。
}

void LBPStereoMatch::Get_Msg(vector<vector<msg_inf>>& target,vector<vector<msg_inf> >msg ,const DIRECTION dir, int i, int j,const TASK task)
{
	//这部分根据目标进行,更新msg或者是初始化datacost
	int source_h, source_w;
	//在更新msg时，应当排除目标点对来源点的影响
	vector<double> exclude;
	
	
	switch (dir)
	{//根据方向确定变量来源
	case UP:
		source_h = i + 1;
		source_w = j;

		if (task == updata_msg)
		{
			for (int i = 0; i < pro.labs[source_h*w+j].len; i++)
			{
				exclude.push_back(msg[source_h*w + j][i].Down);
			}
		}

		break;
	case DOWN:
		source_h = i - 1;
		source_w = j;

		if (task == updata_msg)
		{
			for (int i = 0; i < pro.labs[source_h*w + j].len; i++)
			{
				exclude.push_back(msg[source_h*w + j][i].Up);
			}
		}
		break;
	case LEFT:
		source_h = i ;
		source_w = j-1;

		if (task == updata_msg)
		{
			for (int i = 0; i < pro.labs[source_h*w + source_w].len; i++)
			{
				exclude.push_back(msg[source_h*w + source_w][i].Right);
			}
		}
		break;
	case RIGHT:
		source_h = i;
		source_w = j + 1;

		if (task == updata_msg)
		{
			for (int i = 0; i < pro.labs[source_h*w + source_w].len; i++)
			{
				exclude.push_back(msg[source_h*w + source_w][i].Left);
			}
		}
		break;
	default:
		break;
	}

	//Ed
	vector<double> temp_1(pro.labs[source_h*w + source_w].len, 0);
	for (int k = 0; k < pro.labs[i*w + j].len; k++)
	{
		temp_1[k] = sqrt(pow(pro.labs[source_h*w + source_w].pointpairs[k].stereoleftx - pro.labs[source_h*w + source_w].pointpairs[k].stereorightx, 2) +
			pow(pro.labs[source_h*w + source_w].pointpairs[k].stereolefty - pro.labs[source_h*w + source_w].pointpairs[k].stereorighty, 2) +
			pow(pro.labs[source_h*w + source_w].pointpairs[k].stereoleftz - pro.labs[source_h*w + source_w].pointpairs[k].stereorightz, 2));
	}

	//msg和信息
	vector<double> temp_3(pro.labs[source_h*w + source_w].len,0);
	if (task == updata_msg)
	{
		for (int i = 0; i < pro.labs[source_h*w+source_w].len; i++)
		{
			temp_3[i] = msg[source_h*w + source_w][i].Down + msg[source_h*w + source_w][i].Up +
				msg[source_h*w + source_w][i].Right + msg[source_h*w + source_w][i].Left - exclude[i];
		}
	}

	//更新target
	double temp_2;
	for (int n = 0; n < pro.labs[n*w+j].len; n++)
	{
		double min_temp = 0xffffffff;
		for (int m = 0; m < pro.labs[source_h*w + source_w].len; m++)
		{
			int k;
			temp_2 = min(Alpha, sqrt(pow(pro.labs[i*w + j].pointpairs[n].flatleftx - pro.labs[source_h*w + source_w].pointpairs[m].flatleftx, 2) +
				pow(pro.labs[i*w + j].pointpairs[n].flatlefty - pro.labs[source_h*w + source_w].pointpairs[m].flatlefty, 2))) +
				min(Alpha, sqrt(pow(pro.labs[i*w + j].pointpairs[n].flatrightx - pro.labs[source_h*w + source_w].pointpairs[m].flatrightx, 2) +
					pow(pro.labs[i*w + j].pointpairs[n].flatrighty - pro.labs[source_h*w + source_w].pointpairs[m].flatrighty, 2)));
			min_temp = min(min_temp, temp_2 + temp_3[m] + temp_1[m]);
		}
		if (target[i*w+j].size()<n+1)
		{
			target[i*w + j].push_back({ 0,0,0,0 });
		}
		switch (dir)
		{
		case UP:
			target[i*w + j][n].Up = min_temp;
			break;
		case DOWN:
			target[i*w + j][n].Down = min_temp;
			break;
		case LEFT:
			target[i*w + j][n].Left = min_temp;
			break;
		case RIGHT:
			target[i*w + j][n].Right = min_temp;
			break;
		default:
			break;
		}
	}
}

void LBPStereoMatch::Computephase()
{
	//两个和左右图像等大的部分相位图
	vector<double> L_temp_phase(pic_w*pic_h, 100);
	vector<double> R_temp_phase(pic_w*pic_h, 100);

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			//寻找最佳label
			int bestlab = bel[i*w + j];
			//计算最佳label对应的图像坐标
			int L_x = pro.labs[i*w + j].pointpairs[bestlab].flatleftx;
			int L_y = pro.labs[i*w + j].pointpairs[bestlab].flatlefty;
			//根据坐标为部分绝对相位图赋值
			L_temp_phase[L_x*pic_h + L_y] = pro.abs_phase[i*w + j];

			//计算最佳label对应的图像坐标
			int R_x = pro.labs[i*w + j].pointpairs[bestlab].flatrightx;
			int R_y = pro.labs[i*w + j].pointpairs[bestlab].flatrighty;
			//根据坐标为部分绝对相位图赋值
			R_temp_phase[R_x*pic_h + R_y] = pro.abs_phase[i*w + j];

		}
	}
	L_absphase.assign(L_temp_phase.begin(), L_temp_phase.end());
	R_absphase.assign(R_temp_phase.begin(), R_temp_phase.end());

}


vector<double> LBPStereoMatch::Showphase(DIRECTION dir)
{
	switch (dir)
	{
	case LEFT:
		return L_absphase;
		break;
	case RIGHT:
		return R_absphase;
		break;
	default:
		cerr << "只有左右两幅图" << endl;
		break;
	}
}