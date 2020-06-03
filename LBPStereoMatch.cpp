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
{//�鷳,��Ҫ�����������鲢��ʵ�ֺ��Ż����ݽṹ
	//��һ���ֵ��������TASK���У���datacost����Ŀ����г�ʼ���������Ը���msgΪĿ��
	//��lbp��ʼ��ʱ���������и�����msg���̵�ǰ�벿��
	//����msg


	//��������ͶӰ����һ����a����n��label��,Ҫ�����õ�һ��nά���ݣ�ÿһά��Ӧ����������õ�������ֵ
	//by wh
	
	vector<vector<msg_inf> >temp_target(h*w);

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{//����Ϣ
			if (i != 0)	Get_Msg(temp_target, msg, UP, i, j, task);
			if (i != h-1)	Get_Msg(temp_target, msg, DOWN, i, j, task);
			if (j != 0)	Get_Msg(temp_target, msg, LEFT, i, j, task);
			if (j != w-1)	Get_Msg(temp_target, msg, RIGHT, i, j, task);


		}
	}

	//����Ŀ����и�ֵ
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
//	//����msg
//	double npqU, npqD, npqR, npqL;
//	//��ô��ݵ���Ϣ
//	vector<vector<double> > msg_from_U(h*w), msg_from_D(h*w), msg_from_R(h*w), msg_from_L(h*w);
//
//	for (int i = 0; i < h; i++)
//	{
//		for (int j = 0; j < w; j++)
//		{
//			if (i == 0)
//			{
//				//for msg_from_U
//				//��һ��û������
//				int len = msgD
//				msg_from_U[i*w + j] ;
//			}
//		}
//	}
//
//	//����npqU
//
//}

void LBPStereoMatch::PartPhaseUnwrapping()
{
	//��ʼ�����ۣ�Ϊmsgȷ����С
	vector<vector<msg_inf> > datacost;
	vector<vector<msg_inf> > placeholder;
	//ע�������placeholderû�����壬����������ռλ��
	Msg_Recompute(datacost,initial,placeholder);

	//��ʼ��msg
	//����������������������ķ����Բ�������
	//����һάֻ�������й����еõ�ֵ
	//by wh
	
	vector<vector<msg_inf> > msg(datacost);


	for (int i = 0; i < Iteration_nums; i++)
	{//������Ϣ
		//����msg
		Msg_Recompute(datacost, updata_msg, msg);
	}

	//�������С����label
	//���ݣ�msg֮����С�ı���Ŀ��
	//�Ե�һ��python����Ϊ����
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

	//����lbp������ɡ�
}

void LBPStereoMatch::Get_Msg(vector<vector<msg_inf>>& target,vector<vector<msg_inf> >msg ,const DIRECTION dir, int i, int j,const TASK task)
{
	//�ⲿ�ָ���Ŀ�����,����msg�����ǳ�ʼ��datacost
	int source_h, source_w;
	//�ڸ���msgʱ��Ӧ���ų�Ŀ������Դ���Ӱ��
	vector<double> exclude;
	
	
	switch (dir)
	{//���ݷ���ȷ��������Դ
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

	//msg����Ϣ
	vector<double> temp_3(pro.labs[source_h*w + source_w].len,0);
	if (task == updata_msg)
	{
		for (int i = 0; i < pro.labs[source_h*w+source_w].len; i++)
		{
			temp_3[i] = msg[source_h*w + source_w][i].Down + msg[source_h*w + source_w][i].Up +
				msg[source_h*w + source_w][i].Right + msg[source_h*w + source_w][i].Left - exclude[i];
		}
	}

	//����target
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
	//����������ͼ��ȴ�Ĳ�����λͼ
	vector<double> L_temp_phase(pic_w*pic_h, 100);
	vector<double> R_temp_phase(pic_w*pic_h, 100);

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			//Ѱ�����label
			int bestlab = bel[i*w + j];
			//�������label��Ӧ��ͼ������
			int L_x = pro.labs[i*w + j].pointpairs[bestlab].flatleftx;
			int L_y = pro.labs[i*w + j].pointpairs[bestlab].flatlefty;
			//��������Ϊ���־�����λͼ��ֵ
			L_temp_phase[L_x*pic_h + L_y] = pro.abs_phase[i*w + j];

			//�������label��Ӧ��ͼ������
			int R_x = pro.labs[i*w + j].pointpairs[bestlab].flatrightx;
			int R_y = pro.labs[i*w + j].pointpairs[bestlab].flatrighty;
			//��������Ϊ���־�����λͼ��ֵ
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
		cerr << "ֻ����������ͼ" << endl;
		break;
	}
}