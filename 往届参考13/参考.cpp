#include <stdio.h>
#include <windows.h>
#include <math.h>
#include <graphics.h>
#include <conio.h>

void rw();
void initialize();
void calculate();
void SOR();
void equal();

///////////////////////主函数///////////////////////
void main()
{
	rw();              //读写文件
	initialize();      //空间电场初始化
	calculate();       //空间电场分布计算
	equal();           //绘制等位线
}

///////////////////////变量定义///////////////////////
FILE *fp;                 //读文件指针    
FILE *p;                  //写文件指针
double d;                 //电极厚度
int n;                    //电极总数
double *dz;               //相邻电极之间距离
int *N;                   //相邻电极之间要划分的步长数
double *V;                //电极电位
double r1;                //电极内孔径半径
int M1;                   //电极内孔径半径范围内等步长划分的网格数
double r2;                //从电极内孔边沿到封闭边界处的径向距离
int M2;                   //从电极内孔边沿到封闭边界处的径向距离范围内等步长划分的网格数
double e;                 //迭代控制精度
int NST;                  //输出打印空间电位时网格点间隔数
						  //int INS;                  //轴上电位作等距插值时的步长数
double dV;                //要求扫描搜索等电位线的电位间隔值
int row = 0, col = 0;          //网格点行、列
int *np;                  //电极位置格点
double **zl;              //z轴相邻格点步长
double **rl;              //r轴相邻格点步长
double *z;                //z轴网格点坐标
double *r;                //r轴网格点坐标
double **V1;              //格点电位
double c1, c2, c3, c4, c0;    //系数
double V2;                //迭代后的格点电位
int T;                    //迭代轮次
int Time;                 //迭代次数
int Timeall;              //迭代总次数
double w;                 //迭代因子
double **Vd;              //格点电位残差
double Vds;               //格点电位残差和
double Vdm;               //格点电位最大残差
double Vda;               //格点电位平均残差
double Vda1;              //格点电位平均残差
double Vda2;              //格点电位平均残差
double x;                 //残差比值
double ux;                //中间变量
double wx;                //中间变量
double wm;                //修正后的w
double wm1;               //修正后的w
int m = 1;                  //等电位数量
double *Ve;               //等电位值
double *zp, *rp;           //等电位点坐标
double *Vzp, *Vrp;         //图形窗口中的等电位点坐标（原坐标扩大10倍）

						   ///////////////////////读写文件///////////////////////
void rw()
{
	int i, j;     //循环变量

				  ///////////////读文件///////////////
	fp = fopen("C:\\Users\\axin\\Desktop\\1120131087_樊阿馨\\1120131087_樊阿馨.dat", "r");//打开文件
	if (fp == NULL)
	{
		printf("Cannot open file.\n");         //如果文件出错显示提示信息
		exit(0);                               //调用exit函数终止程序运行
	}
	fscanf(fp, "电极厚度（占1个步长，所有电极厚度相同）=%lfmm；\n", &d);
	fscanf(fp, "电极总数（包括中间电极与荧光屏，但不包括电位为0V的阴极）=%ld；\n", &n);
	fscanf(fp, "相邻电极之间的距离=");
	dz = (double*)calloc(n * sizeof(double), sizeof(double));
	for (i = 0; i<n; i++)
		fscanf(fp, "%lfmm；", &dz[i]);
	fscanf(fp, "\n相邻电极之间要划分的步长数=");
	N = (int*)calloc(n * sizeof(int), sizeof(int));
	for (i = 0; i<n; i++)
		fscanf(fp, "%ld；", &N[i]);
	fscanf(fp, "\n电极电位=");
	V = (double*)calloc(n * sizeof(double), sizeof(double));
	for (i = 0; i<n; i++)
		fscanf(fp, "%lfV；", &V[i]);
	fscanf(fp, "\n电极内孔径半径（所有电极内孔径半径相同）=%lfmm；\n", &r1);
	fscanf(fp, "电极内孔径半径范围内等步长划分的网格数=%ld；\n", &M1);
	fscanf(fp, "从电极内孔边沿到封闭边界处的径向距离=%lfmm；\n", &r2);
	fscanf(fp, "从电极内孔边沿到封闭边界处的径向距离范围内的等步长划分的网格数=%ld；\n", &M2);
	fscanf(fp, "迭代控制精度=%lfV；\n", &e);
	fscanf(fp, "输出打印空间电位时网格点间隔数=%ld；\n", &NST);
	//fscanf(fp,"轴上电位作等距插值时的步长数=%ld；\n",&INS);
	fscanf(fp, "要求扫描搜索等电位线的电位间隔值或电位值=");
	Ve = (double*)calloc(100 * sizeof(double), sizeof(double));
	while (fscanf(fp, "%lfV；", &Ve[m]) != EOF)
		m++;
	if (m == 2)
		dV = Ve[1];
	fclose(fp);

	///////////////写文件///////////////
	p = fopen("C:\\Users\\axin\\Desktop\\1120131087_樊阿馨\\1120131087_樊阿馨.res", "w");//打开文件
	if (p == NULL)
	{
		printf("Cannot open file.\n");         //如果文件出错显示提示信息
		exit(0);                               //调用exit函数终止程序运行
	}
	fprintf(p, "\t\t\t\t原始数据\n");
	fprintf(p, "电极厚度（占1个步长，所有电极厚度相同）=%lfmm；\n", d);
	fprintf(p, "电极总数（包括中间电极与荧光屏，但不包括电位为0V的阴极）=%ld；\n", n);
	fprintf(p, "相邻电极之间的距离=");
	for (i = 0; i<n; i++)
		fprintf(p, "%lfmm；", dz[i]);
	fprintf(p, "\n相邻电极之间要划分的步长数=");
	for (i = 0; i<n; i++)
		fprintf(p, "%ld；", N[i]);
	fprintf(p, "\n电极电位=");
	for (i = 0; i<n; i++)
		fprintf(p, "%lfV；", V[i]);
	fprintf(p, "\n电极内孔径半径（所有电极内孔径半径相同）=%lfmm；\n", r1);
	fprintf(p, "电极内孔径半径范围内等步长划分的网格数=%ld；\n", M1);
	fprintf(p, "从电极内孔边沿到封闭边界处的径向距离=%lfmm；\n", r2);
	fprintf(p, "从电极内孔边沿到封闭边界处的径向距离范围内的等步长划分的网格数=%ld；\n", M2);
	fprintf(p, "迭代控制精度=%lfV；\n", e);
	fprintf(p, "输出打印空间电位时网格点间隔数=%ld；\n", NST);
	//fprintf(p,"轴上电位作等距插值时的步长数=%ld；\n",INS);
	fprintf(p, "要求扫描搜索等电位线的电位间隔值或电位值=");
	for (i = 1; i<m; i++)
		fprintf(p, "%lfV；", Ve[i]);
	fprintf(p, "\n\n\t\t\t\t网格点坐标(z、r坐标方向网格步长划分)\n");
}

///////////////////////空间电场初始化///////////////////////
void initialize()
{
	int i, j, k;    //循环变量

					///////////////行列计算///////////////
	row = M1 + M2 + 1;
	for (i = 0; i<n; i++)
		col = col + N[i] + 1;

	///////////////网格点坐标计算///////////////
	np = (int*)calloc((n + 1) * sizeof(int), sizeof(int));  //电极位置格点分配内存
	np[0] = 0;
	for (i = 1; i<n + 1; i++)
	{
		if (i == n)
			np[i] = np[i - 1] + N[i - 1];    //荧光屏格点位置
		else
			np[i] = np[i - 1] + N[i - 1] + 1;  //电极右边界格点位置
	}
	zl = (double**)calloc((row - 1) * sizeof(double*), sizeof(double));  //z轴相邻格点步长分配内存
	for (i = 0; i<row - 1; i++)
		zl[i] = (double*)calloc((col - 1) * sizeof(double), sizeof(double));

	///////////////z轴相邻格点步长计算///////////////
	j = 1;      //电极变量
	for (i = 0; i<col - 1; i++)
	{
		if (((i + 1) != np[j]) || (j == n))
		{
			for (k = 0; k<row - 1; k++)
				zl[k][i] = dz[j - 1] / N[j - 1];
		}
		else
		{
			for (k = 0; k<row - 1; k++)
				zl[k][i] = d;
			j++;
		}
	}

	///////////////输出z方向格点位置及其坐标///////////////
	fprintf(p, "z坐标方向\n");
	for (i = 0; i<col; i++)
		fprintf(p, "%ld\t\t", i);
	z = (double*)calloc(col * sizeof(double), sizeof(double));
	for (i = 0; i<col; i++)
	{
		if (i == 0)
		{
			z[i] = 0;
			fprintf(p, "\n%.2lf\t", z[i]);
		}
		else
		{
			z[i] = z[i - 1] + zl[0][i - 1];
			fprintf(p, "%.2lf\t", z[i]);
		}
	}
	rl = (double**)calloc((row - 1) * sizeof(double*), sizeof(double));  //r轴相邻格点步长分配内存
	for (i = 0; i<row - 1; i++)
		rl[i] = (double*)calloc((col - 1) * sizeof(double), sizeof(double));

	///////////////r轴相邻格点步长计算///////////////
	for (i = 0; i<col - 1; i++)
	{
		for (j = 0; j<row - 1; j++)
		{
			if (j<M1)
				rl[j][i] = r1 / M1;
			else
				rl[j][i] = r2 / M2;
		}
	}

	///////////////输出r方向格点位置及其坐标///////////////
	fprintf(p, "\nr坐标方向\n");
	for (i = 0; i<row; i++)
		fprintf(p, "%ld\t\t", i);
	r = (double*)calloc(row * sizeof(double), sizeof(double));
	for (i = 0; i<row; i++)
	{
		if (i == 0)
		{
			r[i] = 0;
			fprintf(p, "\n%.2lf\t", r[i]);
		}
		else
		{
			r[i] = r[i - 1] + rl[i - 1][0];
			fprintf(p, "%.2lf\t", r[i]);
		}
	}
	fprintf(p, "\n\n");

	///////////////格点电位初始化///////////////
	V1 = (double**)calloc(row * sizeof(double*), sizeof(double));  //格点电位分配内存
	for (i = 0; i<row; i++)
		V1[i] = (double*)calloc(col * sizeof(double), sizeof(double));
	for (i = 0; i<row; i++)  //阴极位置电位初始化赋值
		V1[i][0] = 0;
	for (i = 1; i<n + 1; i++)  //电极位置电位初始化赋值
	{
		if (i == n)
		{
			for (j = 0; j<row; j++)
				V1[j][np[i]] = V[i - 1];
		}
		else
		{
			for (j = M1; j<row; j++)
			{
				V1[j][np[i] - 1] = V[i - 1];
				V1[j][np[i]] = V[i - 1];
			}
		}
	}
	j = 1;      //电极变量
	for (k = 1; k<col - 1; k++)  //上边界电极间电位初始化赋值
	{
		if (((k<np[n])&(k>np[n - 1])) || ((k<np[j] - 1)&(k>np[j - 1])))
		{
			if (j == 1)
				V1[row - 1][k] = V1[row - 1][k - 1] + V[j - 1] / N[j - 1];
			else
				V1[row - 1][k] = V1[row - 1][k - 1] + (V[j - 1] - V[j - 2]) / N[j - 1];
		}
		else
		{
			if ((k != np[j])&(k != np[j - 1]))
				j++;
		}
	}
}

//////////////////////////空间电场分布计算/////////////////////////
void calculate()
{
	int i, j, k, rule;
	T = 1;
	Time = 1;
	Timeall = 0;
	w = 1;
	SOR();   //第一轮迭代
	Timeall++;
	fprintf(p, "迭代轮次：%-3ld\t\t迭代次数：%-3ld\t\t每轮迭代因子值:%lf\t\t", T, Timeall, w);
	fprintf(p, "每轮迭代结束时平均残差：%lf\t\t最大残差：%lf\n", Vda, Vdm);
	w = 1.375;
	T++;
	for (Time = 1; Time<13; Time++)
	{
		SOR();   //第二轮迭代
		Timeall++;
		if (Time == 11)
			Vda1 = Vda;
		else if (Time == 12)
			Vda2 = Vda;
	}
	fprintf(p, "迭代轮次：%-3ld\t\t迭代次数：%-3ld\t\t每轮迭代因子值:%lf\t\t", T, Timeall, w);
	fprintf(p, "每轮迭代结束时平均残差：%lf\t\t最大残差：%lf\n", Vda, Vdm);
	x = Vda2 / Vda1;
	ux = (x + w - 1) / (sqrt(x)*w);
	wx = 2 / (1 + sqrt(1 - ux * ux));
	wm = 1.25*wx - 0.5;
	do
	{
		w = wm;
		T++;
		for (Time = 1; Time<13; Time++)
		{
			SOR();   //第三轮迭代
			Timeall++;
			if (Time == 11)
				Vda1 = Vda;
			else if (Time == 12)
				Vda2 = Vda;
		}
		fprintf(p, "迭代轮次：%-3ld\t\t迭代次数：%-3ld\t\t每轮迭代因子值:%lf\t\t", T, Timeall, w);
		fprintf(p, "每轮迭代结束时平均残差：%lf\t\t最大残差：%lf\n", Vda, Vdm);
		wm1 = wm;
		x = Vda2 / Vda1;
		ux = (x + w - 1) / (sqrt(x)*w);
		wx = 2 / (1 + sqrt(1 - ux * ux));
		wm = 1.25*wx - 0.5;
	} while (fabs((wm - wm1) / (2 - wm1)) >= 0.05);
	w = wm;
	do
	{
		T++;
		SOR();   //最佳迭代因子确定后继续迭代
		Timeall++;
		fprintf(p, "迭代轮次：%-3ld\t\t迭代次数：%-3ld\t\t每轮迭代因子值:%lf\t\t", T, Timeall, w);
		fprintf(p, "每轮迭代结束时平均残差：%lf\t\t最大残差：%lf\n", Vda, Vdm);
		rule = 0;
		for (i = 0; i<row; i++)
			for (j = 0; j<col; j++)
			{
				if (Vd[i][j] >= e)
					rule++;
			}
	} while (rule>0);

	//////////////////输出空间电场迭代后结果//////////////////
	fprintf(p, "\n\t\t\t\t网格点上的电位值\n\t\t");
	for (i = 0; i<col; i++)
	{
		fprintf(p, "%ld\t\t\t", i);
		i = i + NST - 1;
	}
	fprintf(p, "\n");
	for (i = row - NST; i >= 0; i--)
	{
		fprintf(p, "%ld\t\t", i);
		for (j = 0; j<col; j++)
		{
			fprintf(p, "%lf\t", V1[i][j]);
			j = j + NST - 1;
		}
		fprintf(p, "\n");
		i = i - NST + 1;
	}
}

////////////////////////SOR算法/////////////////////////
void SOR()
{
	int i, j, k;   //循环变量
	Vd = (double**)calloc(row * sizeof(double*), sizeof(double));  //格点电位残差分配内存
	for (i = 0; i<row; i++)
		Vd[i] = (double*)calloc(col * sizeof(double), sizeof(double));
	for (i = 0; i<row; i++)
	{
		for (j = 0; j<col; j++)
		{
			int judge = 1;   //边界点判断指标
			if ((i != row - 1)&(j != 0)&(j != col - 1))   //边界点判断
			{
				if (i>M1 - 1)
				{
					for (k = 1; k<n; k++)
					{
						judge = 0;
						if ((j == np[k] - 1) || (j == np[k]))
						{
							judge = 1;
							k = n;
						}
					}
				}
				else
					judge = 0;
			}
			if (judge == 0)   //非边界点计算
			{
				c1 = 2 / (zl[i][j - 1] * (zl[i][j - 1] + zl[i][j]));   //系数计算
				c2 = 2 / (zl[i][j] * (zl[i][j - 1] + zl[i][j]));
				if (i == 0)
				{
					c3 = 0;
					c4 = 4 / (rl[i][j] * rl[i][j]);
				}
				else
				{
					c3 = (2 * r[i] - rl[i][j]) / (r[i] * rl[i - 1][j] * (rl[i - 1][j] + rl[i][j]));
					c4 = (2 * r[i] + rl[i - 1][j]) / (r[i] * rl[i][j] * (rl[i - 1][j] + rl[i][j]));
				}
				c0 = c1 + c2 + c3 + c4;
				if (i != 0)
					V2 = (1 - w)*V1[i][j] + w * (c1*V1[i][j - 1] + c2 * V1[i][j + 1] + c3 * V1[i - 1][j] + c4 * V1[i + 1][j]) / c0;
				else
					V2 = (1 - w)*V1[i][j] + w * (c1*V1[i][j - 1] + c2 * V1[i][j + 1] + c4 * V1[i + 1][j]) / c0;
				Vd[i][j] = fabs(V2 - V1[i][j]);
				V1[i][j] = V2;
			}
		}
	}
	Vds = 0;
	Vdm = 0;
	for (i = 0; i<row; i++)
	{
		for (j = 0; j<col; j++)
		{
			Vds = Vds + Vd[i][j];
			if (Vd[i][j]>Vdm)
				Vdm = Vd[i][j];
		}
	}
	Vda = Vds / (row*col);
}

////////////////////////绘制等位线////////////////////////
void equal()
{
	int i, i1, j, j1, k, l, c, b, pan, f, f1, f2, f3, f4, f5, f6, *hi;
	double a = 0, zlm, rlm, zrlm;
	if (m == 2)
	{
		for (m = 1; a<(V[n - 1] - dV); m++)   //计算电位数量
			a = a + dV;
		Ve = (double*)calloc((m - 1) * sizeof(double), sizeof(double));
		for (i = 0; Ve[i]<(V[n - 1] - dV); i++)
			Ve[i + 1] = Ve[i] + dV;         //计算各个电位值
	}

	//////////////////计算等位线坐标//////////////////
	initgraph(100 + int(10 * z[col - 1]), 100 + int(10 * r[row - 1]));   //画面距窗口边缘四周均50，坐标值均扩大10倍
	setbkcolor(BLACK);   //采用黑色背景
	cleardevice();
	rectangle(50, 50, 50 + int(10 * z[col - 1]), 50 + int(10 * r[row - 1]));   //距窗口边缘四周均50画矩形
	for (i = 1; i<n; i++)
		bar(50 + int(10 * z[np[i] - 1]), 50, 50 + int(10 * (z[np[i]])), 50 + int(10 * r2));   //画电极
	zp = (double*)calloc(row*col * sizeof(double), sizeof(double));   //等电位点z坐标值分配内存
	rp = (double*)calloc(row*col * sizeof(double), sizeof(double));   //等电位点r坐标值分配内存
	Vzp = (double*)calloc(row*col * sizeof(double), sizeof(double));  //扩大10倍的等电位点z坐标值分配内存
	Vrp = (double*)calloc(row*col * sizeof(double), sizeof(double));  //扩大10倍的等电位点r坐标值分配内存
	zlm = zl[0][0];
	for (i = 0; i<col - 1; i++)
	{
		if (zlm<zl[0][i])
			zlm = zl[0][i];        //计算z方向最大格点间距
	}
	rlm = rl[0][0];
	if (rlm<rl[row - 2][0])
		rlm = rl[row - 2][0];        //计算r方向最大格点间距
	zrlm = sqrt(zlm*zlm + rlm * rlm);	 //计算对角线方向最大间距
	for (k = 1; k<m; k++)
	{
		f1 = 0;    //记录顶行等电位点数量初始化
		l = 0;     //等电位点数量初始化
		for (i1 = row - 1; i1 >= 0; i1--)
		{
			for (j1 = 0; j1<col - 1; j1++)   //z轴方向线性插值计算等电位点坐标
			{
				if (((V1[i1][j1]<Ve[k])&(Ve[k]<V1[i1][j1 + 1])) || ((V1[i1][j1]>Ve[k])&(Ve[k]>V1[i1][j1 + 1])) || (V1[i1][j1] == Ve[k]))
				{
					rp[l] = r[i1];
					if (rp[l] == r[row - 1])
						f1++;          //顶行标记值记录顶行等电位点数量用于控制等位线绘制
					if (V1[i1][j1] == Ve[k])
						zp[l] = z[j1];
					else
						zp[l] = z[j1] + (Ve[k] - V1[i1][j1])*zl[0][j1] / (V1[i1][j1 + 1] - V1[i1][j1]);
					Vrp[l] = 50 + 10 * (r[row - 1] - rp[l]);  //扩大10倍后的图形坐标 
					Vzp[l] = 50 + 10 * zp[l];             //扩大10倍后的图形坐标 
					l++;
				}
			}
			for (j = 0; j<col; j++)        //r轴方向线性插值计算等电位点坐标
			{
				if (i1>0)
				{
					i = i1;
					if (((V1[i][j]<Ve[k])&(Ve[k]<V1[i - 1][j])) || ((V1[i][j]>Ve[k])&(Ve[k]>V1[i - 1][j])))
					{
						zp[l] = z[j];
						rp[l] = r[i] - (Ve[k] - V1[i][j])*rl[i - 1][0] / (V1[i - 1][j] - V1[i][j]);
						Vrp[l] = 50 + 10 * (r[row - 1] - rp[l]);   //扩大10倍后的图形坐标
						Vzp[l] = 50 + 10 * zp[l];              //扩大10倍后的图形坐标
						l++;
					}
				}
			}
		}

		//////////////////绘制等位线//////////////////
		f = 0;   //连接线段起点位置变量初始化
		c = 0;   //储存连接过的等电位点位置循环变量初始化
		f6 = 0;  //用于记录电位等于电极电位时的电极位置初始化
		hi = (int*)calloc(l * sizeof(int), sizeof(int));   //储存连接过的等电位点位置
		do
		{
			fprintf(p, "\n电位值为%lf的等位线的各点坐标值：\n", Ve[k]);
			setcolor(20 * k*RED + 10 * k*GREEN + BLUE);    //确定等位线颜色
			for (i = 0; i<l; i++)
			{
				pan = 1;
				for (b = 0; b<c; b++)    //判断是否为已连接过的点
				{
					if (i == hi[b])
					{
						pan = 0;
						break;
					}
				}
				if ((pan == 1)&(rp[i] == r[row - 1]))    //用于有鞍点时找出下一条等位线的起点位置
				{
					f = i;
					break;
				}
			}
			for (i = f; i<l; i++)
			{
				for (j = 0; j<l; j++)
				{
					pan = 1;
					for (b = 0; b<c; b++)    //判断是否为已连接过的点
					{
						if (j == hi[b])
						{
							pan = 0;
							break;
						}
					}
					if ((j != i)&(pan == 1)&(sqrt((Vzp[j] - Vzp[i])*(Vzp[j] - Vzp[i]) + (Vrp[j] - Vrp[i])*(Vrp[j] - Vrp[i]))<12 * zrlm)&(fabs(Vzp[j] - Vzp[i])<12 * zlm)&(fabs(Vrp[j] - Vrp[i])<12 * rlm))
					{    //如果第二点不是第一点 且 不是已连接过的点 且 两点之间距离满足要求
						f3 = 0;
						for (f2 = 0; f2<n; f2++)    //如果等于电极电位则记录电极位置并做标记
						{
							if (Ve[k] == V[f2])
							{
								f3 = 1;    //做标记
								f6 = f2;   //记录电极位置
								break;
							}
						}
						if ((f3 == 1)&((z[np[f6 + 1]] == zp[j]) || (z[np[f6 + 1] - 1] == zp[j]))&(r[row - 1] == rp[j]))
						{    //如果是电极电位且第二点为电极顶部所在位置则从上往下“Z”字型连接电极电位点
							for (f2 = 0; f2<M2; f2++)
							{
								fprintf(p, "(%lf\t,\t%lf\t)\n", rp[i], zp[i]);   //输出第一点坐标
								if (rp[i] == r[row - 1])
									f1--;     //如果第一点在顶行则标记值减1
								line(int(Vzp[i]), int(Vrp[i]), int(Vzp[j]), int(Vrp[j]));
								hi[c] = i;      //记录第一点位置
								c++;
								fprintf(p, "(%lf\t,\t%lf\t)\n", rp[j], zp[j]);   //输出第二点坐标
								if (rp[j] == r[row - 1])
									f1--;	  //如果第二点在顶行则标记值减1			
								for (f4 = 0; f4<l; f4++)
								{
									if ((zp[f4] == zp[i])&(rp[f4] == rp[i] - r2 / M2))
									{
										i = f4;      //寻找下一个第一点
										break;
									}
								}
								line(int(Vzp[j]), int(Vrp[j]), int(Vzp[i]), int(Vrp[i]));
								hi[c] = j;
								c++;
								for (f5 = 0; f5<l; f5++)
								{
									if ((zp[f5] == zp[j])&(rp[f5] == rp[j] - r2 / M2))
									{
										j = f5;      //寻找下一个第二点
										break;
									}
								}
							}
							fprintf(p, "(%lf\t,\t%lf\t)\n", rp[i], zp[i]);  //输出电极底端左侧电位点坐标
							if (rp[i] == r[row - 1])
								f1--;			//如果此点在顶行则标记值减1		
							line(int(Vzp[i]), int(Vrp[i]), int(Vzp[j]), int(Vrp[j]));   //连接电极底端水平线
							hi[c] = i;            //记录电极底端左侧电位点位置
							c++;
							break;
						}
						else if ((f3 == 1)&((z[np[f6 + 1]] == zp[j]) || (z[np[f6 + 1] - 1] == zp[j]))&(r[row - 1] != rp[j]))
						{   //如果是电极电位且第二点不是电极顶部所在位置则从下往上倒“Z”字型连接电极电位点
							fprintf(p, "(%lf\t,\t%lf\t)\n", rp[i], zp[i]);
							if (rp[i] == r[row - 1])
								f1--;
							line(int(Vzp[i]), int(Vrp[i]), int(Vzp[j]), int(Vrp[j]));
							hi[c] = i;
							c++;
							fprintf(p, "(%lf\t,\t%lf\t)\n", rp[j], zp[j]);
							if (rp[j] == r[row - 1])
								f1--;
							for (f4 = 0; f4<l; f4++)
							{
								if ((zp[f4] == zp[np[f6 + 1]])&(rp[f4] == rp[M1]))
								{
									i = f4;
									j = f4 - 1;
									break;
								}
							}
							for (f2 = 0; f2<M2; f2++)
							{
								fprintf(p, "(%lf\t,\t%lf\t)\n", rp[i], zp[i]);
								if (rp[i] == r[row - 1])
									f1--;
								line(int(Vzp[i]), int(Vrp[i]), int(Vzp[j]), int(Vrp[j]));
								hi[c] = i;
								c++;
								fprintf(p, "(%lf\t,\t%lf\t)\n", rp[j], zp[j]);
								if (rp[j] == r[row - 1])
									f1--;
								for (f4 = 0; f4<l; f4++)
								{
									if ((zp[f4] == zp[i])&(rp[f4] == rp[i] + r2 / M2))
									{
										i = f4;
										break;
									}
								}
								line(int(Vzp[j]), int(Vrp[j]), int(Vzp[i]), int(Vrp[i]));
								hi[c] = j;
								c++;
								for (f5 = 0; f5<l; f5++)
								{
									if ((zp[f5] == zp[j])&(rp[f5] == rp[j] + r2 / M2))
									{
										j = f5;
										break;
									}
								}
							}
							fprintf(p, "(%lf\t,\t%lf\t)\n", rp[i], zp[i]);
							if (rp[i] == r[row - 1])
								f1--;
							line(int(Vzp[i]), int(Vrp[i]), int(Vzp[j]), int(Vrp[j]));
							hi[c] = i;
							c++;
							break;
						}
						else
						{   //如果第二点不是电极所在位置则输出第一点坐标并连接第一二点
							fprintf(p, "(%lf\t,\t%lf\t)\n", rp[i], zp[i]);
							if (rp[i] == r[row - 1])
								f1--;
							line(int(Vzp[i]), int(Vrp[i]), int(Vzp[j]), int(Vrp[j]));
							hi[c] = i;
							c++;
							break;
						}
					}
				}
				if ((rp[j] == r[0]) || (rp[j] == r[row - 1]))
				{   //如果第二点位于顶行或底行则输出第二点坐标并跳出循环判断是否有下一条等位线
					fprintf(p, "(%lf\t,\t%lf\t)\n", rp[j], zp[j]);
					if (rp[j] == r[row - 1])
						f1--;
					hi[c] = j;
					c++;
					break;
				}
				else   //如果第二点不位于顶行或底行则将第二点作为下一次循环的第一点
					i = j - 1;
			}
		} while (f1>0);      //如果顶行标记值不为0则继续搜索连接下一条等位线
	}
	fclose(p);
	system("pause");
	closegraph();
}
