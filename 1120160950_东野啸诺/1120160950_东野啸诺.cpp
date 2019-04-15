#include<stdio.h>
#include<graphics.h>
#include<math.h>
#include<stdlib.h>
#include <conio.h>
#include <stdbool.h>

//////////////////1.读入静态变量/////////////////////////////
static int n;					  //电极总数
static double D;				  //电极厚度
static double *dz;				  //相邻电极之间的距离
static int *N;					  //相邻电极之间要划分的步长数
static double *V;				  //电极电位
static double dr1;				  //电极内孔径半径
static int M1;					  //dr1范围内等步长划分的网格数
static double dr2;				  //从电极内孔边沿到封闭边界处的径向距离
static int M2;					  //dr2范围内等步长划分的网格数
static double e;				  //迭代控制精度e
static int NST;					 //输出打印空间电位时网格点间隔数
static int INS;                  //轴上电位作等距插值时的步长数
static double dV;				 //要求扫描搜索等电位线的电位间隔值
static double *Ve;               //要求扫描搜索等电位线的电位值
static int x;					 //电位个数
int M;							 //x的最大值
void ReadFile()									//读入数据
{
	static FILE *infile;						//设置文件指针
	fopen_s(&infile, "X:\\CADtest\\1120160950\\1120160950_东野啸诺ad.dat", "r");
	if (infile == NULL)
	{
		printf("Opps, Cannot Open File!\n");
		return;
	}
	fscanf_s(infile, "电极总数n=%d;\n", &n, sizeof(int));
	fscanf_s(infile, "电极厚度D=%lfmm;\n", &D, sizeof(double));
	fscanf_s(infile, "相邻电极之间的距离dz=", sizeof(double));
	dz = (double*)calloc(n * sizeof(double), sizeof(double));
	for (int i = 0; i < n; i++)
		fscanf_s(infile, "%lfmm;", &dz[i], sizeof(double));
	fscanf_s(infile, "\n相邻电极之间要划分的步长数N=", sizeof(int));
	N = (int*)calloc(n * sizeof(int), sizeof(int));
	for (int i = 0; i < n; i++)
		fscanf_s(infile, "%d;", &N[i], sizeof(int));
	fscanf_s(infile, "\n电极电位V=", sizeof(double));
	V = (double*)calloc(n * sizeof(double), sizeof(double));
	for (int i = 0; i < n; i++)
		fscanf_s(infile, "%lfV;", &V[i], sizeof(double));
	fscanf_s(infile, "\n电极内孔径半径dr1=%lfmm;\n", &dr1, sizeof(double));
	fscanf_s(infile, "dr1范围内等步长划分的网格数M1=%d;\n", &M1, sizeof(int));
	fscanf_s(infile, "从电极内孔边沿到封闭边界处的径向距离dr2=%lfmm;\n", &dr2, sizeof(double));
	fscanf_s(infile, "dr2范围内等步长划分的网格数M2=%d;\n", &M2, sizeof(int));
	fscanf_s(infile, "迭代控制精度e=%lfV;\n", &e, sizeof(double));
	fscanf_s(infile, "输出打印空间电位时网格点间隔数NST=%d;\n", &NST, sizeof(int));
	fscanf_s(infile, "轴上电位作等距插值时的步长数INS=%d;\n", &INS, sizeof(int));
	fscanf_s(infile, "等电位线的电位间隔或电位值Ve=");
	Ve = (double*)calloc(100 * sizeof(double), sizeof(double));
	while (fscanf_s(infile, "%lfV;", &Ve[x], sizeof(double)) != EOF)
		x += 1;
	if (x == 1)
	{
		dV = Ve[0];
		M = 100 / dV;
		for (int i = 0; i < M; i++)
		{
			Ve[i] = (i + 1)*Ve[0];
		}
	}
	fclose(infile);
}

////////////////////2.检查输入正确性////////////////////////////
static FILE *outfile;			//定义输出文件指针
void CheckInput()													//检查数据是否正确读入
{
	fopen_s(&outfile, "X:\\CADtest\\1120160950\\1120160950_东野啸诺.res", "w");       //生成res文件
	if (outfile == NULL)
	{
		printf("Opps, Cannot Write This File!\n");
		return;
	}
	//
	fprintf(outfile, "\\\\\\\\\\\\\\\\\\\\\\\\\\\\读入数据检查/////////////////////////////\n");
	fprintf(outfile, "\n");
	fprintf(outfile, "                        .::::.\n");
	fprintf(outfile, "                      .::::::::.\n");
	fprintf(outfile, "                     :::::::::::\n");
	fprintf(outfile, "                  ..:::::::::::'\n");
	fprintf(outfile, "               '::::::::::::'\n");
	fprintf(outfile, "              .::::::::::::::\n");
	fprintf(outfile, "            '::::::::::::::..\n");
	fprintf(outfile, "            :::: ..::::::::::::.\n");
	fprintf(outfile, "               ``::::::::::::::::\n");
	fprintf(outfile, "                ::::``:::::::::'        .:::.\n");
	fprintf(outfile, "               ::::'   ':::::'       .::::::::.\n");
	fprintf(outfile, "             .::::'      ::::     .:::::::'::::.\n");
	fprintf(outfile, "            .:::'       :::::  .:::::::::' ':::::.\n");
	fprintf(outfile, "           .::'        :::::.:::::::::'      ':::::.\n");
	fprintf(outfile, "          .::'         ::::::::::::::'         ``::::.\n");
	fprintf(outfile, "      ...:::           ::::::::::::'              ``::.\n");
	fprintf(outfile, "     ```` ':.          ':::::::::'                  ::::..\n");
	fprintf(outfile, "                        '.:::::'                    ':'````..\n");
	fprintf(outfile, "=====================================================================\n");
	fprintf(outfile, "                     美女保佑 永无BUG\n");
	fprintf(outfile, "电极总数n=%ld;\n", n);
	fprintf(outfile, "电极厚度D=%lfmm;\n", D);
	fprintf(outfile, "相邻电极之间的距离dz=");
	for (int i = 0; i < n; i++)
		fprintf(outfile, "%lfmm;", dz[i]);
	fprintf(outfile, "\n相邻电极之间要划分的步长数N=");
	for (int i = 0; i < n; i++)
		fprintf(outfile, "%ld;", N[i]);
	fprintf(outfile, "\n电极电位V=");
	for (int i = 0; i < n; i++)
		fprintf(outfile, "%lfV;", V[i]);
	fprintf(outfile, "\n电极内孔径半径r1=%lfmm;\n", dr1);
	fprintf(outfile, "dr1范围内等步长划分的网格数M1%ld;\n", M1);
	fprintf(outfile, "从电极内孔边沿到封闭边界处的径向距离dr2=%lfmm;\n", dr2);
	fprintf(outfile, "dr2范围内等步长划分的网格数M2=%ld;\n", M2);
	fprintf(outfile, "迭代控制精度e=%lfV;\n", e);
	fprintf(outfile, "输出打印空间电位时网格点间隔数NST=%ld;\n", NST);
	fprintf(outfile, "轴上电位作等距插值时的步长数INS=%ld;\n", INS);
	fprintf(outfile, "等电位线的电位间隔或电位值Ve=");
	for (int i = 0; i < x; i++)
		fprintf(outfile, "%lfV;", Ve[i]);
	fprintf(outfile, "\n\t^\t^\n");
	fprintf(outfile, "╭┘└┘└╮              牛\n");
	fprintf(outfile, "└┐．．┌┘────╮    头\n");
	fprintf(outfile, "╭┴──┤          ├╮  酋\n");
	fprintf(outfile, "│ｏ　ｏ│          │|● 长\n");
	fprintf(outfile, "╰─┬─╯          │|   阿\n");
	fprintf(outfile, "    || ||________|| ||    利\n");
	fprintf(outfile, "    || ||        || ||    斯\n");
	fprintf(outfile, "   |_||_|       |_||_|    塔\n");
	fprintf(outfile, "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\检查完毕，读入文件正确/////////////////////////////\n");
};

////////////////////3.初始化电场空间////////////////////////////
static int Vertical;
static int Horizontal;							//定义行列数
static double **E;							    //定义电场分布指针
												//注明：定义为静态变量的好处开辟了新的内存空间，避免了循环赋值的过程
void InitializeField()
{
	for (int i = 0; i < n; i++)			//纵向
	{
		Horizontal += N[i] + 1;			//按照步长数初始化
	}
	Vertical = M1 + M2 + 1;				  //横向
	fprintf(outfile, "\n网格划分成功。 %d 行%d 列；\n", Vertical, Horizontal);
	E = (double**)calloc(Vertical * sizeof(double*), sizeof(double));		 //电场大小
	for (int i = 0; i < Vertical; i++)
		E[i] = (double*)calloc(Horizontal * sizeof(double), sizeof(double));

	//初值0
	for (int i = 0; i < Vertical; i++)
	{
		for (int j = 0; j < Horizontal; j++)
		{
			E[i][j] = 0;
		}
	}
	//为电极左右电压赋值
	static int k;
	for (int i = 0; i <= M2; i++)
	{
		for (int j = 0; j < n; j++)
		{
			E[i][k + N[j]] = V[j];
			E[i][k + N[j] + 1] = V[j];
			k += N[j] + 1;
		}
	}

	//电极之间电压
	for (int i = 0; i < M2; i++)
	{
		static int k;
		static int l;										//循环变量
		for (int j = k + 1; j < k + N[l]; j++)
		{
			E[i][j] = E[i][j - 1] + V[l] / N[l];		//V[]/N[]表现了相邻电极之间插值相同	          
		}
		k += N[l] + 1;
		for (l = 1; l < n; l++)
		{
			for (int j = k + 1; j<k + N[l]; j++)
			{
				E[i][j] = E[i][j - 1] + (V[l] - V[l - 1]) / N[l];
			}
			k += N[l] + 1;
		}
		for (int j = 0; j < Horizontal; ++j)				//靠近电极边缘赋值
		{
			E[i + 1][j] = E[i][j];
		}	
	}

	//荧光屏赋值
	for (int i = 0; i < Vertical; i++)
	{
		E[i][Horizontal - 1] = V[n - 1];
	}

	//显示赋值结果
	for (int i = 0; i < Vertical; i++)
	{
		for (int j = 0; j < Horizontal; j++)
		{
			printf(" %lf ", E[i][j]);
		}
		printf("\n已经成功赋值，电场初始化完成\n");
	}
}

///////////////////////4.初始化网格////////////////////////////
double *GridLenth;
double *GridWidth;				//格点长宽
double *X;
double *Y;						//格点坐标
int *Rank;

void InitializeGrid()
{
	GridLenth = (double*)calloc((Horizontal - 1) * sizeof(double), sizeof(double));
	GridWidth = (double*)calloc((Vertical - 1) * sizeof(double), sizeof(double));
	X = (double*)calloc(Horizontal, sizeof(double));
	Y = (double*)calloc(Vertical, sizeof(double));		//定义大小
	int k = 0;
	for (int i = 0; i < n; i++)						//纵向计算
	{
		for (int j = k; j < (k + N[i]); j++)
		{
			GridLenth[j] = dz[i] / N[i];
		}
		k = k + N[i];
		if (k >= Horizontal - 1)
		{
			break;
		}
		GridLenth[k] = D;
		k++;
	}

	for (int i = 0; i < Vertical - 1; i++)			//横向计算
	{
		if (i < M2)
			GridWidth[i] = dr2 / M2;
		else
			GridWidth[i] = dr1 / M1;
	}
	for (int i = 0; i < Horizontal; i++)			//纵向赋值
	{
		if (i == 0)
			X[i] = 0;
		else
			X[i] = X[i - 1] + GridLenth[i - 1];
	}
	for (int i = 0; i < Vertical; i++)				//横向赋值
	{
		if (i == 0)
			Y[i] = 0;
		else
			Y[i] = Y[i - 1] + GridWidth[i - 1];
	}

	fprintf(outfile, "\n格点坐标:\n");
	fprintf(outfile, "格点X坐标:\n");
	for (int i = 0; i < Horizontal; i++)
		fprintf(outfile, "%ld\t", i);
	fprintf(outfile, "\n");
	for (int i = 0; i < Horizontal; i++)
		fprintf(outfile, "%lf\t", X[i]);
	fprintf(outfile, "\n格点Y坐标\n");
	for (int i = 0; i < Vertical; i++)
		fprintf(outfile, "%ld\t", i);
	fprintf(outfile, "\n");
	for (int i = 0; i < Vertical; i++)
		fprintf(outfile, "%lf\t", Y[i]);
}

////////////////////////5.主函数////////////////////////////////
void ComputeElectricField();
void EqupotentialLine();
void SOR();

void main()														  //主函数
{	
	ReadFile();
	CheckInput();
	InitializeField();
	InitializeGrid();
	ComputeElectricField();
	EqupotentialLine();
	system("pause");
}

////////////////////////6.计算电场////////////////////////////////
/////SOR初始化参数/////
double **Res;					//残差
double ResSum;					//残差和
double ResMax;					//最大残差
double Va1;
double Va2;
double ResAverage;				//平均残差
double c1, c2, c3, c4, c0;
double ω;						//超张弛迭代因子
double U;						//中心点电位

//SOR
void SOR()
{
	Res = (double**)calloc(Vertical * sizeof(double*), sizeof(double));	 //残差内存大小
	for (int i = 0; i < Vertical; i++)
		Res[i] = (double*)calloc(Horizontal * sizeof(double), sizeof(double));

	for (int i = 0; i < Vertical; i++)
	{
		for (int j = 0; j < Horizontal; j++)
		{
			bool Flag = TRUE;
			int K = 0;
			if ((i != 0) && (j != 0) && (j != Horizontal - 1))        //判断边界点
			{
				Flag = FALSE;
				for (int k = 0; k < n; k++)
				{
					if (i <= M2)
					{
						if (j >= (K + N[k]) && j <= (K + N[k] + 1))
						{
							Flag = TRUE;
							break;
						}
					}
					K += N[k] + 1;
				}
			}
			if (Flag == FALSE)						//非边界点
			{
				c1 = 2 / (GridLenth[j - 1] * (GridLenth[j - 1] + GridLenth[j]));
				c2 = 2 / (GridLenth[j] * (GridLenth[j - 1] + GridLenth[j]));			//计算c1,c2

				if (i == (Vertical - 1))			//计算c3,c4（需分类讨论）
				{
					c3 = 0;
					c4 = 4 / (pow(GridWidth[i - 1], 2));
				}
				else
				{
					c3 = (2 * Y[i] - GridWidth[i - 1]) / (Y[i] * GridWidth[i] * (GridWidth[i - 1] + GridWidth[i]));
					c4 = (2 * Y[i] + GridWidth[i]) / (Y[i] * GridWidth[i - 1] * (GridWidth[i - 1] + GridWidth[i]));
				}
				c0 = c1 + c2 + c3 + c4;

				if (i == (Vertical - 1))				//计算电位U
				{
					U = (1 - ω)*E[i][j] + ω * (c1*E[i][j - 1] + c2 * E[i][j + 1] + c4 * E[i - 1][j]) / c0;
				}
				else
				{
					U = (1 - ω)*E[i][j] + ω * (c1*E[i][j - 1] + c2 * E[i][j + 1] + c3 * E[i + 1][j] + c4 * E[i - 1][j]) / c0;
				}
				Res[i][j] = fabs(U - E[i][j]);
				E[i][j] = U;
			}
		}
	}
	ResSum = 0; ResMax = 0;
	for (int i = 0; i < Vertical; i++)				//残差和
	{
		for (int j = 0; j < Horizontal; j++)
		{
			ResSum = ResSum + Res[i][j];
			if (Res[i][j] > ResMax)
				ResMax = Res[i][j];
		}
	}
	ResAverage = ResSum / (Vertical*Horizontal);     //平均残差	
}

//计算电场分布
void ComputeElectricField()
{
	/////计算电场参数/////
	double λ;					 //残差平均值之比λ
	double acc;				  //加速因子
	double ωm;					//最佳迭代因子
	double ωmo;				  //上次迭代因子
	int Round = 1;                //迭代轮次
	int Times = 0;				 //迭代次数

								 //第一轮迭代
	ω = 1;
	printf("一轮迭代数据：\n第%d轮迭代\n", Round);
	fprintf(outfile, "一轮迭代数据：\n第%d轮迭代\n", Round);
	SOR();
	Times++;
	printf("迭代次数：%d\t迭代因子:%lf\t平均残差：%lf\t最大残差：%lf\n", Times, ω, ResAverage, ResMax);
	fprintf(outfile, "迭代次数：%d\t迭代因子:%lf\t平均残差：%lf\t最大残差：%lf\n", Times, ω, ResAverage, ResMax);
	Round++;

	//第二轮迭代
	ω = 1.375;
	printf("二轮迭代数据：\n第%d轮迭代\n", Round);
	fprintf(outfile, "二轮迭代数据：\n第%d轮迭代\n", Round);
	Times = 0;
	for (int t = 1; t <= 12; t++)			//12次迭代
	{
		SOR();
		Times++;
		if (t == 11)
			Va1 = ResAverage;
		else if (t == 12)
			Va2 = ResAverage;
	}
	printf("迭代次数：%d\t迭代因子:%lf\t平均残差：%lf\t最大残差：%lf\n", Times, ω, ResAverage, ResMax);
	fprintf(outfile, "迭代次数：%d\t迭代因子:%lf\t平均残差：%lf\t最大残差：%lf\n", Times, ω, ResAverage, ResMax);
	λ = Va2 / Va1;														//求残差平均值之比λ
	acc = 2 / (1 + sqrt(1 - pow((λ + ω - 1) / (sqrt(λ)*ω), 2)));    //加速因子
	ωm = 1.25*acc - 0.5;												//最佳迭代因子
	Round++;

	//第三轮迭代
	do
	{
		ω = ωm;
		printf("三轮迭代数据：\n第%d轮迭代\n", Round);
		fprintf(outfile, "三轮迭代数据：\n第%d轮迭代\n", Round);
		Times = 0;
		for (int t = 1; t <= 12; t++)
		{
			SOR();
			Times++;
			if (t == 11)
				Va1 = ResAverage;
			else if (t == 12)
				Va2 = ResAverage;

		}
		printf("迭代次数：%d\t迭代因子:%lf\t平均残差：%lf\t最大残差：%lf\n", Times, ω, ResAverage, ResMax);
		fprintf(outfile, "迭代次数：%d\t迭代因子:%lf\t平均残差：%lf\t最大残差：%lf\n", Times, ω, ResAverage, ResMax);
		λ = Va2 / Va1;
		ωmo = ωm;
		acc = 2 / (1 + sqrt(1 - pow((λ + ω - 1) / (sqrt(λ)*ω), 2)));
		ωm = 1.25*acc - 0.5;
	} while (fabs((ωm - ωmo) / (2 - ωmo)) >= e);               //判断ωm是否已经趋于稳定

															   //第Round轮迭代
	ω = ωm;
	int s;
	do
	{
		Round++;
		printf("第%d轮迭代", Round);
		fprintf(outfile, "\n第%d轮迭代\n", Round);
		SOR();
		Times++;
		printf("迭代次数：%d\t迭代因子:%lf\t平均残差：%lf\t最大残差：%lf\n", Times, ω, ResAverage, ResMax);
		fprintf(outfile, "迭代次数：%d\t迭代因子:%lf\t平均残差：%lf\t最大残差：%lf\n", Times, ω, ResAverage, ResMax);
		s = 0;
		for (int i = 0; i < Vertical; i++)
			for (int j = 0; j < Horizontal; j++)
			{
				if (Res[i][j] >= e)//迭代达到控制精度ε则停止
					s++;
			}
	} while (s > 0);

	//输出空间电场分布
	printf("\n格点上电位U:\n");
	fprintf(outfile, "\n格点上电位U:\n");
	for (int i = 0; i <= Horizontal; i++)
	{
		i += NST - 1;
		printf("%d\t", i);
		fprintf(outfile, "%d\t", i);
	}
	printf("\n");
	fprintf(outfile, "\n");
	for (int i = 0; i <= Vertical; i++)
	{
		printf("%d\t", i);
		fprintf(outfile, "%d\t", i);
		for (int j = 0; j <= Horizontal; j++)
		{
			fprintf(outfile, "%lf\t", E[i][j]);
			printf("%lf\t", E[i][j]);
			j += NST - 1;
		}
		i += NST - 1;
		printf("\n");
		fprintf(outfile, "\n");
	}
	printf("\n");
	fprintf(outfile, "\n");
}

////////////////////////7.绘制等位线////////////////////////////////
double *xx;				//x坐标
double *yy;             //y坐标
double *XX;				//放大后X坐标
double *YY;             //放大后Y坐标

void EqupotentialLine()
{
	int i, j, k, i1, j1, l, b, f1, f2, f3, f4, f5, *st;
	static double para, Xm, Ym, XYm;
	if (x == 1)
	{
		for (x = 1; para < (V[n - 1] - dV); x++)			//电位数
			para += dV;
		Ve = (double*)calloc(x - 1, sizeof(double));
		for (i = 0; Ve[i] < (V[n - 1] - dV); i++)
			Ve[i + 1] = Ve[i] + dV;							//电位值
	}

	//以下定义图形并计算等位线坐标
	initgraph(100 + int(10 * X[Horizontal - 1]), 100 + int(10 * Y[Vertical - 1]));		//初始图像
	setbkcolor(BLACK);																	//黑色背景
	cleardevice();
	rectangle(50, 50, 50 + int(10 * X[Horizontal - 1]), 50 + int(10 * Y[Vertical - 1]));//边框
	
	//计算电极位置
	Rank = (int*)calloc((n + 1) * sizeof(int), sizeof(int));
	Rank[0] = 0;
	for (i = 1; i < n + 1; i++)
	{
		if (i == n)
			Rank[i] = Rank[i - 1] + N[i - 1];
		else
			Rank[i] = Rank[i - 1] + N[i - 1] + 1;
	}
	for (i = 1; i < n; i++)
		bar(50 + int(10 * X[Rank[i] - 1]), 50, 50 + int(10 * (X[Rank[i]])), 50 + int(10 * dr2));   //画电极
	xx = (double*)calloc(Vertical*Horizontal, sizeof(double)); 
	yy = (double*)calloc(Vertical*Horizontal, sizeof(double));
	XX = (double*)calloc(Vertical*Horizontal, sizeof(double)); 
	YY = (double*)calloc(Vertical*Horizontal, sizeof(double));
	Xm = GridLenth[0];
	for (i = 0; i <= Horizontal; i++)
	{
		if (Xm < GridLenth[i])
			Xm = GridLenth[i];
	}
	Ym = GridWidth[0];
	if (Ym < GridWidth[Vertical - 2])
		Ym = GridWidth[Vertical - 2];
	XYm = sqrt(Xm*Xm + Ym * Ym);	 //计算对角线方向最大间距

	for (k = 0; k < x; k++)
	{
		int f1 = 0;									//顶行等电位点数量初始化
		l = 0;									    //等电位点数量初始化
		for (i1 = 0; i1<Vertical; i1++)
		{
			for (j1 = 0; j1 < Horizontal - 1; j1++)   //X方向上扫描赋电位
			{
				if (((E[i1][j1] < Ve[k])&(Ve[k] < E[i1][j1 + 1])) || ((E[i1][j1] > Ve[k])&(Ve[k] > E[i1][j1 + 1])) || (E[i1][j1] == Ve[k]))
				{
					yy[l] = Y[i1];
					if (yy[l] == Y[0])
						f1++;         
					if (E[i1][j1] == Ve[k])
						xx[l] = X[j1];
					else
						xx[l] = X[j1] + (Ve[k] - E[i1][j1])*GridLenth[j1] / (E[i1][j1 + 1] - E[i1][j1]);
					YY[l] = 50 + 10 * yy[l]; 
					XX[l] = 50 + 10 * xx[l];             //扩大10倍
					l++;
				}
			}

			//！！！！！！！！注意，在无鞍点情况下不需要纵向扫描，重复扫描会导致等位线无法连接！！！！！！！！！！！！！//
			for (j = 0; j < Horizontal; j++)        //Y方向上扫描赋电位
			{
				if (i1<Vertical - 1)
				{
				i = i1;
				if (((E[i][j] < Ve[k])&(Ve[k] < E[i + 1][j])) || ((E[i][j] > Ve[k])&(Ve[k] > E[i + 1][j])))
				{
					xx[l] = X[j];
					yy[l] = Y[i] + (Ve[k] - E[i][j])*GridWidth[i] / (E[i + 1][j] - E[i][j]);
					YY[l] = 50 + 10 * yy[l];
					XX[l] = 50 + 10 * xx[l];              //扩大10倍
					l++;
				}
			}

			}
		}
		//绘制等电位线
		int f = 0;   //连接线段起点位置变量初始化
		int c = 0;   //储存连接过的等电位点位置循环变量初始化
		int f6 = 0;  //用于记录电位等于电极电位时的电极位置初始化
		bool J = FALSE;
		st = (int*)calloc(l * sizeof(int), sizeof(int));   //储存连接过的等电位点位置
		do
		{
			fprintf(outfile, "\n电位为%lf坐标：\n", Ve[k]);
			setcolor(20 * k*RED + 10 * k*GREEN + BLUE);    //等位线颜色
			for (i = 0; i < l; i++)
			{
				J = TRUE;
				for (b = 0; b < c; b++)    //判断点状态
				{
					if (i == st[b])
					{
						J = FALSE;
						break;
					}
				}
				if ((J == TRUE)&(yy[i] == Y[0]))    
				{
					f = i;
					break;
				}
			}
			for (i = f; i < l; i++)
			{
				for (j = 0; j < l; j++)
				{
					J = TRUE;
					for (b = 0; b < c; b++)     //判断点状态
					{
						if (j == st[b])
						{
							J = FALSE;
							break;
						}
					}
					if ((j != i)&(J == TRUE)&(sqrt((XX[j] - XX[i])*(XX[j] - XX[i]) + (YY[j] - YY[i])*(YY[j] - YY[i])) < 12 * XYm)&(fabs(XX[j] - XX[i]) < 12 * Xm)&(fabs(YY[j] - YY[i]) < 12 * Ym))
					{    
						f3 = 0;
						for (f2 = 0; f2 < n; f2++)    //标记电极
						{
							if (Ve[k] == V[f2])
							{
								f3 = 1;
								f6 = f2;
								break;
							}
						}
						if ((f3 == 1)&((X[Rank[f6 + 1]] == xx[j]) || (X[Rank[f6 + 1] - 1] == xx[j]))&(Y[0] == yy[j]))	//Z字连接
						{    
							for (f2 = 0; f2 < M2; f2++)
							{
								fprintf(outfile, "(%lf,%lf)\n", yy[i], xx[i]);  //第一点坐标
								if (yy[i] == Y[0])
									f1--;   
								line(int(XX[i]), int(YY[i]), int(XX[j]), int(YY[j]));
								st[c] = i;      
								c++;
								fprintf(outfile, "(%lf,%lf)\n", yy[j], xx[j]);   //第二点坐标
								if (yy[j] == Y[0])
									f1--;		
								for (f4 = 0; f4 < l; f4++)
								{
									if ((xx[f4] == xx[i])&(yy[f4] == yy[i] + dr2 / M2))
									{
										i = f4;					//重复寻找
										break;
									}
								}

								line(int(XX[j]), int(YY[j]), int(XX[i]), int(YY[i]));
								st[c] = j;
								c ++;
								for (f5 = 0; f5 < l; f5++)
								{
									if ((xx[f5] == xx[j])&(yy[f5] == yy[j] + dr2 / M2))
									{
										j = f5;     //重复寻找
										break;
									}
								}
							}
							fprintf(outfile, "(%lf,%lf)\n", yy[i], xx[i]);  //电极底端电位点坐标
							if (yy[i] == Y[0])
								f1--;		
							line(int(XX[i]), int(YY[i]), int(XX[j]), int(YY[j]));   //连接底端
							st[c] = i;
							c++;
							break;
						}
						//下面循环方式相同，以下略
						else if ((f3 == 1)&((X[Rank[f6 + 1]] == xx[j]) || (X[Rank[f6 + 1] - 1] == xx[j]))&(Y[0] != yy[j]))	//自下而上Z字连接
						{   
							fprintf(outfile, "(%lf,%lf)\n", yy[i], xx[i]);
							if (yy[i] == Y[0])
								f1--;
							line(int(XX[i]), int(YY[i]), int(XX[j]), int(YY[j]));
							st[c] = i;
							c++;
							fprintf(outfile, "(%lf,%lf)\n", yy[j], xx[j]);
							if (yy[j] == Y[0])
								f1--;
							for (f4 = 0; f4 < l; f4++)
							{
								if ((xx[f4] == xx[Rank[f6 + 1]])&(yy[f4] == yy[M2]))
								{
									i = f4;
									j = f4 - 1;
									break;
								}
							}
							for (f2 = 0; f2 < M2; f2++)
							{
								fprintf(outfile, "(%lf,%lf)\n", yy[i], xx[i]);
								if (yy[i] == Y[0])
									f1--;
								line(int(XX[i]), int(YY[i]), int(XX[j]), int(YY[j]));
								st[c] = i;
								c++;
								fprintf(outfile, "(%lf,%lf)\n", yy[j], xx[j]);
								if (yy[j] == Y[0])
									f1--;
								for (f4 = 0; f4 < l; f4++)
								{
									if ((xx[f4] == xx[i])&(yy[f4] == yy[i] - dr2 / M2))
									{
										i = f4;
										break;
									}
								}
								line(int(XX[j]), int(YY[j]), int(XX[i]), int(YY[i]));
								st[c] = j;
								c++;
								for (f5 = 0; f5 < l; f5++)
								{
									if ((xx[f5] == xx[j])&(yy[f5] == yy[j] - dr2 / M2))
									{
										j = f5;
										break;
									}
								}
							}
							fprintf(outfile, "(%lf,%lf)\n", yy[i], xx[i]);
							if (yy[i] == Y[0])
								f1--;
							line(int(XX[i]), int(YY[i]), int(XX[j]), int(YY[j]));
							st[c] = i;
							c++;
							break;
						}
						else
						{  
							fprintf(outfile, "(%lf,%lf)\n", yy[i], xx[i]);
							if (yy[i] == Y[0])
								f1--;
							line(int(XX[i]), int(YY[i]), int(XX[j]), int(YY[j]));
							st[c] = i;
							c++;
							break;
						}
					}
				}
				if ((yy[j] == Y[0]) || (yy[j] == Y[Vertical - 1]))
				{
					fprintf(outfile, "(%lf,%lf)\n", yy[j], xx[j]);
					if (yy[j] == Y[0])
						f1--;
					st[c] = j;
					c++;
					break;
				}
				else
					i = j - 1;
			}
		} while (f1 > 0);
	}
}
