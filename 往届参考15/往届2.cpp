#include<StdAfx.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<graphics.h>

void rwfile();
void initialize();
void SOR();
void calculate();
void equpotential_line();

//主函数
void main()
{
	rwfile();//读写子函数
	initialize();//初始化子函数
	calculate();//计算子函数（其中用到SOR子函数）
	equpotential_line();//等位线子函数
}

//参数定义
FILE *infile;				//原始数据文件
FILE *outfile;				//输出结果文件
double δ;				//电极厚度	
int n;					//电极总数	
double *dz;				//电极间距
int *N;					//相邻电极间要划分的步长数
double *dd;				//电极长度
int *L;					//电极要划分的步长数
double *V;				//电极电位
double r;				//电极内孔径半径
int M;					//r范围内等步长划分网格数
double e;				//迭代控制精度ε
int NST;				//输出打印空间电位时网格点间隔数
double dv;				//要求扫描搜索等电位线的电位间隔值
int row=0,col=0;		//行数row，列数col
double **E;				//格点电位值
double *GridLenth;				//每个小格子长度
double *X;				//格点z坐标
double *GridWidth;				//每个小格子宽度
double *Y;				//格点r坐标
double c1,c2,c3,c4,c0;	//五点差分公式中系数
double ω;				//SOR中超张弛迭代因子
double U;				//中心点电位
double **Res;			//格点电位残差
double ResSum;            //格点电位残差和
double ResMax;            //格点电位最大残差
double ResAverage;            //格点电位平均残差
double Va1;           //格点电位平均残差
double Va2;           //格点电位平均残差
int T;                  //迭代轮次
int Time;               //迭代次数
int TimeS;              //迭代总次数
double λ;              //残差平均值之比λ
double ωλ;            //加速因子ωλ
double μλ;            //ωλ中的中间变量
double ωm;             //修正的最佳迭代因子
double ωm1;            //修正ω因子
int x=0;				//要求扫描搜索等电位线的电位值个数
double *Ve;				//要求扫描搜索等电位线的电位值
double *ze;				//等电位点的z坐标
double *re;             //等电位点的r坐标
double *zeb;			//图中放大的等电位点坐标z坐标
double *reb;            //图中放大的等电位点坐标r坐标
/*
double **VeZ;
double **VeR;
double **VeV;
*/
int *Rank;				//电极右边界格点位置（列数）

//读写函数
void rwfile()
{
	int i,j;
	infile=fopen("C:\\CADtest\\1120151052\\1120151052_李中石.dat","r");     //读入数据文件
	if (infile==NULL)
	{
		printf("数据读取失败！\n");//读入失败就报错
		exit(0);
	}
	fscanf(infile,"电极厚度δ=%lfmm;\n",&δ);//从数据文件里读取数据
	fscanf(infile,"电极总数n=%ld;\n",&n);
	fscanf(infile,"电极间距dz=");
	dz=(double*)calloc(n*sizeof(double),sizeof(double));
	for (i=0;i<n;i++)
		fscanf(infile,"%lfmm;",&dz[i]);
	fscanf(infile,"\n相邻电极间要划分的步长数N=");
	N=(int*)calloc(n*sizeof(int),sizeof(int));
	for (i=0;i<n;i++)
		fscanf(infile,"%ld;",&N[i]);
	fscanf(infile,"\n电极长度dd=");
	dd=(double*)calloc((n-1)*sizeof(double),sizeof(double));
	for (i=0;i<n-1;i++)
		fscanf(infile,"%lfmm;",&dd[i]);
	fscanf(infile,"\n电极要划分的步长数L=");
	L=(int*)calloc((n-1)*sizeof(int),sizeof(int));
	for (i=0;i<n-1;i++)
		fscanf(infile,"%ld;",&L[i]);
	fscanf(infile,"\n电极电位V=");
	V=(double*)calloc(n*sizeof(double),sizeof(double));
	for (i=0;i<n;i++)
		fscanf(infile,"%lfV;",&V[i]);
	fscanf(infile,"\n电极内孔径半径r=%lfmm;\n",&r);
	fscanf(infile,"r范围内等步长划分网格数M=%ld;\n",&M);
	fscanf(infile,"迭代控制精度ε=%lfV;\n",&e);
	fscanf(infile,"输出打印空间电位时网格点间隔数NST=%ld;\n",&NST);
	fscanf(infile,"等电位线的电位间隔或电位值=");
	Ve=(double*)calloc(100*sizeof(double),sizeof(double));//给要求扫描搜索等电位线的电位值分配空间
	while (fscanf(infile,"%lfV;",&Ve[x])!=EOF)
		x++;//要求扫描的电位值个数
	if (x==1)//若该项只有一个数则输入的是dv
		dv=Ve[0];
	fclose(infile);//关掉读入的目标文件

	outfile=fopen("C:\\CADtest\\1120151052\\1120151052_李中石.res","w");//生成输出文件
	if (outfile==NULL)//报错
	{
		printf("Cannot open file.\n");
		exit(0);
	}
	fprintf(outfile,"\t\t\t\t ==========\n");
	fprintf(outfile,"\t\t\t\t【原始数据】\n");//将原始数据重新写出来
	fprintf(outfile,"\t\t\t\t ==========\n");
	fprintf(outfile,"电极厚度δ=%lfmm;\n",δ);
	fprintf(outfile,"电极总数n=%ld;\n",n);
	fprintf(outfile,"相邻电极间距离∆zi=");
	for (i=0;i<n;i++)
		fprintf(outfile,"%lfmm;",dz[i]);
	fprintf(outfile,"\n相邻电极间要划分的步长数Ni=");
	for (i=0;i<n;i++)
		fprintf(outfile,"%ld;",N[i]);
	fprintf(outfile,"\n电极长度∆di=");
	for (i=0;i<n-1;i++)
		fprintf(outfile,"%lfmm;",dd[i]);
	fprintf(outfile,"\n电极要划分的步长数Li=");
	for (i=0;i<n-1;i++)
		fprintf(outfile,"%ld;",L[i]);
	fprintf(outfile,"\n电极电位Vi=");
	for (i=0;i<n;i++)
		fprintf(outfile,"%lfV;",V[i]);
	fprintf(outfile,"\n电极内孔径半径r=%lfmm;\n",r);
	fprintf(outfile,"r范围内等步长划分网格数M=%ld;\n",M);
	fprintf(outfile,"迭代控制精度ε=%lfV;\n",e);
	fprintf(outfile,"输出打印空间电位时网格点间隔数NST=%ld;\n",NST);
	fprintf(outfile,"等电位线的电位间隔∆V或电位值EVi=");
	for (i=0;i<x;i++)
		fprintf(outfile,"%lfV;",Ve[i]);
	fprintf(outfile,"\n\n\t\t\t\t ============\n");
	fprintf(outfile,"\t\t\t\t【网格点划分】\n");
	fprintf(outfile,"\t\t\t\t ============\n");
}

//初始化函数
void initialize()
{
	int i,j;//定义变量
	int k=0;
	row=M+2;//r方向的行数=网格数+电极厚度（1）+1
	for (i=0;i<n;i++)
		col=col+N[i];
	for (i=0;i<n-1;i++)
		col=col+L[i];
	col=col+1;//z方向的列数=网格数+1
	printf("步长划分完毕。共 %d 行 ， %d 列；\n电位运算结果已输出至原目录下Output1.dat。\n",row,col);

	Rank=(int*)calloc((n+1)*sizeof(int),sizeof(int));//为电极右侧点位置（对应的列数）分配内存
	Rank[0]=0;//令dj[0]=0
	for(i=1;i<n+1;i++)//为dj的每个元素赋值
	{
		if (i==n)
			Rank[i]=Rank[i-1]+N[i-1];//荧光屏处（相当于第n个电极右边界）格点位置
		else
			Rank[i]=Rank[i-1]+N[i-1]+L[i-1];//一般电极右边界格点位置
	}

	E=(double**)calloc(row*sizeof(double*),sizeof(double));  //为格点的电场分配内存（二维）
	for (i=0;i<row;i++)
		E[i]=(double*)calloc(col*sizeof(double),sizeof(double));

	for (i=0;i<row;i++)//赋E初值全部为零
	{
		for (j=0;j<col;j++)
		{
			E[i][j]=0;
		}
	}

	for (i=0;i<n;i++)//电极（等势体）上第0行第1行应该相等，都赋值成电极电压
	{
		for (j=(k+N[i]);j<=(k+N[i]+L[i]);j++)
		{
			E[0][j]=V[i];
			E[1][j]=V[i];
		}
		k=k+N[i]+L[i];
	}

	for (i=0;i<row;i++)//荧光屏（整列）全部赋值V[n-1]
	{
		E[i][col-1]=V[n-1];
	}
	k=1;

	for (i=0;i<n;i++)//为电极之间电压赋值（插值）
	{
		for (j=k;j<(k+N[i]);j++)
		{
			if (i==0)
			{
				E[0][j]=E[0][j-1]+V[i]/N[i];
			}
			else
			{
				E[0][j]=E[0][j-1]+(V[i]-V[i-1])/N[i];
			}
		}
		k=k+N[i]+L[i];
	}
	//至此，电场电位初始化over

	GridLenth=(double*)calloc((col-1)*sizeof(double),sizeof(double));//为每个格子的（z向）长度分配内存
	k=0;
	for (i=0;i<n;i++)//计算每个格子的（z向）长度
	{
		for (j=k;j<(k+N[i]);j++)//计算电极之间的格子的（z向）长度
		{
			GridLenth[j]=dz[i]/N[i];
		}
		k=k+N[i];
		if (i==n-1)//若怼到屏上即跳出
			break;
		for (j=k;j<(k+L[i]);j++)//计算电极内部的格子的（z向）长度
		{
			GridLenth[j]=dd[i]/L[i];
		}
		k=k+L[i];
	}

	fprintf(outfile,"[z坐标方向]\n");//输出每个格子的（z向）长度
	for (i=0;i<col;i++)
		fprintf(outfile,"%ld\t\t",i);
	fprintf(outfile,"\n");
	for (i=0;i<(col-1);i++)
		fprintf(outfile,"%lf\t",GridLenth[i]);

	X=(double*)calloc(col*sizeof(double),sizeof(double));//为每个位置的（z向）横坐标分配内存
	X[0]=0;
	for (i=1;i<col;i++)
	{
		X[i]=X[i-1]+GridLenth[i-1];
	}

	GridWidth=(double*)calloc((row-1)*sizeof(double),sizeof(double));//为每个格子的（r向）宽度分配内存
	for (i=(row-2);i>0;i--)//计算每个格子的（r向）宽度
	{
		GridWidth[i]=r/M;
	}
	GridWidth[0]=δ;
	fprintf(outfile,"\n[r坐标方向]\n");//输出每个格子的（r向）宽度
	for (i=0;i<row;i++)
		fprintf(outfile,"%ld\t\t",i);
	fprintf(outfile,"\n");
	for (i=0;i<(row-1);i++)
		fprintf(outfile,"%lf\t",GridWidth[i]);
	fprintf(outfile,"\n\n");
	Y=(double*)calloc(row*sizeof(double),sizeof(double));//为每个位置的（r向）纵坐标分配内存
	Y[row-1]=0;
	for (i=(row-2);i>=0;i--)
	{
		Y[i]=Y[i+1]+GridWidth[i];
	}
	//至此坐标计算over，像管初始化over
}

//SOR迭代算法函数
void SOR()
{
	int i,j,k;
	double h1,h2,h3,h4,r0;

	Res=(double**)calloc(row*sizeof(double*),sizeof(double));  //为每个格子的电位残差分配内存（二维）
	for (i=0;i<row;i++)
		Res[i]=(double*)calloc(col*sizeof(double),sizeof(double));

	for (i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			int flag=1;//边界点判断参量：边界1，非边界0
			int K=0;//计数器
			if ((i!=0)&&(j!=0)&&(j!=col-1))//上、左、右边界点判断
			{
				flag=0;
					for (k=0;k<n;k++)//电极位置边界点也不参与计算
					{
						if (i==1)
						{
						if (j>=(K+N[k])&&j<=(K+N[k]+L[k]))
						{
							flag=1;
							break;
						}
						}
						K=K+N[k]+L[k];
					}
			}
			if (flag==0)//其余的非边界点参与计算
			{
				c1=2/(GridLenth[j-1]*(GridLenth[j-1]+GridLenth[j]));//SOR中c系数计算 c1=2/h1(h1+h2)
		        	c2=2/(GridLenth[j]*(GridLenth[j-1]+GridLenth[j]));//c2=2/h2(h1+h2)

		        	if (i==(row-1))//对于c3,c4需要分类讨论
		     		{
			     		c3=0;//最下面一行c3=0
			        	c4=4/(GridWidth[i-1]*GridWidth[i-1]);//最下面一行c4=4/(h4*h4)
		        	}
		        	else
		        	{
			        	c3=(2*Y[i]-GridWidth[i-1])/(Y[i]*GridWidth[i]*(GridWidth[i-1]+GridWidth[i]));//其余情况的c3
			        	c4=(2*Y[i]+GridWidth[i])/(Y[i]*GridWidth[i-1]*(GridWidth[i-1]+GridWidth[i]));//其余情况的c4
		        	}
		        	c0=c1+c2+c3+c4;

				if (i==(row-1))//对于电位φ0需要分类讨论
				{
					U=(1-ω)*E[i][j]+ω*(c1*E[i][j-1]+c2*E[i][j+1]+c4*E[i-1][j])/c0;
				}
				else 
				{
					U=(1-ω)*E[i][j]+ω*(c1*E[i][j-1]+c2*E[i][j+1]+c3*E[i+1][j]+c4*E[i-1][j])/c0;
				}
				Res[i][j]=fabs(U-E[i][j]);//计算每个点的电位残差（求绝对值）
				E[i][j]=U;
			}	
		}
	}
	ResSum=0;//初始化残差和
	ResMax=0;//初始化最大残差
	for (i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			ResSum=ResSum+Res[i][j];
			if (Res[i][j]>ResMax)
				ResMax=Res[i][j];
		}
	}
	ResAverage=ResSum/(row*col);//计算平均残差	
}

//计算电场分布函数
void calculate()
{
	int i,j,k,sign;//计数器及迭代标志
	fprintf(outfile,"\t\t\t\t ==========\n\t\t\t\t【迭代过程】\n\t\t\t\t ==========\n");
	T=1;//迭代轮次
	Time=1;
	TimeS=0;//总迭代次数
	ω=1;//超张弛迭代因子
	SOR();//第一轮迭代
	TimeS=TimeS+1;//总迭代次数+1
	fprintf(outfile,"迭代轮次：%-3ld\t\t迭代次数：%-3ld\t\t迭代因子值:%lf\t\t",T,TimeS,ω);
	fprintf(outfile,"平均残差：%lf\t\t最大残差：%lf\n",ResAverage,ResMax);

	ω=1.375;//此即“卡瑞建议值”
	T=T+1;//迭代轮次+1
	for (Time=1;Time<13;Time++)
	{
		SOR();//第二轮迭代
		TimeS=TimeS+1;
		if (Time==11)
			Va1=ResAverage;
		else if (Time==12)
			Va2=ResAverage;
	}
	fprintf(outfile,"迭代轮次：%-3ld\t\t迭代次数：%-3ld\t\t迭代因子值:%lf\t\t",T,TimeS,ω);
	fprintf(outfile,"平均残差：%lf\t\t最大残差：%lf\n",ResAverage,ResMax);
	λ=Va2/Va1;//求残差平均值之比λ
	μλ=(λ+ω-1)/(sqrt(λ)*ω);//求中间变量μλ
	ωλ=2/(1+sqrt(1-μλ*μλ));//求加速因子ωλ
	ωm=1.25*ωλ-0.5;//求得最佳迭代因子

	do
	{
		ω=ωm;
	    T=T+1;
	    for (Time=1;Time<13;Time++)
	    {
	    	SOR();//第三轮迭代
	    	TimeS=TimeS+1;
		    if (Time==11)
			    Va1=ResAverage;
	    	else if (Time==12)
			    Va2=ResAverage;
	    }
	    fprintf(outfile,"迭代轮次：%-3ld\t\t迭代次数：%-3ld\t\t迭代因子值:%lf\t\t",T,TimeS,ω);
	    fprintf(outfile,"平均残差：%lf\t\t最大残差：%lf\n",ResAverage,ResMax);
		ωm1=ωm;
	    λ=Va2/Va1;
	    μλ=(λ+ω-1)/(sqrt(λ)*ω);
	    ωλ=2/(1+sqrt(1-μλ*μλ));
	    ωm=1.25*ωλ-0.5;
	}
	while(fabs((ωm-ωm1)/(2-ωm1))>=0.05);//连续两轮求得ωm1和ωm满足该式，则固定ωm1
	ω=ωm;

	do
	{
		T=T+1;
		SOR();//最佳迭代因子确定后继续迭代
		TimeS=TimeS+1;
		fprintf(outfile,"迭代轮次：%-3ld\t\t迭代次数：%-3ld\t\t迭代因子值:%lf\t\t",T,TimeS,ω);
	    fprintf(outfile,"平均残差：%lf\t\t最大残差：%lf\n",ResAverage,ResMax);
		sign=0;
		for (i=0;i<row;i++)
			for(j=0;j<col;j++)
			{
				if (Res[i][j]>=e)//迭代达到控制精度ε则停止
					sign=sign+1;
			}
	}
	while(sign>0);
	//至此，像管中每格点的电位计算完毕
	
	//输出计算得到的空间电场分布
	fprintf(outfile,"\n\t\t\t\t ==================\n");
	fprintf(outfile,"\t\t\t\t【网格点上的电位值】\n");
	fprintf(outfile,"\t\t\t\t ==================\n\t\t");
	for (i=0;i<col;i++)
	{
		fprintf(outfile,"%ld\t\t",i);
		i=i+NST-1;
	}
	fprintf(outfile,"\n");
	for (i=row-NST;i>=0;i--)
	{
		fprintf(outfile,"%ld\t",i);
		for (j=0;j<col;j++)
		{
			fprintf(outfile,"%lf\t",E[i][j]);
			j=j+NST-1;
		}
		fprintf(outfile,"\n");
		i=i-NST+1;
	}
}

//绘制等位线图
void equpotential_line()
{
	int i,i1,j,j1,k,l,c,b,pan,f,f1,f2,f3,f4,f5,f6,*hi;
	double a=0,zlm,rlm,zrlm;
	if (x==1)//如果输入的是间隔值（只有一个数）则计算等电位值
	{
		for (x=1;a<(V[n-1]-dv);x++)//计算电位数量
	    	a=a+dv;
	    Ve=(double*)calloc((x-1)*sizeof(double),sizeof(double));
    	for (i=0;Ve[i]<(V[n-1]-dv);i++)
	        Ve[i+1]=Ve[i]+dv;//计算各等位线的电位值
	}
	
	//以下定义图形并计算等位线坐标
	initgraph(100+int(10*X[col-1]),100+int(10*Y[row-1]));//画面距窗口边缘四周均50，坐标值均扩大10倍
	setbkcolor(BLACK);//设定黑色背景
	cleardevice();
    rectangle(50,50,50+int(10*X[col-1]),50+int(10*Y[row-1]));//距窗口边缘四周均50画矩形
	for (i=1;i<n;i++)
		bar(50+int(10*X[Rank[i]-1]),50,50+int(10*(X[Rank[i]])),50+int(10*δ));   //画电极
	ze=(double*)calloc(row*col*sizeof(double),sizeof(double));   //等电位点z坐标值分配内存
	re=(double*)calloc(row*col*sizeof(double),sizeof(double));   //等电位点r坐标值分配内存
	zeb=(double*)calloc(row*col*sizeof(double),sizeof(double));  //扩大10倍的等电位点z坐标值分配内存
	reb=(double*)calloc(row*col*sizeof(double),sizeof(double));  //扩大10倍的等电位点r坐标值分配内存
	zlm=GridLenth[0];
	for (i=0;i<col-1;i++)
	{
		if (zlm<GridLenth[i])
			zlm=GridLenth[i];        //计算z方向最大格点间距
	}
	rlm=GridWidth[0];
	if (rlm<GridWidth[row-2])
		rlm=GridWidth[row-2];        //计算r方向最大格点间距
	zrlm=sqrt(zlm*zlm+rlm*rlm);	 //计算对角线方向最大间距
	for (k=1;k<x;k++)
	{
		f1=0;    //记录顶行等电位点数量初始化
		l=0;     //等电位点数量初始化
		for (i1=row-1;i1>=0;i1--)
		{
			for (j1=0;j1<col-1;j1++)   //z轴方向线性插值计算等电位点坐标
			{
				if (((E[i1][j1]<Ve[k])&(Ve[k]<E[i1][j1+1]))||((E[i1][j1]>Ve[k])&(Ve[k]>E[i1][j1+1]))||(E[i1][j1]==Ve[k]))
				{	
					re[l]=Y[i1];
					if (re[l]==Y[row-1])
						f1++;          //顶行标记值记录顶行等电位点数量用于控制等位线绘制
					if (E[i1][j1]==Ve[k])
						ze[l]=X[j1];
					else
						ze[l]=X[j1]+(Ve[k]-E[i1][j1])*GridLenth[j1]/(E[i1][j1+1]-E[i1][j1]);
					reb[l]=50+10*(Y[row-1]-re[l]);  //扩大10倍后的图形坐标 
					zeb[l]=50+10*ze[l];             //扩大10倍后的图形坐标 
					l++;
				}
			}
			for (j=0;j<col;j++)        //r轴方向线性插值计算等电位点坐标
			{
				if (i1>0)
				{
					i=i1;
					if (((E[i][j]<Ve[k])&(Ve[k]<E[i-1][j]))||((E[i][j]>Ve[k])&(Ve[k]>E[i-1][j])))
					{
						ze[l]=X[j];
						re[l]=Y[i]-(Ve[k]-E[i][j])*GridWidth[i-1]/(E[i-1][j]-E[i][j]);
						reb[l]=50+10*(Y[row-1]-re[l]);   //扩大10倍后的图形坐标
						zeb[l]=50+10*ze[l];              //扩大10倍后的图形坐标
						l++;
					}
				}
			}
		}
		
		//////////////////绘制等位线//////////////////
		f=0;   //连接线段起点位置变量初始化
		c=0;   //储存连接过的等电位点位置循环变量初始化
		f6=0;  //用于记录电位等于电极电位时的电极位置初始化
		hi=(int*)calloc(l*sizeof(int),sizeof(int));   //储存连接过的等电位点位置
		do
		{		
			fprintf(outfile,"\n电位值为%lf的各点坐标值：\n",Ve[k]);
		    setcolor(20*k*RED+10*k*GREEN+BLUE);    //确定等位线颜色
			for (i=0;i<l;i++)
			{
				pan=1;
				for (b=0;b<c;b++)    //判断是否为已连接过的点
				{
					if (i==hi[b])
					{
						pan=0;
						break;
					}
				}
				if ((pan==1)&(re[i]==Y[row-1]))    //用于有鞍点时找出下一条等位线的起点位置
				{
					f=i;
					break;
				}
			}
			for (i=f;i<l;i++)
			{
				for (j=0;j<l;j++)
				{
					pan=1;
					for (b=0;b<c;b++)    //判断是否为已连接过的点
					{
						if (j==hi[b])
						{
							pan=0;
							break;
						}
					}	
					if ((j!=i)&(pan==1)&(sqrt((zeb[j]-zeb[i])*(zeb[j]-zeb[i])+(reb[j]-reb[i])*(reb[j]-reb[i]))<12*zrlm)&(fabs(zeb[j]-zeb[i])<12*zlm)&(fabs(reb[j]-reb[i])<12*rlm))
					{    //如果第二点不是第一点且不是已连接过的点且两点之间距离满足要求
						f3=0;
						for (f2=0;f2<n;f2++)    //如果等于电极电位则记录电极位置并做标记
						{
							if (Ve[k]==V[f2])
							{
								f3=1;    //做标记
								f6=f2;   //记录电极位置
								break;
							}				
						}				
						if ((f3==1)&((X[Rank[f6+1]]==ze[j])||(X[Rank[f6+1]-1]==ze[j]))&(Y[row-1]==re[j]))
						{    //如果是电极电位且第二点为电极顶部所在位置则从上往下“Z”字型连接电极电位点
							for (f2=0;f2<1;f2++)
							{
								fprintf(outfile,"(%lf,%lf)\n",re[i],ze[i]);   //输出第一点坐标
								if (re[i]==Y[row-1])
									f1--;     //如果第一点在顶行则标记值减1
								line(int(zeb[i]),int(reb[i]),int(zeb[j]),int(reb[j])); 
								hi[c]=i;      //记录第一点位置
								c++;
								fprintf(outfile,"(%lf,%lf)\n",re[j],ze[j]);   //输出第二点坐标
								if (re[j]==Y[row-1])						
									f1--;	  //如果第二点在顶行则标记值减1			
							    for (f4=0;f4<l;f4++)
								{
									if ((ze[f4]==ze[i])&(re[f4]==re[i]-δ/1))
									{
										i=f4;      //寻找下一个第一点
										break;
									}
								}
								line(int(zeb[j]),int(reb[j]),int(zeb[i]),int(reb[i]));
								hi[c]=j;
								c++;
								for (f5=0;f5<l;f5++)
								{
									if ((ze[f5]==ze[j])&(re[f5]==re[j]-δ/1))
									{
										j=f5;      //寻找下一个第二点
										break;
									}
								}					    
							}
							fprintf(outfile,"(%lf,%lf)\n",re[i],ze[i]);  //输出电极底端左侧电位点坐标
							if (re[i]==Y[row-1])							
								f1--;			//如果此点在顶行则标记值减1		
							line(int(zeb[i]),int(reb[i]),int(zeb[j]),int(reb[j]));   //连接电极底端水平线
							hi[c]=i;            //记录电极底端左侧电位点位置
							c++;
							break;					
						}
						else if ((f3==1)&((X[Rank[f6+1]]==ze[j])||(X[Rank[f6+1]-1]==ze[j]))&(Y[row-1]!=re[j]))
						{   //如果是电极电位且第二点不是电极顶部所在位置则从下往上倒“Z”字型连接电极电位点
							fprintf(outfile,"(%lf,%lf)\n",re[i],ze[i]);
							if (re[i]==Y[row-1])
								f1--;
							line(int(zeb[i]),int(reb[i]),int(zeb[j]),int(reb[j]));
							hi[c]=i;
							c++;
							fprintf(outfile,"(%lf,%lf)\n",re[j],ze[j]);
							if (re[j]==Y[row-1])
								f1--;																		
							for (f4=0;f4<l;f4++)
							{
								if ((ze[f4]==ze[Rank[f6+1]])&(re[f4]==re[M]))
								{
									i=f4;
									j=f4-1;
									break;
								}
							}											    
					    	for (f2=0;f2<1;f2++)
							{
								fprintf(outfile,"(%lf,%lf)\n",re[i],ze[i]);
								if (re[i]==Y[row-1])					
						    		f1--;
								line(int(zeb[i]),int(reb[i]),int(zeb[j]),int(reb[j]));
								hi[c]=i;
								c++;
								fprintf(outfile,"(%lf,%lf)\n",re[j],ze[j]);
								if (re[j]==Y[row-1])			
					 		    	f1--;					 
					    		for (f4=0;f4<l;f4++)
								{
									if ((ze[f4]==ze[i])&(re[f4]==re[i]+δ))
									{
										i=f4;
										break;
									}
								}
								line(int(zeb[j]),int(reb[j]),int(zeb[i]),int(reb[i]));
								hi[c]=j;
								c++;
								for (f5=0;f5<l;f5++)
								{
									if ((ze[f5]==ze[j])&(re[f5]==re[j]+δ))
									{
										j=f5;
										break;
									}
								}				    
							}
							fprintf(outfile,"(%lf,%lf)\n",re[i],ze[i]);
							if (re[i]==Y[row-1])					
								f1--;
							line(int(zeb[i]),int(reb[i]),int(zeb[j]),int(reb[j]));
							hi[c]=i;
							c++;
							break;
						}
						else
						{   //如果第二点不是电极所在位置则输出第一点坐标并连接第一二点
							fprintf(outfile,"(%lf,%lf)\n",re[i],ze[i]);
							if (re[i]==Y[row-1])							
								f1--;						
							line(int(zeb[i]),int(reb[i]),int(zeb[j]),int(reb[j]));
							hi[c]=i;
							c++;
							break;
						}
				 }
				}
				if ((GridWidth[j]==Y[0])||(GridWidth[j]==Y[row-1]))
				{   //如果第二点位于顶行或底行则输出第二点坐标并跳出循环判断是否有下一条等位线
					fprintf(outfile,"(%lf,%lf)\n",re[j],ze[j]);
					if (re[j]==Y[row-1])
						f1--;
					hi[c]=j;
					c++;				
		     		break;
				}
				else   //如果第二点不位于顶行或底行则将第二点作为下一次循环的第一点
					i=j-1;
	 }
	}
	while (f1>0);      //如果顶行标记值不为0则继续搜索连接下一条等位线
	}
	fclose(outfile);
	system("pause");
	closegraph();

}