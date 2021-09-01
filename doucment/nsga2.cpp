#include<stdio.h>
#include<stdlib.h>
#include<Windows.h>
#include<math.h>
#include<time.h>
#include<iostream>
#define Dimension 30//����ά���������ＴZDT1����xi��i�����ֵ
#define popsize 100//��Ⱥ��С
#define generation 500 //���ܴ���
#define URAND (rand()/(RAND_MAX+1.0))//���������

int temp1[popsize];//��ʱ����
int mark[popsize];//�������
//���������������ڲ����µ��Ӵ�
using namespace std;
//�����������
class individual
{
public:
    double value[Dimension];//xi��ֵ
    int sp[2*popsize];
    //��֧����弯��SP�������ǿ��н�ռ������б�����p֧��ĸ�����ɵļ��ϡ�
    int np;
    //֧�����np���������ڿ��н�ռ��п���֧�����p�����Ը����������
    int is_dominated;//����sp�ĸ���
    void init();//��ʼ������
    int rank;//���ȼ���Pareto����Ϊ��ǰ��߼�
    double crowding_distance;//ӵ������
    double fvalue[2];//ZDT1����Ŀ�꺯����ֵ
    void f_count();//����fvalue��ֵ
};
//Ⱥ���������
class population
{
public:
    population();//���ʼ��
    individual P[popsize];
    individual Q[popsize];
    individual R[2*popsize];
    void set_p_q();
    //�������һ����ʼ����P���ڴ˻����ϲ��ö�Ԫ������ѡ��
    //����ͱ�����������Ӵ�Q��P��QȺ���ģ��Ϊpopsize
    //��Pt��Qt���뵽Rt�У���ʼʱt=0������Rt���п��ٷ�֧�������
    //���������в�ͬ�ȼ��ķ�֧��⼯F1��F2........
    int Rnum;
    int Pnum;
    int Qnum;
    //P,Q,R��Ԫ�صĸ���
    void make_new_pop();//�����µ��Ӵ�
    void fast_nondominated_sort();//���ٷ�֧������
    void calu_crowding_distance(int i);//ӵ���������
    void f_sort(int i);//��ӵ�����뽵������
    void maincal();//��Ҫ����
    int choice(int a,int b);
    //�����������ڲ�ͬ�ȼ��ķ�֧��⼯�����ȿ��ǵȼ���Ž�С��
    //��������������ͬһ�ȼ��ķ�֧��⼯�����ȿ���ӵ������ϴ��
    int len[2*popsize];//�������콻����Ⱥ��Fi�ĳ��ȵļ���
    int len_f;//����Ⱥ��rankֵ
};
//ȫ�ֱ��������ֺ�������
individual F[2*popsize][2*popsize];

double rand_real(double low,double high)
//�������ʵ��
{
    double h;
    h=(high-low)*URAND+low+0.001;
    if(h>=high)
        h=high-0.001;
    return h;
}

int rand_int(int low,int high)
//�����������
{
    return int((high-low+1)*URAND)+low;
}
/*
����������qsort
void qsort( void *base, size_t num, size_t width, int (__cdecl *compare )
����qsort��F[i]���鰴��cmp3����
*/
int cmp1(const void *a,const void *b)
//Ŀ�꺯��f1����������
{
    const individual *e=(const individual *)a;
    const individual *f=(const individual *)b;
    if(e->fvalue[0]==f->fvalue[0])
        return 0;
    else if(e->fvalue[0]<f->fvalue[0])
        return -1;
    else return 1;
}

int cmp2(const void *a,const void *b)
//Ŀ�꺯��f2����������
{
    const individual *e=(const individual *)a;
    const individual *f=(const individual *)b;
    if(e->fvalue[1]==f->fvalue[1])
        return 0;
    else if(e->fvalue[1]<f->fvalue[1])
        return -1;
    else return 1;
}
int cmp_c_d(const void *a,const void *b)
//��ӵ�����뽵������
{
    const individual *e=(const individual *)a;
    const individual *f=(const individual *)b;
    if(e->crowding_distance==f->crowding_distance)
        return 0;
    else if(e->crowding_distance<f->crowding_distance)
        return 1;
    else
        return -1;
}

void population::f_sort(int i)
{
 int n;
 n=len[i];
 qsort(F[i],n,sizeof(individual),cmp_c_d);
}

population::population()
{
    int i;
    for(i=0;i<popsize;i++)
    {
        P[i].init();
    }
    for(i=0;i<popsize;i++)
    {
        P[i].f_count();
    }
    Pnum=popsize;
    Qnum=0;
    Rnum=0;
}
void individual::init()
{
    for(int i=0;i<Dimension;i++)
        value[i]=rand_real(0.0,1.0);
}
/*
���ö����ƽ����������Ӵ���

1���������һ����ʼ����Po���ڴ˻����ϲ��ö�Ԫ������ѡ�񡢽���ͱ�����������Ӵ�Qo�� Po ��QoȺ���ģ��ΪN
2����Pt��Qt���뵽Rt�У���ʼʱt=0������Rt���п��ٷ�֧������򣬹��������в�ͬ�ȼ��ķ�֧��⼯F1��F2����..
3��������Ҫ����Fi�����и����ӵ�����룬������ӵ���Ƚ����������Pt+1��ֱ��Pt+1��ģΪN��ͼ�е�FiΪF3
*/
void population::make_new_pop()
{
    int i,j,x,y,t1,t2,t3;
    double s,u,b;
    memset(mark,0,sizeof(mark));
    t3=0;
    while(t3<popsize/2)
    {
        while(t1=t2=rand_int(0,popsize-1),mark[t1]);
        while(t1==t2||mark[t2])
        {
            t2=rand_int(0,popsize-1);
        }
        t1=choice(t1,t2);
        temp1[t3++]=t1;
        mark[t1]=1;
    }
    for(i=0;i<popsize;i++)
    {
        s=rand_real(0.0,1.0);
        if(s<=0.9)
        {
            for(j=0;j<Dimension;j++)
            {
                u=rand_real((0.0+1e-6),(1.0-1e-6));
                if(u<=0.5)
                    b=pow(2*u,1.0/21);
                else
                    b=1.0/pow(2*(1-u),1.0/21);
                x=y=rand_int(0,popsize/2-1);
                while(x==y)
                    y=rand_int(0,popsize/2-1);
                Q[i].value[j]=1.0/2*((1-b)*P[temp1[x]].value[j]+(1+b)*P[temp1[y]].value[j]);
                if(Q[i].value[j]<0)
                    Q[i].value[j]=1e-6;
                else if(Q[i].value[j]>1)
                    Q[i].value[j]=1.0-(1e-6);
                if(i+1<popsize)
                {
                    Q[i+1].value[j]=1.0/2*((1+b)*P[temp1[x]].value[j]+(1-b)*P[temp1[y]].value[j]);
                    if(Q[i+1].value[j]<=0)
                        Q[i+1].value[j]=1e-6;
                    else if(Q[i+1].value[j]>1)
                        Q[i+1].value[j]=(1-1e-6);
                }
            }
            i++;
        }
        else
        {
            for(j=0;j<Dimension;j++)
            {
                x=rand_int(0,popsize/2-1);
                u=rand_real(0.0+(1e-6),1.0-(1e-6));
                if(u<0.5)
                    u=pow(2*u,1.0/21)-1;
                else
                    u=1-pow(2*(1-u),1.0/21);
                Q[i].value[j]=P[temp1[x]].value[j]+(1.0-0.0)*u;
                if(Q[i].value[j]<0)
                    Q[i].value[j]=1e-6;
                else if(Q[i].value[j]>1)
                    Q[i].value[j]=1-(1e-6);
            }
        }
    }
    Qnum=popsize;
    for(i=0;i<popsize;i++)
        Q[i].f_count();
}

void population::set_p_q()
{
    Rnum=0;
    Qnum=popsize;
    int i;
    for(i=0;i<Pnum;i++)
        R[Rnum++]=P[i];
    for(i=0;i<Qnum;i++)
        R[Rnum++]=Q[i];
    for(i=0;i<2*popsize;i++)
        R[i].f_count();
}

void population::set_p_q()
{
    Rnum=0;
    Qnum=popsize;
    int i;
    for(i=0;i<Pnum;i++)
        R[Rnum++]=P[i];
    for(i=0;i<Qnum;i++)
        R[Rnum++]=Q[i];
    for(i=0;i<2*popsize;i++)
        R[i].f_count();
}
//ZDT1���⺯��ֵ�ļ���

void individual::f_count()
{
    fvalue[0]=value[0];
    int i;
    double g=1,sum=0;
    for(i=1;i<Dimension;i++)
    {
        sum+=value[i];
    }
    sum=9*(sum/(Dimension-1));
    g+=sum;
    fvalue[1]=g*(1-sqrt(value[0]/g));
}
//�ж�Ŀ�꺯��ֵ�Ƿ�֧��
bool e_is_dominated(const individual &a,const individual &b)
{
    if((a.fvalue[0]<=b.fvalue[0])&&(a.fvalue[1]<=b.fvalue[1]))
    {
        if((a.fvalue[0]==b.fvalue[0])&&a.fvalue[1]==b.fvalue[1])
            return false;
        else
            return true;
    }
    else
        return false;
}
//���ٷ�֧�����򷨣��ص㣡����
void population::fast_nondominated_sort()
{  
    int i,j,k;
    individual H[2*popsize];
    int h_len=0;
    for(i=0;i<2*popsize;i++)
    {
        R[i].np=0;
        R[i].is_dominated=0;
        len[i]=0;
    }
    for(i=0;i<2*popsize;i++)
    {
        for(j=0;j<2*popsize;j++)
        {
            if(i!=j)
            {
                if(e_is_dominated(R[i],R[j]))
                    R[i].sp[R[i].is_dominated++]=j;
                else if(e_is_dominated(R[j],R[i]))
                    R[i].np+=1;
            }
        }
        if(R[i].np==0)
        {
            len_f=1;
            F[0][len[0]++]=R[i];
        }

    }
    i=0;
    while(len[i]!=0)
    {
        h_len=0;
        for(j=0;j<len[i];j++)
        {
            for(k=0;k<F[i][j].is_dominated;k++)
            {
                R[F[i][j].sp[k]].np--;
                if(R[F[i][j].sp[k]].np==0)
                {
                    H[h_len++]=R[F[i][j].sp[k]];
                    R[F[i][j].sp[k]].rank=i+2;
                }
            }
        }
        i++;
        len[i]=h_len;
        if(h_len!=0)
        {
            len_f++;
            for(j=0;j<len[i];j++)
                F[i][j]=H[j];
        }
    }
}
//����ӵ�����룺�ص㣡����������ͼ��������£�����
void population::calu_crowding_distance(int i)
{
    int n=len[i];
    double m_max,m_min;
    int j;
    for(j=0;j<n;j++)
        F[i][j].crowding_distance=0;
    F[i][0].crowding_distance=F[i][n-1].crowding_distance=0xffffff;
    qsort(F[i],n,sizeof(individual),cmp1);
    m_max=-0xfffff;
    m_min=0xfffff;
    for(j=0;j<n;j++)
    {
        if(m_max<F[i][j].fvalue[0])
            m_max=F[i][j].fvalue[0];
        if(m_min>F[i][j].fvalue[0])
            m_min=F[i][j].fvalue[0];
    }
    for(j=1;j<n-1;j++)
        F[i][j].crowding_distance+=(F[i][j+1].fvalue[0]-F[i][j-1].fvalue[0])/(m_max-m_min);
    F[i][0].crowding_distance=F[i][n-1].crowding_distance=0xffffff;
    qsort(F[i],n,sizeof(individual),cmp2);
    m_max=-0xfffff;
    m_min=0xfffff;
    for(j=0;j<n;j++)
    {
        if(m_max<F[i][j].fvalue[1])
            m_max=F[i][j].fvalue[1];
        if(m_min>F[i][j].fvalue[1])
            m_min=F[i][j].fvalue[1];
    }
    for(j=1;j<n-1;j++)
        F[i][j].crowding_distance+=(F[i][j+1].fvalue[1]-F[i][j-1].fvalue[1])/(m_max-m_min);
}
//�ɼ������Ե�ѡ��
int population::choice(int a,int b)
{
    if(P[a].rank<P[b].rank)
        return a;
    else if(P[a].rank==P[b].rank)
    {
        if(P[a].crowding_distance>P[b].crowding_distance)
            return a;
        else
            return b;
    }
    else
        return b;
}
//��Ҫ��������
void population::maincal()
{
    int s,i,j;
    s=generation;
    make_new_pop();
    while(s--)
    {
        printf("The %d generation\n",s);
        set_p_q();
        fast_nondominated_sort();
        Pnum=0;
        i=0;
        while(Pnum+len[i]<=popsize)
        {
            calu_crowding_distance(i);
            for(j=0;j<len[i];j++)
                P[Pnum++]=F[i][j];
            i++;
            if(i>=len_f)break;
        }
        if(i<len_f)
        {
            calu_crowding_distance(i);
            f_sort(i);
        }
        for(j=0;j<popsize-Pnum;j++)
            P[Pnum++]=F[i][j];
        make_new_pop();
    }
}

//������
int main()
{
    FILE *p;
    p=fopen("d:\\My_NSGA2.txt","w+");
    srand((unsigned int)(time(0)));
    population pop;
    pop.maincal();
    int i,j;
    fprintf(p,"XuYi All Rights Reserved.\nWelcome to OmegaXYZ: www.omegaxyz.com\n");
    fprintf(p,"Problem ZDT1\n");
    fprintf(p,"\n");
    for(i=0;i<popsize;i++)
    {
        fprintf(p,"The %d generation situation:\n",i);
        for(j=1;j<=Dimension;j++)
        {
            fprintf(p,"x%d=%e  ",j,pop.P[i].value[j]);
        }
        fprintf(p,"\n");
        fprintf(p,"f1(x)=%f   f2(x)=%f\n",pop.P[i].fvalue[0],pop.P[i].fvalue[1]);
    }
    fclose(p);
    return 1;
}