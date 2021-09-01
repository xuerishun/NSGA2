#include<math.h>
#include<time.h>
#include<iostream>
#define Dimension 30//����ά���������ＴZDT1����xi��i�����ֵ
#define popsize 100//��Ⱥ��С
#define mu 20 //
#define mum 20 //

#define min_range 0.0 //�½�
#define max_range 1.0 //�Ͻ�
#define half_popsize popsize/2
#define generation 500 //���ܴ���
#define URAND (rand()/(RAND_MAX+1.0))//���������

using namespace std;
namespace NSGA2
{
    //�����������
    class individual
    {
       public:
        double value[Dimension];//xi��ֵ
        int sp[2 * popsize];//��֧����弯��SP�������ǿ��н�ռ������б�����p֧��ĸ�����ɵļ��ϡ�
        int np;//֧�����np���������ڿ��н�ռ��п���֧�����p�����Ը����������
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
       private:
        //ȫ�ֱ��������ֺ�������
        individual F[2 * popsize][2 * popsize];
       protected:
        individual P[popsize];//����
        individual Q[popsize];//�Ӵ�
        individual R[2 * popsize];//�ϲ�
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
        void print(string fileName = "My_NSGA2.txt");//��ӡ������ļ�
        void exportCsv(string fileName = "NSGA2.csv");
        int choice(int a, int b);
        //�����������ڲ�ͬ�ȼ��ķ�֧��⼯�����ȿ��ǵȼ���Ž�С��
        //��������������ͬһ�ȼ��ķ�֧��⼯�����ȿ���ӵ������ϴ��
        int len[2 * popsize];//�������콻����Ⱥ��Fi�ĳ��ȵļ���
        int len_f;//����Ⱥ��rankֵ
       public:
        population();//���ʼ��
        void maincal();//��Ҫ����
    };
    /// <summary>
    /// �������ʵ��
    /// </summary>
    /// <param name="low"></param>
    /// <param name="high"></param>
    /// <returns></returns>
    double rand_real(double low, double high);
    /// <summary>
    /// �����������
    /// </summary>
    /// <param name="low"></param>
    /// <param name="high"></param>
    /// <returns></returns>
    int rand_int(int low, int high);
}