#include<math.h>
#include<time.h>
#include<iostream>
#define Dimension 30//基因维数，在这里即ZDT1问题xi的i的最大值
#define popsize 100//种群大小
#define mu 20 //
#define mum 20 //

#define min_range 0.0 //下界
#define max_range 1.0 //上界
#define half_popsize popsize/2
#define generation 500 //繁衍代数
#define URAND (rand()/(RAND_MAX+1.0))//产生随机数

using namespace std;
namespace NSGA2
{
    //个体的类声明
    class individual
    {
       public:
        double value[Dimension];//xi的值
        int sp[2 * popsize];//被支配个体集合SP。该量是可行解空间中所有被个体p支配的个体组成的集合。
        int np;//支配个数np。该量是在可行解空间中可以支配个体p的所以个体的数量。
        int is_dominated;//集合sp的个数
        void init();//初始化个体
        int rank;//优先级，Pareto级别为当前最高级
        double crowding_distance;//拥挤距离
        double fvalue[2];//ZDT1问题目标函数的值
        void f_count();//计算fvalue的值
    };
    //群体的类声明
    class population
    {
       private:
        //全局变量及部分函数声明
        individual F[2 * popsize][2 * popsize];
       protected:
        individual P[popsize];//父代
        individual Q[popsize];//子代
        individual R[2 * popsize];//合并
        void set_p_q();
        //随机产生一个初始父代P，在此基础上采用二元锦标赛选择、
        //交叉和变异操作产生子代Q。P和Q群体规模均为popsize
        //将Pt和Qt并入到Rt中（初始时t=0），对Rt进行快速非支配解排序，
        //构造其所有不同等级的非支配解集F1、F2........
        int Rnum;
        int Pnum;
        int Qnum;
        //P,Q,R中元素的个数
        void make_new_pop();//产生新的子代
        void fast_nondominated_sort();//快速非支配排序
        void calu_crowding_distance(int i);//拥挤距离计算
        void f_sort(int i);//对拥挤距离降序排列
        void print(string fileName = "My_NSGA2.txt");//打印结果到文件
        void exportCsv(string fileName = "NSGA2.csv");
        int choice(int a, int b);
        //两个个体属于不同等级的非支配解集，优先考虑等级序号较小的
        //若两个个体属于同一等级的非支配解集，优先考虑拥挤距离较大的
        int len[2 * popsize];//各个变异交叉后的群体Fi的长度的集合
        int len_f;//整个群体rank值
       public:
        population();//类初始化
        void maincal();//主要操作
    };
    /// <summary>
    /// 产生随机实数
    /// </summary>
    /// <param name="low"></param>
    /// <param name="high"></param>
    /// <returns></returns>
    double rand_real(double low, double high);
    /// <summary>
    /// 产生随机整数
    /// </summary>
    /// <param name="low"></param>
    /// <param name="high"></param>
    /// <returns></returns>
    int rand_int(int low, int high);
}