#pragma once
#include "../../common/nsga2.h"
#include<iostream>

constexpr int ynum = 2; //目标函数数量,多目标函数数量

constexpr double min_range = 1e-6; //基因下界
constexpr double max_range = 1.0 - min_range; //基因上界
constexpr int generation = 500; //繁衍代数
constexpr int popsize = 500;  //种群大小(population size)
constexpr int xnum = 30;//基因维度Dimension，在这里即ZDT1问题xi的i的最大值(自变量数量)

namespace NSGA2
{
    using namespace std;
    class populationOfZDT1;
    //个体的类声明
    class IndividualOfZDT1
        :public Individual<double>
    {
    public:
        double y[ynum];//目标函数的值
        void eval() override;//计算y值
        void constraints(int genNo, double& g) override;
        bool is_dominated(const Individual* oth) override;
        double RandGene(int genNo) override;
        double HalfGene(int genNo) override;
    public:
        IndividualOfZDT1(populationOfZDT1 * Parent);
    };
    typedef IndividualOfZDT1* pindividual;

    //群体的类声明
    class populationOfZDT1 :public population<IndividualOfZDT1>
    {
    private:
        void print(string fileName = "My_NSGA2.txt");//打印结果到文件
        void exportCsv(string fileName = "NSGA2.csv");
    protected:
        /// <summary>
        /// 调用结束
        /// </summary>
        void maincalled() override;
        /// <summary>
        /// 排序
        /// </summary>
        /// <param name="Fi"></param>
        /// <param name="num">函数号 1-f1升序,2-f2升序,n-fn升序,0-拥挤距离降序</param>
        bool f_sort(vector<PIndividual> Fi, int num) override;//排序
        IndividualOfZDT1 * NewIndividual() override;
        void onexec(int iternum) override;
        void calu_crowding_distance(vector<PIndividual> Fi) override;
    public:
        populationOfZDT1();
    };
}
