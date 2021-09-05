#pragma once
#include "../../common/nsga2.h"
#include<iostream>

constexpr int ynum = 2; //Ŀ�꺯������,��Ŀ�꺯������

constexpr double min_range = 1e-6; //�����½�
constexpr double max_range = 1.0 - min_range; //�����Ͻ�
constexpr int generation = 500; //���ܴ���
constexpr int popsize = 500;  //��Ⱥ��С(population size)
constexpr int xnum = 30;//����ά��Dimension�������ＴZDT1����xi��i�����ֵ(�Ա�������)

namespace NSGA2
{
    using namespace std;
    class populationOfZDT1;
    //�����������
    class IndividualOfZDT1
        :public Individual<double>
    {
    public:
        double y[ynum];//Ŀ�꺯����ֵ
        void eval() override;//����yֵ
        void constraints(int genNo, double& g) override;
        bool is_dominated(const Individual* oth) override;
        double RandGene(int genNo) override;
        double HalfGene(int genNo) override;
    public:
        IndividualOfZDT1(populationOfZDT1 * Parent);
    };
    typedef IndividualOfZDT1* pindividual;

    //Ⱥ���������
    class populationOfZDT1 :public population<IndividualOfZDT1>
    {
    private:
        void print(string fileName = "My_NSGA2.txt");//��ӡ������ļ�
        void exportCsv(string fileName = "NSGA2.csv");
    protected:
        /// <summary>
        /// ���ý���
        /// </summary>
        void maincalled() override;
        /// <summary>
        /// ����
        /// </summary>
        /// <param name="Fi"></param>
        /// <param name="num">������ 1-f1����,2-f2����,n-fn����,0-ӵ�����뽵��</param>
        bool f_sort(vector<PIndividual> Fi, int num) override;//����
        IndividualOfZDT1 * NewIndividual() override;
        void onexec(int iternum) override;
        void calu_crowding_distance(vector<PIndividual> Fi) override;
    public:
        populationOfZDT1();
    };
}
