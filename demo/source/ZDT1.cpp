#pragma once
#include "../include/ZDT1.h"
#include <algorithm>
namespace NSGA2
{
    //=========================IndividualOfZDT1=======================================
    //ZDT1问题函数值的计算
    void IndividualOfZDT1::eval()
    {
        y[0] = x[0];
        int i;
        double g = 1, sum = 0;
        for (i = 1; i < parent->Dimension(); i++)
        {
            sum += x[i];
        }
        sum = 9 * (sum / (parent->Dimension() - 1));
        g += sum;
        y[1] = g * (1 - sqrt(x[0] / g));
    }
    void  IndividualOfZDT1::constraints(int genNo, double & g)
    {
        //当前各基因约束条件一致，所以基因号未使用
        if (g < min_range)
            g = min_range;
        else if (g > max_range)
            g = max_range;
    }
    //判断目标函数值是否被支配    
    bool  IndividualOfZDT1::is_dominated(const Individual* oth)
    {
        const IndividualOfZDT1* e = (const IndividualOfZDT1*)oth;
        if ((y[0] <= e->y[0]) && (y[1] <= e->y[1]))
        {
            if ((y[0] == e->y[0]) && y[1] == e->y[1])
                return false;
            else
                return true;
        }
        else
            return false;
    }
    double IndividualOfZDT1::RandGene(int genNo)
    {
       return rand_real(min_range, max_range);
    }
    double IndividualOfZDT1::HalfGene(int genNo)
    {
        return (max_range - min_range) / 2;
    }
    IndividualOfZDT1::IndividualOfZDT1(populationOfZDT1* Parent)
        : Individual(reinterpret_cast<population<Individual<double>>*>(Parent))
    {
        init();
    }
    //===============================IndividualOfZDT1 end==================================
    /*
    关于排序函数qsort
    void qsort( void *base, size_t num, size_t width, int (__cdecl *compare))
    利用qsort对F[i]数组按照cmp3排序
    */
    inline bool cmp1(const void* a, const void* b)
        //目标函数f1的升序排序
    {
        const IndividualOfZDT1* e = (const IndividualOfZDT1*)a;
        const IndividualOfZDT1* f = (const IndividualOfZDT1*)b;
        return e->y[0] >= f->y[0];
    }
    inline bool cmp2(const void* a, const void* b)
        //目标函数f2的升序排序
    {
        const IndividualOfZDT1* e = (const IndividualOfZDT1*)a;
        const IndividualOfZDT1* f = (const IndividualOfZDT1*)b;
        return e->y[1] >= f->y[1];
    }
    /// <summary>
    /// 排序
    /// </summary>
    /// <param name="Fi"></param>
    /// <param name="count"></param>
    /// <param name="num">函数号 0-f1升序,1-f2升序,2-拥挤距离降序</param>
    bool populationOfZDT1::f_sort(vector<PIndividual> Fi, int num)
    {
        bool result = population<IndividualOfZDT1>::f_sort(Fi,num);
        if (!result)
        {
            switch (num)
            {
            case 1:
                sort(Fi.begin(), Fi.end(), cmp1);
                result = true;
                break;
            case 2:
                sort(Fi.begin(), Fi.end(), cmp2);
                result = true;
                break;
            }
        }
        return result;
    }
    IndividualOfZDT1* populationOfZDT1::NewIndividual()
    {
        return new IndividualOfZDT1(this);
    }
    void populationOfZDT1::calu_crowding_distance(vector<PIndividual> Fi)
        {
            const float MaxD = 0xffffff;
            const float MinD = -0xffffff;
            int n = Fi.size();
            //所有个体的拥挤度清零
            for (const auto& fj : Fi)
            {
                fj->crowding_distance = 0;
            }
            /*根据目标函数计算个体拥挤度*/
            for (int i = 0; i < ynum; i++)
            {
                f_sort(Fi, i + 1);//根据目标函数排序
                Fi[0]->crowding_distance = Fi[n - 1]->crowding_distance = MaxD;//两个边界个体拥挤度​置为∞        
                double fmax = MinD;
                double fmin = MaxD;
                for (int j = 0; j < n; j++) //查找f1最大最小值
                {
                    if (fmax < Fi[j]->y[i])
                        fmax = Fi[j]->y[i];
                    if (fmin > Fi[j]->y[i])
                        fmin = Fi[j]->y[i];
                }
                for (int j = 1; j < n - 1; j++)//计算个体拥挤度
                    Fi[j]->crowding_distance += (Fi[j + 1]->y[i] - Fi[j - 1]->y[i]) / (fmax - fmin);
            }
        }
    populationOfZDT1::populationOfZDT1()
        : population<IndividualOfZDT1>(popsize, generation, xnum)
    {
    }
    void populationOfZDT1::onexec(int iternum)
    {
        if ((iternum % 100) == 0)
            printf("The %d generation\n", iternum);
    }
    /// <summary>
    /// 调用结束
    /// </summary>
    void populationOfZDT1::maincalled()
    {
        string pFileName = "My_NSGA2.txt";
        string eFileName = "NSGA2.csv";
        string mainFile = __argv[0];
        string directory;
        const size_t last_slash_idx = mainFile.rfind('\\');
        if (std::string::npos != last_slash_idx)
        {
            directory = mainFile.substr(0, last_slash_idx+1);
            pFileName.insert(0,directory);
            eFileName.insert(0, directory);
        }
        print(pFileName);
        exportCsv(eFileName);
    }
    void populationOfZDT1::print(string fileName)
    {
        FILE* p;
        errno_t err_no = fopen_s(&p, fileName.c_str(), "w+");
        fprintf(p, "XuYi All Rights Reserved.\nWelcome to OmegaXYZ: www.omegaxyz.com\n");
        fprintf(p, "Problem ZDT1\n");
        fprintf(p, "\n");
        for (int i = 0; i < Popsize(); i++)
        {
            fprintf(p, "The %d generation situation:\n", i);
            for (int j = 1; j < Dimension(); j++)
            {
                fprintf(p, "x%d=%e  ", j, P[i]->x[j]);
            }
            fprintf(p, "\n");
            double f2c = 1 - sqrt(P[i]->y[0]);
            fprintf(p, "f1(x)=%0.3f   f2(x)=%0.3f   f2c(x)=%0.3f\n", P[i]->y[0], P[i]->y[1], f2c);
        }
        fclose(p);
        cout << "优化结果保存在[" << fileName << "]，优化结束" << endl;
    }
    void populationOfZDT1::exportCsv(string fileName)
    {
        FILE* p;
        errno_t err_no = fopen_s(&p, fileName.c_str(), "w+");
        int i;
        fprintf(p, "XuYi All Rights Reserved.\nWelcome to OmegaXYZ: www.omegaxyz.com\n");
        fprintf(p, "Problem ZDT1\n");
        fprintf(p, "F1,F2,F2C\n");
        for (i = 0; i < Popsize(); i++)
        {
            double f2c = 1 - sqrt(P[i]->y[0]);
            fprintf(p, "%0.3f,%0.3f,%0.3f\n", P[i]->y[0], P[i]->y[1], f2c);
        }
        fclose(p);
        cout << "优化结果导出到[" << fileName << "]，优化结束" << endl;
    }
}