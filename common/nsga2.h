#pragma once
#include<math.h>
#include<time.h>
#include<iostream>
#include<vector>
#include <algorithm>

#ifndef NSGA2CPP
#define NSGA2CPP

constexpr int  mu = 20; //交叉算法的分布指数(distribution index for crossoverand)
constexpr int  mum = 20; //变异算法的分布指数(distribution index for mutation)
constexpr double pm = 0.1; //变异概率(probability for mutation)0.0001-0.1
constexpr double pc = 1 - pm; //交叉概率(probability for crossover)

#define URAND (rand()/(RAND_MAX+1.0))//产生随机数
namespace NSGA2
{
    using namespace std;
    typedef void* VoidPtr;
    template<typename Gene> class Individual;
    template<typename TIndividual> class population;
    //个体的类声明
    template<typename Gene> class Individual
    {
    protected:
        population<Individual<Gene>>* parent;
    protected:
        void init();//初始化个体
    public:
        vector<Gene> x;//xi的值
    public:
        vector<int> sp;//被支配个体集合SP。该量是(可行解空间中)所有被个体p支配的个体组成的集合。
        int np;//支配个数np。该量是(可行解空间中)可以支配个体p的所有个体的数量。
        int dominatedNum;//集合sp的个数
        int rank;//优先级，Pareto级别为当前最高级
        double crowding_distance;//拥挤距离
    public:
        virtual void eval() = 0;//计算y值
        /// <summary>
        /// 基因约束条件
        /// </summary>
        /// <param name="genNo">基因号</param>
        /// <param name="value">基因值</param>
        virtual void constraints(int genNo, Gene& value) = 0;
        /// <summary>
        /// 判断oth是否被当前个体支配
        /// </summary>
        /// <param name="oth"></param>
        /// <returns></returns>
        virtual bool is_dominated(const Individual* oth) = 0;
        /// <summary>
        /// 随机生成基因
        /// </summary>
        /// <param name="genNo"></param>
        /// <returns></returns>
        virtual Gene RandGene(int genNo) = 0;
        /// <summary>
        /// 基因半值
        /// </summary>
        /// <param name="genNo"></param>
        /// <returns></returns>
        virtual Gene HalfGene(int genNo) = 0;
    public:
        /// <summary>
        /// 基因变异
        /// </summary>
        /// <param name="genNo">基因号</param>
        /// <returns>基因约束值</returns>
        Gene mutation(int genNo);
        /// <summary>
        /// 交叉
        /// </summary>
        /// <param name="genNo">基因号</param>
        /// <param name="oth">其他基因值</param>
        /// <returns>交叉生成的新基因</returns>
        Gene* crossover(int genNo, Gene oth);
        Individual(population<Individual<Gene>>* Parent);
    };
    //群体的类声明
    template<typename TIndividual> class population
    {
    public:
        typedef TIndividual* PIndividual;
    private:
        int popsize_t;
        int generation_t;
        int half_popsize;
        int xnum_t;
        //全局变量及部分函数声明
        std::vector <std::vector<PIndividual>>F;//[2 * popsize] ;
        void clear(vector<PIndividual>& Populate);
    protected:
        std::vector<PIndividual> P;//父代popsize
        std::vector<PIndividual> Q;//子代popsize
        std::vector<PIndividual> R;//合并2 * popsize
        //int* Lengths; //各个变异交叉后的群体Fi的长度的集合 [2 * popsize]
        /// <summary>
        /// 合并两代
        /// </summary>
        int merge_p_q();
        /// <summary>
        /// 产生新的子代,随机产生一个初始父代P，在此基础上采用二元锦标赛选择、
        /// 交叉和变异操作产生子代Q。P和Q群体规模均为popsize
        /// 将Pt和Qt并入到Rt中（初始时t=0），对Rt进行快速非支配解排序，
        /// 构造其所有不同等级的非支配解集F1、F2........
        /// </summary>
        void make_new_pop();//产生新的子代
        void fast_nondominated_sort(int Rnum);//快速非支配排序
        int choice(int a, int b);
        //两个个体属于不同等级的非支配解集，优先考虑等级序号较小的
        //若两个个体属于同一等级的非支配解集，优先考虑拥挤距离较大的
        int rank;//整个群体rank值
        void init();
        void final();
    protected:
        /// <summary>
        /// 调用结束
        /// </summary>
        virtual void maincalled();
        /// <summary>
        /// 交叉
        /// </summary>
        virtual void crossover(std::vector<int>& markNum);
        /// <summary>
        /// 变异
        /// </summary>
        virtual void mutation(std::vector<int>& markNum);
        /// <summary>
        /// 目标排序函数，非常重要！！
        /// </summary>
        /// <param name="Fi"></param>
        /// <param name="count"></param>
        /// <param name="num">函数号 1-f1升序,2-f2升序,n-fn升序,0-拥挤距离降序</param>
        virtual bool f_sort(vector<PIndividual> Fi, int num);//排序
        /// <summary>
        /// 新建个体
        /// </summary>
        /// <returns></returns>
        virtual TIndividual* NewIndividual() = 0;
        /// <summary>
        /// 拥挤距离Dn算法
        /// 1.所有个体的拥挤度清零Dn=0，n∈1…N
        /// 2.遍历每个目标函数fm
        ///    ① 根据该目标函数对该等级的个体进行排序，记fmmax​为个体目标函数值fm​的最大值，fmmin​为个体目标函数值fm​的最小值；
        ///    ② 对于排序后两个边界个体(n=1,n=N)的拥挤度D1​和Dn​置为∞；
        ///    ③ 计算Dn​=Dn​+(fm​(i+1)−fm​(i−1))/(fmmax​−fmmin​),其中fm​(i+1)是该个体排序后后一位的目标函数值。
        /// </summary>
        /// <param name="Fi"></param>
        virtual void calu_crowding_distance(vector<PIndividual> Fi) = 0;
        virtual void onexec(int iternum);
    public:
        void maincall();//主要操作
        population(int popsize_t, int generation_t, int xnum_t);
        /// <summary>
        /// 种群大小
        /// </summary>
        /// <returns></returns>
        int Popsize();
        /// <summary>
        /// 基因维度
        /// </summary>
        /// <returns></returns>
        int Dimension();
    };
    double rand_real(double low, double high);
    int rand_int(int low, int high);//产生随机整数  

    //=========================================================================================================
    // 实现
    // =======================================================================================================
    
    //产生随机实数
    inline double rand_real(double low, double high)
    {
        double h;
        h = (high - low) * URAND + low + 0.001;
        if (h >= high)
            h = high - 0.001;
        return h;
    }
    //产生随机整数   
    inline int rand_int(int low, int high)     
    {
        return int((high - low + 1) * URAND) + low;
    }
    //=========================TIndividual=======================================
    template<typename Gene>
        void Individual<Gene>::init()
    {
        for (int i = 0; i < parent->Dimension(); i++)
            x.push_back(RandGene(i));
        eval();
    }
    template<typename Gene>Gene Individual<Gene>
        ::mutation(int genNo)
    {
        Gene u = RandGene(genNo);
        if (u < HalfGene(genNo))
            u = pow(2 * u, 1.0 / (mum + 1)) - 1;
        else
            u = 1 - pow(2 * (1 - u), 1.0 / (mum + 1));
        Gene value = x[genNo] + u;
        constraints(genNo, value);
        return value;
    }
    template<typename Gene> Gene * Individual<Gene>
        ::crossover(int genNo, Gene oth)
    {
        //默认模拟二进制交叉
        Gene u = RandGene(genNo);
        Gene b;
        if (u <= HalfGene(genNo))
            b = pow(2 * u, 1.0 / (mu + 1));
        else
            b = 1.0 / pow(2 * (1 - u), 1.0 / (mu + 1));
        Gene g1 = x[genNo];
        Gene result[2];
        result[0] = 1.0 / 2 * ((1 + b) * g1 + (1 - b) * oth);
        constraints(genNo, result[0]);
        result[1] = 1.0 / 2 * ((1 - b) * g1 + (1 + b) * oth);
        constraints(genNo, result[1]);
        return result;
    }
    template<typename Gene>Individual<Gene>
        ::Individual(population<Individual<Gene>>* Parent)
    {
        parent = Parent;
        np = 0;
        dominatedNum = 0;
        rank = 0;
        crowding_distance = 0;
    }
    //===============================TIndividual end=============================
    template<typename TIndividual>int population<TIndividual>
        ::merge_p_q()
    {
        clear(R);
        for (auto& p : P)
        {
            R.push_back(p);
            p = nullptr;
        }
        int maxCount=Q.size()<=popsize_t ? Q.size():popsize_t;
        for (int i = 0; i < maxCount; i++)
        {
            R.push_back(Q[i]);
            Q[i] = nullptr;
        }
        clear(P);
        clear(Q);
        return R.size();
    }
    //精英保留策略选择,原则：选择优先级高的(数值小)，相同则选择拥挤度小的(数值大)
    template<typename TIndividual>int population<TIndividual>
        ::choice(int a, int b)
    {
        if (P[a]->rank < P[b]->rank)
            return a;
        else if (P[a]->rank == P[b]->rank)
        {
            if (P[a]->crowding_distance > P[b]->crowding_distance)
                return a;
            else
                return b;
        }
        else
            return b;
    }
    template<typename TIndividual>void population<TIndividual>
        ::crossover(std::vector<int>& markNum)
    {
        PIndividual Q0 = NewIndividual();
        PIndividual Q1 = NewIndividual();
        Q.push_back(Q0);
        Q.push_back(Q1);
        for (int j = 0; j < xnum_t; j++)
        {
            int mark0, mark1;
            mark0 = mark1 = rand_int(0, half_popsize - 1);
            while (mark0 == mark1)
                mark1 = rand_int(0, half_popsize - 1);
            mark0 = markNum[mark0];
            mark1 = markNum[mark1];
            //默认模拟二进制交叉
            auto cross = P[mark0]->crossover(j, P[mark1]->x[j]);
            Q0->x[j] = cross[0];
            Q1->x[j] = cross[1];
        }
        Q0->eval();
        Q1->eval();
    }
    template<typename TIndividual>void population<TIndividual>
        ::mutation(std::vector<int> & markNum)
    {
        //多项式变异
        PIndividual Qi = NewIndividual();
        Q.push_back(Qi);
        for (int j = 0; j < xnum_t; j++)
        {
            int mark0 = markNum[rand_int(0, half_popsize - 1)];
            Qi->x[j] = P[mark0]->mutation(j);
        }
        Qi->eval();
    }
    /*
    利用二进制锦标赛产生子代：
    1、随机产生一个初始父代Po，在此基础上采用二元锦标赛选择、交叉和变异操作产生子代Qo， Po 和Qo群体规模均为N
    2、将Pt和Qt并入到Rt中（初始时t=0），对Rt进行快速非支配解排序，构造其所有不同等级的非支配解集F1、F2……..
    3、按照需要计算Fi中所有个体的拥挤距离，并根据拥挤比较运算符构造Pt+1，直至Pt+1规模为N，图中的Fi为F3

    实数编码的交叉操作（SBX）
        模拟二进制交叉：
            x1j​(t)=0.5×[(1+γj​)x1j​(t)+(1−γj​)x2j​(t)]
            x2j​(t)=0.5×[(1−γj​)x1j​(t)+(1+γj​)x2j​(t)]
           其中：
            γj​=(2uj​)^1/(η+1),当uj​<0.5
            γj​=(1/(2(1−uj​))^1/(η+1)​,当uj>=0.5
    多项式变异（polynomial mutation）
        x1j​(t)=x1j​(t)+∆j​
        其中：
            ∆j​=(2uj​)^1/(η+1)−1,当uj​<0.5
            ∆j​=1−(2(1−uj​))^1/(η+1)​,当uj>=0.5
    */
    template<typename TIndividual> void population<TIndividual>::make_new_pop()
    {
        //两个用于产生新子代的数组
        std::vector<bool> mark(popsize_t);//标记数组
        std::vector<int> markNum(half_popsize);//标记序号
        int t3 = 0;
        //二元锦标赛选择
        while (t3 < half_popsize)
        {
            int t1=0;
            int t2=0;
            while (t1 = t2 = rand_int(0, popsize_t - 1), mark[t1]);
            while (t1 == t2 || mark[t2])
            {
                t2 = rand_int(0, popsize_t - 1);
            }
            t1 = choice(t1, t2);
            markNum[t3++] = t1;
            mark[t1] = true;
        }
        clear(Q);
        for (int i = 0; i < half_popsize; i++)
        {
            double s = rand_real(0.0, 1.0);//交叉概率值大于pc时变异，否则交叉
            if (s > pc)
            {
                mutation(markNum);
            }
            else 
            {
                crossover(markNum);
            }
        }
    }
    template<typename TIndividual> void population<TIndividual>
        ::clear(vector<PIndividual>& Populate)
    {
        for (auto& pop : Populate)
        {
            if(pop!=nullptr)
                delete pop;
        }
        Populate.clear();
    }
    //快速非支配排序法：重点！！！
    template<typename TIndividual> void population<TIndividual>
        ::fast_nondominated_sort(int Rnum)
        {
            vector<PIndividual> H;
            for (auto& f : F)
            {
                f.clear();
            }
            F.clear();
            F.resize(Rnum);
            int Len = 0;
            vector<bool> Remove(Rnum);
            for (int i = 0; i < Rnum; i++)
            {
                R[i]->np = 0;//支配个数np
                R[i]->dominatedNum = 0;//被支配的个数
                R[i]->sp.clear();
            }
            for (int i = 0; i < Rnum; i++)
            {
                for (int j = 0; j < Rnum; j++)
                {
                    if (i != j)//自己不能支配自身
                    {
                        if (R[i]->is_dominated(R[j]))
                        {
                            R[i]->dominatedNum++;
                            R[i]->sp.push_back(j);//如果i支配j，把j添加到i的被支配列表中
                        }
                        else if (R[j]->is_dominated(R[i]))
                            R[i]->np += 1;//如果i被j支配，则把np加1
                    };
                }
                if (R[i]->np == 0)//如果该个体的np为0，则该个体为Pareto第一级
                {
                    //将R[i]归并
                    rank = 1;
                    F[0].push_back(R[i]);
                    Remove[i] = true;
                    Len++;
                }
            }
            int i = 0;
            while (Len != 0)
            {
                H.clear();
                for (int j = 0; j < Len; j++)
                {
                    for (int k = 0; k < F[i][j]->dominatedNum; k++)//对所有在dominatedNum集合中的个体进行排序
                    {
                        int index = F[i][j]->sp[k];
                        R[index]->np--;
                        if (R[index]->np == 0) //如果该个体的支配个数为0，则该个体是非支配个体
                        {
                            R[index]->rank = i + 2;
                            H.push_back(R[index]);
                            Remove[index] = true;
                        }
                    }
                }
                i++;
                Len = H.size();
                if (Len > 0)
                {
                    rank++;                    
                    for (const auto& h : H)
                        F[i].push_back(h);
                }
            }
            for (int j = 0; j < Rnum; j++)
            {
                if (!Remove[j])
                    delete R[j];
                R[j] = nullptr;
            }
            R.clear();
        }        
    template<typename TIndividual> inline bool cmp_c_d(void const* a, void const* b)
        //对拥挤距离降序排序
    {
        const TIndividual* e = (const TIndividual*)a;
        const TIndividual* f = (const TIndividual*)b;
        return e->crowding_distance < f->crowding_distance;
    }
     template<typename TIndividual> bool population<TIndividual>
         ::f_sort(vector<PIndividual> Fi, int num)
     {
            bool result = false;
            switch (num)
            {
            case 0:
                sort(Fi.begin(), Fi.end(), cmp_c_d<TIndividual>);
            default:
                result = false; 
                break;
            }
            return result;
     }
    //主要操作函数
    template<typename TIndividual> void population<TIndividual>::maincall()
    {
        init();
        try
        {
            int gen;
            gen = generation_t;
            make_new_pop();
            while (gen--)
            {
                onexec(generation_t-gen);
                int Rnum=merge_p_q();
                fast_nondominated_sort(Rnum);
                int Pnum = 0;
                int i = 0;
                while (Pnum + F[i].size() <= popsize_t)
                {
                    calu_crowding_distance(F[i]);
                    for (int j = 0; j < F[i].size(); j++)
                    {
                        Pnum++;
                        P.push_back(F[i][j]);
                    }
                    F[i].clear();
                    i++;
                    if (i >= rank)break;
                }
                if (i < rank)
                {
                    calu_crowding_distance(F[i]);
                    f_sort(F[i],0);
                }
                int count = popsize_t - Pnum;
                for (int j = 0; j < count; j++)
                {
                    P.push_back(F[i][j]);
                    Pnum++;
                }
                for (int j = count; j < F[i].size(); j++)
                    delete F[i][j];
                F[i].clear();

                for (int j = i + 1; j < F.size(); j++)
                {
                    clear(F[j]);
                }
                make_new_pop();
            }
            maincalled();
        }
        catch(...)
        {
        }
        final();
    }
    template<typename TIndividual> void population<TIndividual>::init()
    {
        for (int i = 0; i < popsize_t; i++)
        {
            P.push_back(NewIndividual());
        }
    }
    template<typename TIndividual> void population<TIndividual>::final()
    {
        clear(P);
        clear(Q);
        clear(R);
        for (auto& f : F)
        {
            f.clear();
        }
        F.clear();
    }
    template<typename TIndividual>int population<TIndividual>::Popsize()
    {
        return popsize_t;
    }
    template<typename TIndividual>int population<TIndividual>::Dimension()
    {
        return xnum_t;
    }
    template<typename TIndividual> population<TIndividual>
        ::population(int popsize_t, int generation_t, int xnum_t)
    {
        srand((unsigned int)(time(0)));
        this->generation_t = generation_t;
        this->popsize_t = popsize_t;
        this->xnum_t = xnum_t;
        half_popsize = popsize_t / 2;
    }
    template<typename TIndividual> void population<TIndividual>::maincalled()
    {}
    template<typename TIndividual> void population<TIndividual>::onexec(int iternum)
    {}
}
#endif