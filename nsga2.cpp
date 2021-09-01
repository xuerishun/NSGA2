#include "../include/nsga2.h"

namespace NSGA2
{
    double rand_real(double low, double high)
        //产生随机实数
    {
        double h;
        h = (high - low) * URAND + low + 0.001;
        if (h >= high)
            h = high - 0.001;
        return h;
    }
    int rand_int(int low, int high)
        //产生随机整数
    {
        return int((high - low + 1) * URAND) + low;
    }
    /*
    关于排序函数qsort
    void qsort( void *base, size_t num, size_t width, int (__cdecl *compare )
    利用qsort对F[i]数组按照cmp3排序
    */
    inline int cmp1(const void* a, const void* b)
        //目标函数f1的升序排序
    {
        const individual* e = (const individual*)a;
        const individual* f = (const individual*)b;
        if (e->fvalue[0] == f->fvalue[0])
            return 0;
        else if (e->fvalue[0] < f->fvalue[0])
            return -1;
        else return 1;
    }
    inline int cmp2(const void* a, const void* b)
        //目标函数f2的升序排序
    {
        const individual* e = (const individual*)a;
        const individual* f = (const individual*)b;
        if (e->fvalue[1] == f->fvalue[1])
            return 0;
        else if (e->fvalue[1] < f->fvalue[1])
            return -1;
        else return 1;
    }
    inline int cmp_c_d(const void* a, const void* b)
        //对拥挤距离降序排序
    {
        const individual* e = (const individual*)a;
        const individual* f = (const individual*)b;
        if (e->crowding_distance == f->crowding_distance)
            return 0;
        else if (e->crowding_distance < f->crowding_distance)
            return 1;
        else
            return -1;
    }
    //判断目标函数值是否被支配
    inline bool e_is_dominated(const individual& a, const individual& b)
    {
        if ((a.fvalue[0] <= b.fvalue[0]) && (a.fvalue[1] <= b.fvalue[1]))
        {
            if ((a.fvalue[0] == b.fvalue[0]) && a.fvalue[1] == b.fvalue[1])
                return false;
            else
                return true;
        }
        else
            return false;
    }
    void individual::init()
    {
        for (int i = 0; i < Dimension; i++)
            value[i] = rand_real(0.0, 1.0);
    }
    //ZDT1问题函数值的计算
    void individual::f_count()
    {
        fvalue[0] = value[0];
        int i;
        double g = 1, sum = 0;
        for (i = 1; i < Dimension; i++)
        {
            sum += value[i];
        }
        sum = 9 * (sum / (Dimension - 1));
        g += sum;
        fvalue[1] = g * (1 - sqrt(value[0] / g));
    }
    void population::set_p_q()
    {
        Rnum = 0;
        Qnum = popsize;
        int i;
        for (i = 0; i < Pnum; i++)
            R[Rnum++] = P[i];
        for (i = 0; i < Qnum; i++)
            R[Rnum++] = Q[i];
        for (i = 0; i < 2 * popsize; i++)
            R[i].f_count();
    }
    void population::f_sort(int i)
    {
        int n;
        n = len[i];
        qsort(F[i], n, sizeof(individual), cmp_c_d);
    }
    /*
    利用二进制锦标赛产生子代：
    1、随机产生一个初始父代Po，在此基础上采用二元锦标赛选择、交叉和变异操作产生子代Qo， Po 和Qo群体规模均为N
    2、将Pt和Qt并入到Rt中（初始时t=0），对Rt进行快速非支配解排序，构造其所有不同等级的非支配解集F1、F2……..
    3、按照需要计算Fi中所有个体的拥挤距离，并根据拥挤比较运算符构造Pt+1，直至Pt+1规模为N，图中的Fi为F3
    */
    void population::make_new_pop()
    {
        int i, j, x, y, t1, t2, t3;
        double s, u, b;
        //两个用于产生新子代的数组
        int mark[popsize];//标记数组
        int markNum[half_popsize];//标记序号
        memset(mark, 0, sizeof(mark));
        t3 = 0;
        while (t3 < half_popsize)
        {
            while (t1 = t2 = rand_int(0, popsize - 1), mark[t1]);
            while (t1 == t2 || mark[t2])
            {
                t2 = rand_int(0, popsize - 1);
            }
            t1 = choice(t1, t2);
            markNum[t3++] = t1;
            mark[t1] = 1;
        }
        for (i = 0; i < popsize; i++)
        {
            s = rand_real(0.0, 1.0);
            if (s <= 0.9)
            {
                for (j = 0; j < Dimension; j++)
                {
                    u = rand_real((0.0 + 1e-6), (1.0 - 1e-6));
                    if (u <= 0.5)
                        b = pow(2 * u, 1.0 / (mu+1));
                    else
                        b = 1.0 / pow(2 * (1 - u), 1.0 / (mu + 1));
                    x = y = rand_int(0, half_popsize - 1);
                    while (x == y)
                        y = rand_int(0, half_popsize - 1);
                    //Q[i].value[j] = 1.0 / 2 * ((1 - b) * P[markNum[x]].value[j] + (1 + b) * P[markNum[y]].value[j]);
                    Q[i].value[j] = 1.0 / 2 * ((1 + b) * P[markNum[x]].value[j] + (1 - b) * P[markNum[y]].value[j]);
                    if (Q[i].value[j] < min_range)
                        Q[i].value[j] = min_range + 1e-6;
                    else if (Q[i].value[j] > max_range)
                        Q[i].value[j] = max_range - (1e-6);
                    if (i + 1 < popsize)
                    {
                        //Q[i + 1].value[j] = 1.0 / 2 * ((1 + b) * P[markNum[x]].value[j] + (1 - b) * P[markNum[y]].value[j]);
                        Q[i + 1].value[j] = 1.0 / 2 * ((1 - b) * P[markNum[x]].value[j] + (1 + b) * P[markNum[y]].value[j]);
                        if (Q[i + 1].value[j] <= min_range )
                            Q[i + 1].value[j] = min_range + 1e-6;
                        else if (Q[i + 1].value[j] > max_range)
                            Q[i + 1].value[j] = (max_range - 1e-6);
                    }
                }
                i++;
            }
            else
            {
                for (j = 0; j < Dimension; j++)
                {
                    x = rand_int(0, half_popsize - 1);
                    u = rand_real(0.0 + (1e-6), 1.0 - (1e-6));
                    if (u < 0.5)
                        u = pow(2 * u, 1.0 / 21) - 1;
                    else
                        u = 1 - pow(2 * (1 - u), 1.0 / 21);
                    Q[i].value[j] = P[markNum[x]].value[j] + (1.0 - 0.0) * u;
                    if (Q[i].value[j] < 0)
                        Q[i].value[j] = 1e-6;
                    else if (Q[i].value[j] > 1)
                        Q[i].value[j] = 1 - (1e-6);
                }
            }
        }
        Qnum = popsize;
        for (i = 0; i < popsize; i++)
            Q[i].f_count();
    }
    //快速非支配排序法：重点！！！
    void population::fast_nondominated_sort()
    {
        int i, j, k;
        individual H[2 * popsize];
        int h_len = 0;
        for (i = 0; i < 2 * popsize; i++)
        {
            R[i].np = 0;//支配个数np
            R[i].is_dominated = 0;//被支配的个数
            len[i] = 0;//初始化
        }
        for (i = 0; i < 2 * popsize; i++)
        {
            for (j = 0; j < 2 * popsize; j++)
            {
                if (i != j)//自己不能支配自身
                {
                    if (e_is_dominated(R[i], R[j]))
                        R[i].sp[R[i].is_dominated++] = j;//如果i支配j，把i添加到j的is_domied列表中
                    else if (e_is_dominated(R[j], R[i]))
                        R[i].np += 1;//如果i被j支配，则把np加1
                }
            }
            if (R[i].np == 0)//如果该个体的np为0，则该个体为Pareto第一级
            {
                len_f = 1;
                F[0][len[0]++] = R[i];//将R[i]归并
            }

        }
        i = 0;
        while (len[i] != 0)
        {
            h_len = 0;
            for (j = 0; j < len[i]; j++)
            {
                for (k = 0; k < F[i][j].is_dominated; k++)//对所有在is_dominated集合中的个体进行排序
                {
                    R[F[i][j].sp[k]].np--;
                    if (R[F[i][j].sp[k]].np == 0) //如果该个体的支配个数为0，则该个体是非支配个体
                    {
                        H[h_len++] = R[F[i][j].sp[k]];
                        R[F[i][j].sp[k]].rank = i + 2;
                    }
                }
            }
            i++;
            len[i] = h_len;
            if (h_len != 0)
            {
                len_f++;
                for (j = 0; j < len[i]; j++)
                    F[i][j] = H[j];
            }
        }
    }
    //计算拥挤距离：重点！！！具体解释见其他文章！！！
    void population::calu_crowding_distance(int i)
    {
        int n = len[i];
        double m_max, m_min;
        int j;
        for (j = 0; j < n; j++)
            F[i][j].crowding_distance = 0;
        F[i][0].crowding_distance = F[i][n - 1].crowding_distance = 0xffffff;
        qsort(F[i], n, sizeof(individual), cmp1);
        m_max = -0xfffff;
        m_min = 0xfffff;
        for (j = 0; j < n; j++)
        {
            if (m_max < F[i][j].fvalue[0])
                m_max = F[i][j].fvalue[0];
            if (m_min > F[i][j].fvalue[0])
                m_min = F[i][j].fvalue[0];
        }
        for (j = 1; j < n - 1; j++)
            F[i][j].crowding_distance += (F[i][j + 1].fvalue[0] - F[i][j - 1].fvalue[0]) / (m_max - m_min);
        F[i][0].crowding_distance = F[i][n - 1].crowding_distance = 0xffffff;
        qsort(F[i], n, sizeof(individual), cmp2);
        m_max = -0xfffff;
        m_min = 0xfffff;
        for (j = 0; j < n; j++)
        {
            if (m_max < F[i][j].fvalue[1])
                m_max = F[i][j].fvalue[1];
            if (m_min > F[i][j].fvalue[1])
                m_min = F[i][j].fvalue[1];
        }
        for (j = 1; j < n - 1; j++)
            F[i][j].crowding_distance += (F[i][j + 1].fvalue[1] - F[i][j - 1].fvalue[1]) / (m_max - m_min);
    }
    //采集多样性的选择
    int population::choice(int a, int b)
    {
        if (P[a].rank < P[b].rank)
            return a;
        else if (P[a].rank == P[b].rank)
        {
            if (P[a].crowding_distance > P[b].crowding_distance)
                return a;
            else
                return b;
        }
        else
            return b;
    }
    void population::print(string fileName)
    {
        FILE* p;
        errno_t err_no = fopen_s(&p, fileName.c_str(), "w+");
        int i, j;
        fprintf(p, "XuYi All Rights Reserved.\nWelcome to OmegaXYZ: www.omegaxyz.com\n");
        fprintf(p, "Problem ZDT1\n");
        fprintf(p, "\n");
        for (i = 0; i < popsize; i++)
        {
            fprintf(p, "The %d generation situation:\n", i);
            for (j = 1; j <= Dimension; j++)
            {
                fprintf(p, "x%d=%e  ", j, P[i].value[j]);
            }
            fprintf(p, "\n");
            fprintf(p, "f1(x)=%0.3f   f2(x)=%0.3f\n", P[i].fvalue[0], P[i].fvalue[1]);
        }
        fclose(p);
        cout << "优化结果保存在[" << fileName << "]，优化结束" << endl;
    }
    void population::exportCsv(string fileName)
    {
        FILE* p;
        errno_t err_no = fopen_s(&p, fileName.c_str(), "w+");
        int i, j;
        fprintf(p, "XuYi All Rights Reserved.\nWelcome to OmegaXYZ: www.omegaxyz.com\n");
        fprintf(p, "Problem ZDT1\n");
        fprintf(p, "F1,F2\n");
        qsort(P, popsize, sizeof(individual), cmp1);
        for (i = 0; i < popsize; i++)
        {
            fprintf(p, "%0.3f,%0.3f\n", P[i].fvalue[0], P[i].fvalue[1]);
        }
        fclose(p);
        cout << "优化结果导出到[" << fileName << "]，优化结束" << endl;
    }
    //主要操作函数
    void population::maincal()
    {
        int gen, i, j;
        gen = generation;
        make_new_pop();
        while (gen--)
        {
            printf("The %d generation\n", generation-gen);
            set_p_q();
            fast_nondominated_sort();
            Pnum = 0;
            i = 0;
            while (Pnum + len[i] <= popsize)
            {
                calu_crowding_distance(i);
                for (j = 0; j < len[i]; j++)
                    P[Pnum++] = F[i][j];
                i++;
                if (i >= len_f)break;
            }
            if (i < len_f)
            {
                calu_crowding_distance(i);
                f_sort(i);
            }
            for (j = 0; j < popsize - Pnum; j++)
                P[Pnum++] = F[i][j];
            make_new_pop();
        }
        print();
        exportCsv();
    }
    population::population()
    {
        srand((unsigned int)(time(0)));
        int i;
        for (i = 0; i < popsize; i++)
        {
            P[i].init();
        }
        for (i = 0; i < popsize; i++)
        {
            P[i].f_count();
        }
        Pnum = popsize;
        Qnum = 0;
        Rnum = 0;
    }
}