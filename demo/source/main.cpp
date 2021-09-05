#include "../include/ZDT1.h"

using namespace std;

int main()
{
	cout << "开始执行:" << endl;
	cout << "\t繁衍代数:" << generation << endl;
	cout << "\t种群大小:" << popsize << endl;
	cout << "\t基因维度:" << xnum << endl;
	cout << "\t目标维度:" << ynum << endl;
	cout << "请稍候..." << endl;

	NSGA2::populationOfZDT1* pop = new NSGA2::populationOfZDT1();
	pop->maincall();
	delete pop;
	char a;
	std::cout << "执行结束,输入任意值退出...";
	scanf_s("%c", &a);
	return 0;
}