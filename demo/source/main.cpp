#include "../include/ZDT1.h"

using namespace std;

int main()
{
	cout << "��ʼִ��:" << endl;
	cout << "\t���ܴ���:" << generation << endl;
	cout << "\t��Ⱥ��С:" << popsize << endl;
	cout << "\t����ά��:" << xnum << endl;
	cout << "\tĿ��ά��:" << ynum << endl;
	cout << "���Ժ�..." << endl;

	NSGA2::populationOfZDT1* pop = new NSGA2::populationOfZDT1();
	pop->maincall();
	delete pop;
	char a;
	std::cout << "ִ�н���,��������ֵ�˳�...";
	scanf_s("%c", &a);
	return 0;
}