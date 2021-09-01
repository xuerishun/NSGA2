#include "../include/nsga2.h"
using namespace NSGA2;

int main()
{
	population *pop=new population();
	pop->maincal();
	delete pop;
	return 0;
}