// CLTM.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "estimate.h"
#include "data.h"


int _tmain(int argc, _TCHAR* argv[])
{
	int i=0;
	string vsm_filename="C:\\Users\\hel2\\Desktop\\F-HLda Model\\LDA data set\\LDA-27\\Vocab-6\\VSMFile.txt";
	char dir[100]="C:\\Users\\hel2\\Desktop\\F-HLda Model\\LDA data set\\LDA-27\\Vocab-6\\dTM\\res";
	data d0;
	corpus* cps=d0.read_data(vsm_filename);
	estimate est;
	est.run_vi(cps, dir);
	scanf("%d",i);
	return 0;
}


