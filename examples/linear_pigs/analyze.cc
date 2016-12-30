#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>  
#include <math.h>
using namespace std;

int main()
{
	int nSizeRow = 10000;
	int nSizeCol = 10;
    double myArray[nSizeRow][nSizeCol];

	ifstream file("pigs.eng");
    if(file.is_open())
    {
        for(int i = 0; i < nSizeRow; ++i)
        {
            for(int j = 0; j < nSizeCol; ++j)
            {
        		file >> myArray[i][j];
            }
        }
    }
	
	double sum1Mean   = 0.0;
	double sum2Mean   = 0.0;
	double sum1MeanSq = 0.0;
	double sum2MeanSq = 0.0;

	for (int i = 0; i<nSizeRow; ++i)
	{
		sum1Mean     += myArray[i][2];
		sum2Mean     += myArray[i][3];
		sum1MeanSq   += myArray[i][2]*myArray[i][2];
		sum2MeanSq   += myArray[i][3]*myArray[i][3];
	}

	double meanKin    = sum1Mean/double(nSizeRow);
	double meanPot    = sum2Mean/double(nSizeRow);
    double meanSqKin  = sum1MeanSq/double(nSizeRow);
    double meanSqPot  = sum2MeanSq/double(nSizeRow);
	double errorKin   = sqrt((meanSqKin - meanKin*meanKin)/nSizeRow);
	double errorPot   = sqrt((meanSqPot - meanPot*meanPot)/nSizeRow);
	cout << errorKin <<"   "<<errorPot<<endl;	
	return (0);
}

