# spline3
实现三次样条插值
#include <stdio.h>
#include <math.h>


void spline3(int len, int times, int *pxin, float *pyin, int *pxout, float *pyout)
{
	/*
		函数说明： 该函数输出的数据长度（len-1）*times + 1最多为300。
		参量说明：
			len: 输入数据的长度。
			times: 插值的倍数。 输出的长度为（len-1）*times + 1
			*pxin: 输入数据的横坐标
			*pyin: 输入数据的纵坐标
			*pxout: 输出数据的横坐标
			*pyout: 输出数据的纵坐标
	*/

	int n = len - 1;
	 //---------------------------------------构造方程组
  int ii,jj;
	float h[150] = {0.0};   // 之前一直把这个数据设置为int类型，导致h[ii] / (h[ii]+h[ii-1]) 都成了0！！
  float u[150] = {0.0};
  float t[150] = {0.0};
  float d[150] = {0.0};
  float M[150] = {0.0};

	//计算步长  一共给了len个点，所以步长一共有len-1个

  for(ii=0; ii<n; ii++)
  {
			h[ii] = *(pxin+ii+1) - *(pxin+ii);
  }

	// 长度为len-2, n-1
	for (jj=1; jj<=n-1; jj++)
	{
			u[jj] = h[jj-1] / (h[jj-1] + h[jj]);
			t[jj] = 1 - u[jj];
			d[jj] = 6 / ( h[jj-1] + h[jj] ) * ( (pyin[jj+1]-pyin[jj])/h[jj] - (pyin[jj]-pyin[jj-1])/h[jj-1] );
	}

	// ----------------------求解三弯矩方程
	float c[150] = {0};
	float dtemp[150] = {0};

	c[1] = t[1] / 2;
	for (ii=2; ii<=n-2; ii++)
	{
		c[ii] = t[ii] / (2 - c[ii-1]*u[ii]);
	}

	dtemp[1] = d[1] / 2;
	for (ii=2; ii<=n-1; ii++)
	{
		dtemp[ii] = ( d[ii] - dtemp[ii-1]*u[ii] ) / (2 - c[ii-1]*u[ii]);
	}

	// 自然边界
	M[0] = 0;
	M[n] = 0;
	M[n-1] = dtemp[n-1];
	for (ii=n-2; ii>=1 ; ii--)
	{
		M[ii] = dtemp[ii] - c[ii]*M[ii+1];
	}

	// --------------------------------------------利用获取的函数求出插值后的值
	int tempIdx;
	float A, B, C, D;
	int nTimes = times;
	int xIdxEnd = nTimes * (len-1);
	for (ii=0; ii<=xIdxEnd; ii++)
	{
		*(pxout+ii) = ii;
	}

	// tempIdx的范围从0~n-1
	for (ii=0; ii<xIdxEnd; ii++)
	{
		tempIdx = ii/nTimes;
		A = pyin[tempIdx];
		B = ( pyin[tempIdx+1] - pyin[tempIdx] ) / h[tempIdx] - ( 2*M[tempIdx] + M[tempIdx+1] ) / 6 * h[tempIdx];
		C = M[tempIdx] / 2;
		D = (M[tempIdx+1] - M[tempIdx]) / (6*h[tempIdx]);
		*(pyout+ii) = A + B*(ii-pxin[tempIdx]) + C * powf(ii-pxin[tempIdx], 2) + D * powf(ii - pxin[tempIdx], 3) ;
	}
	*(pyout+xIdxEnd) = pyin[n];
}
