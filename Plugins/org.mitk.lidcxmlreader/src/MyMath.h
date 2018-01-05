/********************************************************************
	created:	2015/10/13
	created:	13:10:2015   10:47
	filename: 	E:\MITK\build\MITK-build\Plugins\org.mitk.lidcxmlreader\MyMath.h
	file path:	E:\MITK\build\MITK-build\Plugins\org.mitk.lidcxmlreader
	file base:	MyMath
	file ext:	h
	author:

	purpose:
	*********************************************************************/

#ifndef MyMath_h__
#define MyMath_h__

#include "Typedefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 


namespace MyMath
{
    int Random(int m, int n)
    {
        int pos, dis;
        if (m == n)
        {
            return m;
        }
        else if (m > n)
        {
            pos = n;
            dis = m - n + 1;
            return rand() % dis + pos;
        }
        else
        {
            pos = m;
            dis = n - m + 1;
            return rand() % dis + pos;
        }
    }

	
	/*
	 *	
	 */
	inline int round(double number){
		number = number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
		return (int)number;
	}

	/*
	 *	
	 */
	template<typename T>
	inline bool isEqual(T num1, T num2)
	{
		if (abs(num1 - num2) < F_POSITIVE_MIN)
			return true;
		else
			return false;
	}



	/**
	* @description 射线法判断点是否在多边形内部
	* @param {Object} p 待判断的点，格式：{ x: X坐标, y: Y坐标 }
	* @param {Array} poly 多边形顶点，数组成员的格式同 p
	* @param {int} szPlys ply数组的大小
	* @return {String} 点 p 和多边形 poly 的几何关系 
	* -1: out; 0: on; 1: in;
	*/
	template<typename T>
	inline int isPointInPolygon(TPoint2<T>& p, TPoint2<T>* poly, int szPlys)
	{
		int px = p.x, py = p.y;
		bool flag = false;

		for (int i = 0, l = szPlys, j = l - 1; i < l; j = i, i++) 
		{
			T sx = poly[i].x;
			T sy = poly[i].y;
			T tx = poly[j].x;
			T ty = poly[j].y;

			// 点与多边形顶点重合
			if ((isEqual(sx, px) && isEqual(sy, py)) || (isEqual(tx, px) && isEqual(ty, py))) 
				return 0;

			// 判断线段两端点是否在射线两侧
			if ((sy < py && ty >= py) || (sy >= py && ty < py)) 
			{
				// 线段上与射线 Y 坐标相同的点的 X 坐标
				double x = sx + (py - sy) * (tx - sx) / (double)(ty - sy);
				int iXRound = round(x);
				
				// 点在多边形的边上
				if (x == px)
					return 0;

				// 射线穿过多边形的边界
				if (x > px)
					flag = !flag;
				}
			}

		// 射线穿过多边形边界的次数为奇数时点在多边形内
		return flag ? 1 : -1;
	}
	template<typename T>
	inline int isPointInPolygon(T x, T y, TPoint2<T>* poly, int szPlys)
	{
		TPoint2<T> p;
		p.x = x;
		p.y = y;
		return isPointInPolygon<T>(p, poly, szPlys);
	}

};









#endif // MyMath_h__











