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
	* @description ���߷��жϵ��Ƿ��ڶ�����ڲ�
	* @param {Object} p ���жϵĵ㣬��ʽ��{ x: X����, y: Y���� }
	* @param {Array} poly ����ζ��㣬�����Ա�ĸ�ʽͬ p
	* @param {int} szPlys ply����Ĵ�С
	* @return {String} �� p �Ͷ���� poly �ļ��ι�ϵ 
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

			// �������ζ����غ�
			if ((isEqual(sx, px) && isEqual(sy, py)) || (isEqual(tx, px) && isEqual(ty, py))) 
				return 0;

			// �ж��߶����˵��Ƿ�����������
			if ((sy < py && ty >= py) || (sy >= py && ty < py)) 
			{
				// �߶��������� Y ������ͬ�ĵ�� X ����
				double x = sx + (py - sy) * (tx - sx) / (double)(ty - sy);
				int iXRound = round(x);
				
				// ���ڶ���εı���
				if (x == px)
					return 0;

				// ���ߴ�������εı߽�
				if (x > px)
					flag = !flag;
				}
			}

		// ���ߴ�������α߽�Ĵ���Ϊ����ʱ���ڶ������
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











