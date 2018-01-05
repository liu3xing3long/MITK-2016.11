#ifndef TypeDefs_h__
#define TypeDefs_h__

#include <vector>
#include <map>
#include <queue>
#include <algorithm>
#include "itkImage.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
#define VDIM							3
typedef int 							TPixelType;
typedef itk::Image<TPixelType, VDIM>	ImageType;
typedef itk::Image<TPixelType, 2>		ImageType2D;
typedef ImageType::IndexType			IndexType;
typedef ImageType2D::IndexType			IndexType2D;

typedef unsigned char						BinaryPixelType;
typedef itk::Image<BinaryPixelType, VDIM>	BinaryImageType;
typedef itk::Image<BinaryPixelType, 2>		BinaryImageType2D;

typedef double								DoublePixelType;
typedef itk::Image<DoublePixelType, VDIM>	DoubleImageType;
typedef itk::Image<DoublePixelType, 2>		DoubleImageType2D;

const int		N_EXPERT_COUNT		= 4;
const int		N_LOCAL_WINDOW_R	= 3;
const int		N_NON_NODULE_R      = 5;
const int		N_NON_NODULE_D      = (2 * N_NON_NODULE_R + 1);
const int		N_NODULE_R          = 15;
const int		N_NODULE_D			= (2 * N_NODULE_R + 1);
const int		N_LOCAL_WINDOW_D	= (2 * N_LOCAL_WINDOW_R + 1);

const int		N_NODULE_BOW_R = 32;
const int		N_NODULE_BOW_D = (2 * N_NODULE_BOW_R + 1);

const float		F_CT_WIDTH			= /*650.0f*/1000.0f;
const float		DW_MATH_PI			= 3.14159265359;
const float		F_POSITIVE_MIN		= 0.00000001f;
const float		DW_POSITIVE_MIN		= 1E-16;
const double	DW_SIGMA = 4.5/*0.39894228040141*/;
const double	DW_2PI_SIGMA2 = 1.0 / (2 * DW_MATH_PI * DW_SIGMA * DW_SIGMA);
const double	SQRT_DW_2PI_SIGMA2 = sqrt(1.0 / (2 * DW_MATH_PI * DW_SIGMA * DW_SIGMA));
const float		F_NORMALIZE_SPACE = 1.0000f;
// const int       N_NODULE_TYPES = 4;

//////////////////////////////////////////////////////////////////////////
/// Genetic Structure
//////////////////////////////////////////////////////////////////////////
#define IMAGE_X_DIM 512
#define IMAGE_Y_DIM 512

#ifndef SAFE_DELETE
#define SAFE_DELETE(p)			{ if(p) { delete (p);		(p)=NULL; } }
#endif

#ifndef SAFE_DELETE_ARRAY
#define SAFE_DELETE_ARRAY(p)	{ if(p) { delete [] (p);		(p)=NULL; } }
#endif

#ifndef SAFE_RELEASE
#define SAFE_RELEASE(p)			{ if(p) { (p)->Release();	(p)=NULL; } }
#endif

/************************************************************************/
/* template point 3                                                     */
/************************************************************************/
template<typename T>
struct TPoint3
{
	/**
	*	x, y ����, ��Ӧ�Ҷ�ֵ
	*/
	T x;
	T y;
	T z;
	float fVal;
	float fPtD;
	TPoint3() :x((T)0.0), y((T)0.0), z((T)0.0), fVal(0.00), fPtD(0.0){
	}

	T Length(TPoint3& ptO)
	{
		T tDeltaX = (x - ptO.x);
		T tDeltaY = (y - ptO.y);
		T tDeltaZ = (z - ptO.z);

		return sqrt(tDeltaX*tDeltaX + tDeltaY*tDeltaY + tDeltaZ*tDeltaZ);
	}

	TPoint3(const TPoint3<T>& pt)
	{
		x = pt.x;
		y = pt.y;
		z = pt.z;
		fVal = pt.fVal;
		fPtD = pt.fPtD;
	}
	TPoint3& operator=(const TPoint3<T>& pt)
	{
		x = pt.x;
		y = pt.y;
		z = pt.z;
		fVal = pt.fVal;
		fPtD = pt.fPtD;
		return *this;
	}
	void resetDummy(){
		x = y = z = fVal = fPtD = -1.0;
	}
};

typedef TPoint3<int>   SPoint3I;
typedef TPoint3<float> SPoint3F;
/************************************************************************/
/* template point 2                                                     */
/************************************************************************/
/************************************************************************/
/* !!!!!!!!!PAY ATTENTATION TO IDX UPDATE!!!!!!!!!!                     */
/************************************************************************/

template<typename T>
struct TPoint2
{
	T x, y;
	int iIdx;
	TPoint2(){
		iIdx = -1;
		x = y = -1;
	}
	TPoint2(T tx, T ty){
		x = tx;
		y = ty;
		iIdx = -1;
	}
	/*
	*	copy construction function
	*/
	TPoint2(const TPoint2<T>& pt)
	{
		x = pt.x;
		y = pt.y;
		iIdx = pt.iIdx;
	}
	TPoint2& operator=(const TPoint2<T>& pt)
	{
		x = pt.x;
		y = pt.y;
		iIdx = pt.iIdx;
		return *this;
	}

	bool operator==(const TPoint2<T>& pt)
	{
		if (abs(pt.x - x)<1E-6 && abs(pt.y - y)<1E-6 && (pt.iIdx == iIdx))
			return true;
		else
			return false;
	}

	bool operator!=(const TPoint2<T>& pt)
	{
		if (abs(pt.x - x)>1E-6 || abs(pt.y - y)>1E-6 || (pt.iIdx != iIdx))
			return true;
		else
			return false;
	}

	bool operator<(const TPoint2<T>& pt)
	{
		return (x<pt.x && y<pt.y) ? true : false;
	}
	bool operator>(const TPoint2<T>& pt)
	{
		return (x>pt.x && y>pt.y) ? true : false;
	}
	TPoint2 operator-(const TPoint2<T>& pt)
	{
		return TPoint2(x - pt.x, y - pt.y);
	}
	TPoint2 operator+(const TPoint2<T>& pt)
	{
		return TPoint2(x + pt.x, y + pt.y);
	}
	TPoint2 operator+=(const TPoint2<T>& pt)
	{
		x = x + pt.x;
		y = y + pt.y;
		return *this;
	}
	/*
	*	dummy reset the current vector
	*/
	void resetDummy(){
		x = y = (T)(-1.0);
		iIdx = -1;
	}

	/// this function is depreated and use operator == to replace it
	// 	bool Equal(const TPoint2& pt)
	// 	{
	// 		if (abs(pt.x-x)<1E-6 && abs(pt.y-y)<1E-6 && (pt.iIdx== iIdx))
	// 			return true;
	// 		else 
	// 			return false;
	// 	}

	/*
	*	return length of the current vector
	*/
	float Distance(TPoint2& ptO)
	{
		float tDeltaX = float(x - ptO.x);
		float tDeltaY = float(y - ptO.y);
		return sqrt(tDeltaX*tDeltaX + tDeltaY*tDeltaY);
	}

};

typedef TPoint2<int>   SPoint2I;
typedef TPoint2<float> SPoint2F;



//////////////////////////////////////////////////////////////////////////
struct map2vector
{
	template<typename T>
	const typename T::second_type& operator()(const T& p)
	{
		return p.second;
	}
};

//////////////////////////////////////////////////////////////////////////
/// points related
typedef vector<SPoint3I>	V3IPoints;
typedef vector<SPoint3F>	V3FPoints;
typedef map<int, SPoint3I>	M3IPoints;
typedef map<int, SPoint3F>	M3FPoints;
typedef vector<V3FPoints>	VV3FPoints;
typedef map<int, V3IPoints>	MV3IPoints;
typedef map<int, V3FPoints>	MV3FPoints;
typedef vector<M3IPoints>	VM3IPoints;
typedef map<int, M3IPoints>	MM3IPoints;
typedef vector<SPoint2I>	V2IPoints;
typedef vector<V2IPoints>	VV2IPoints;
typedef vector<SPoint2F>	V2FPoints;
typedef vector<V2FPoints>	VV2FPoints;
typedef map<int, SPoint2I>	M2IPoints;
typedef vector<int>			VecInt;
typedef vector<VecInt>		VVecInt;
typedef vector<double>		VecDW;
typedef vector<VecDW>		VVecDW;

/**  MuM2IPoints enables duplicate insertion*/
typedef multimap<int, V2IPoints> MuM2IPoints;


//////////////////////////////////////////////////////////////////////////
/// Map related containers
typedef std::multimap<int, int>			MMAP;
typedef std::pair<int, int>				MMAPPair;
typedef MMAP::iterator					MMAPIter;
typedef std::pair <MMAPIter, MMAPIter>	MMAPIterPair;

//////////////////////////////////////////////////////////////////////////
/// line template
template <class T>
struct SLineT
{
	vector<T>	vPts;
	T			sBegin;
	T			sEnd;
	float		fMinEnergy;
	float		fLineLength;
	int			nLineIndex;

	SLineT()
	{
		vPts.clear();
		fMinEnergy = 0.f;
		fLineLength = 0.f;
		nLineIndex = -1;
	}


	SLineT(const SLineT& lOther)
	{
		const vector<T>& vOthers = lOther.vPts;
		int nOthers = vOthers.size();
		for (int i = 0; i<nOthers; i++){
			vPts.push_back(vOthers[i]);
		}
		sBegin = lOther.sBegin;
		sEnd = lOther.sEnd;
		fMinEnergy = lOther.fMinEnergy;
		fLineLength = lOther.fLineLength;
		nLineIndex = lOther.nLineIndex;
	}

	void clear()
	{
		vPts.clear();
		fMinEnergy = 0.f;
		fLineLength = 0.f;
		nLineIndex = -1;
		sBegin.resetDummy();
		sEnd.resetDummy();
	}
};



//////////////////////////////////////////////////////////////////////////
/// Application Specific Structure
//////////////////////////////////////////////////////////////////////////
enum ENoduleType
{
	ENT_Nodule, 
	ENT_NonNodule,
	ENT_Count,
};

struct SNoduleBounds
{
	int iXMin;
	int iYMin;
	int iXMax;
	int iYMax;
};
struct SNoduleInfo
{
	V2IPoints	vContour;
	std::string	strUID;
	double		dwZPosition;
	int			nSlicePos;
	int			nExpertId;
	ENoduleType eType;
	bool		bInclusion;
	SNoduleBounds sBounds;

    SNoduleInfo(const SNoduleInfo& otherInfo)
    {
        vContour.clear();
        vContour.insert(vContour.end(), otherInfo.vContour.begin(), otherInfo.vContour.end());
        strUID = otherInfo.strUID;
        dwZPosition = otherInfo.dwZPosition;
        nSlicePos = otherInfo.nSlicePos;
        nExpertId = otherInfo.nExpertId;
        eType = otherInfo.eType;
        bInclusion = otherInfo.bInclusion;
        sBounds.iXMax = otherInfo.sBounds.iXMax;
        sBounds.iXMin = otherInfo.sBounds.iXMin;
        sBounds.iYMax = otherInfo.sBounds.iYMax;
        sBounds.iYMin = otherInfo.sBounds.iYMin;
    }

	SNoduleInfo()
	{
		eType = ENT_Count;
		bInclusion = false;
		sBounds.iXMax = sBounds.iXMin = sBounds.iYMin = sBounds.iYMax = -1;
		nSlicePos = -1;
		strUID = "";
	}

	~SNoduleInfo()
	{
		clear();
	}

	void clear()
	{
		eType = ENT_Count;
		bInclusion = false;
		vContour.clear();
		dwZPosition = -0.1;
	}
};

typedef vector<SNoduleInfo>		VSNoduleInfo;
typedef vector<VSNoduleInfo>	VVSNoduleInfo;
typedef VSNoduleInfo::iterator	VSNoduleInfoIter;
typedef VVSNoduleInfo::iterator	VVSNoduleInfoIter;



#endif // Typedefs_h__
