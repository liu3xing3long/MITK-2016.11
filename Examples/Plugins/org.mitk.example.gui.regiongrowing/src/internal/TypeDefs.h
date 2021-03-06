#ifndef TypeDefs_h__
#define TypeDefs_h__

#include <itkImage.h>

//////////////////////////////////////////////////////////////////////////
#define VDIM							3

typedef int 							TPixelType;
typedef itk::Image<TPixelType, VDIM>	ImageType;
typedef itk::Image<TPixelType, 2>		ImageType2D;

typedef unsigned char						BinaryPixelType;
typedef itk::Image<BinaryPixelType, VDIM>	BinaryImageType;
typedef itk::Image<BinaryPixelType, 2>		BinaryImageType2D;

typedef double								DoublePixelType;
typedef itk::Image<DoublePixelType, VDIM>	DoubleImageType;
typedef itk::Image<DoublePixelType, 2>		DoubleImageType2D;

typedef float								FloatPixelType;
typedef itk::Image<FloatPixelType, VDIM>	FloatImageType;
typedef itk::Image<FloatPixelType, 2>		FloatImageType2D;

typedef itk::RGBPixel<unsigned char>		RGBPixelType;
typedef itk::Image<RGBPixelType, VDIM>		RGBImageType;
typedef itk::Image<RGBPixelType, 2>		    RGBImageType2D;


const int		N_EXPERT_COUNT = 4;
const int		N_LOCAL_WINDOW_R = 3;
const int		N_NON_NODULE_R = 5;
const int		N_NON_NODULE_D = (2 * N_NON_NODULE_R + 1);
const int		N_NODULE_R = 15;
const int		N_NODULE_D = (2 * N_NODULE_R + 1);
const int		N_LOCAL_WINDOW_D = (2 * N_LOCAL_WINDOW_R + 1);

const float		F_CT_WIDTH = 1000.0f/*650.0f*/;
const float		DW_MATH_PI = 3.14159265359;
const double	F_POSITIVE_MIN = 0.00000001;
const float		DW_POSITIVE_MIN = 1E-16;

// N_EXPORT_DIM MUST be ODD
const int       N_EXPORT_DIM = 65;

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
/* Const Defs                                                           */
/************************************************************************/
#define PROFILE_FUNCTION_BEGIN {\
LARGE_INTEGER litmp; LONGLONG QStart,QEnd;double dfMinus, dfFreq, dfTim;\
QueryPerformanceFrequency(&litmp);\
dfFreq = (double)litmp.QuadPart;\
QueryPerformanceCounter(&litmp);\
QStart = litmp.QuadPart;

#define PROFILE_FUNCTION_H_END(FuncName)\
QueryPerformanceCounter(&litmp);\
QEnd = litmp.QuadPart; dfMinus = (double)(QEnd-QStart);dfTim = dfMinus / dfFreq;\
char strDebug[256];\
sprintf(strDebug, "Profiling %s , Elasped %7.9f ms", #FuncName, dfTim*1000.f);\
MITK_INFO << strDebug;}

#define PROFILE_FUNCTION_D_END(FuncName)\
QueryPerformanceCounter(&litmp);\
QEnd = litmp.QuadPart; dfMinus = (double)(QEnd-QStart);dfTim = dfMinus / dfFreq;\
printf("Profiling %s , Elasped %7.9f ms\n", #FuncName, dfTim*1000.f);}


#endif // TypeDefs_h__
