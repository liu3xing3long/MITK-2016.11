#ifndef MyDrawableNodule_h__
#define MyDrawableNodule_h__



#include <QLabel>
#include <vector>
#include <mitkImage.h>
#include "Typedefs.h"


class vtkImageData;

class CMyDrawableNodule :
	public QLabel
{
public:
	CMyDrawableNodule();
	~CMyDrawableNodule();

private:
	void paintEvent(QPaintEvent *event);

public:
	void addNoduleInfo(int idxExpert, SNoduleInfo& sNoduleInfo);
	void setBaseImageFromVtkData(vtkImageData* pVTKImg);
	void setBaseImageData(mitk::Image* mitkImage, int nPosSlice);

protected:
	QImage vtkImageDataToQImage(vtkImageData* imageData);


private:
	/* nodule infor corresponding to this image*/
	VVSNoduleInfo m_vNoduleInfo;
	QImage* m_originalImage;
};

#endif // MyDrawableNodule_h__
