#include "MyDrawableNodule.h"
#include "QPainter"
#include "mitkImage.h"

//itk&&vtk headers
#include "mitkImageVtkReadAccessor.h"
#include "mitkImageVtkWriteAccessor.h"
#include "mitkImagePixelReadAccessor.h"
#include "mitkImagePixelWriteAccessor.h"

CMyDrawableNodule::CMyDrawableNodule()
{
	m_vNoduleInfo.resize(N_EXPERT_COUNT);
}


CMyDrawableNodule::~CMyDrawableNodule()
{
	m_vNoduleInfo.clear();
	SAFE_DELETE(m_originalImage);
}

//////////////////////////////////////////////////////////////////////////
void CMyDrawableNodule::paintEvent(QPaintEvent *event)
{
	QLabel::paintEvent(event);
	QPainter painter(this);

    const int I_PAINTER_SIZE = 1;

	for (int i = 0; i < N_EXPERT_COUNT; i++)
	{
		VSNoduleInfo& vInfoThisExpert = m_vNoduleInfo[i];
		if (vInfoThisExpert.empty())
			continue;

		int szInfoThisExpert = vInfoThisExpert.size();
		for (int j = 0; j < szInfoThisExpert; j++)
		{
			SNoduleInfo& thisNoduleInfor = vInfoThisExpert[j];
			SNoduleBounds& sBounds = thisNoduleInfor.sBounds;
			int iWidth = sBounds.iXMax - sBounds.iXMin;
			int iHeight = sBounds.iYMax - sBounds.iYMin;
			if (iWidth == 0)
				iWidth = 1;
			if (iHeight == 0)
				iHeight = 1;

			switch (i)
			{
			case 0:
                painter.setPen(QPen(Qt::red, I_PAINTER_SIZE));
				break;
			case 1:
                painter.setPen(QPen(Qt::green, I_PAINTER_SIZE));
				break;
			case 2:
                painter.setPen(QPen(Qt::blue, I_PAINTER_SIZE));
				break;
			case 3:
                painter.setPen(QPen(Qt::yellow, I_PAINTER_SIZE));
				break;
			}

			switch (thisNoduleInfor.eType)
			{
            case ENT_Nodule:
            {
                painter.drawRect(/*QRect*/sBounds.iXMin, sBounds.iYMin, iWidth, iHeight);
                //for (int iPT = 0; iPT < thisNoduleInfor.vContour.size(); iPT++)
                //{
                 //   painter.drawPoint(QPoint(thisNoduleInfor.vContour[iPT].x, thisNoduleInfor.vContour[iPT].y));
                //}
            }
            break;

            case ENT_NonNodule:
            {
                painter.drawEllipse(QPointF(sBounds.iXMin, sBounds.iYMin), N_NODULE_R, N_NODULE_R);
            }
            break;
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////
void CMyDrawableNodule::addNoduleInfo(int idxExpert, SNoduleInfo& sNoduleInfo)
{
	m_vNoduleInfo[idxExpert].push_back(sNoduleInfo);
}

//////////////////////////////////////////////////////////////////////////
void CMyDrawableNodule::setBaseImageData(mitk::Image* mitkImage, int nPosSlice)
{
	if (NULL != mitkImage)
		m_originalImage = NULL;

	int width = mitkImage->GetDimension(0);
	int height = mitkImage->GetDimension(1);
	m_originalImage = new QImage(width, height, QImage::Format_RGB32);

	mitk::Image::ImageDataItemPointer pSlice = mitkImage->GetSliceData(nPosSlice);
	mitk::ImagePixelReadAccessor<TPixelType, 2> readAccess(mitkImage, pSlice);
	itk::Index<2> read_idx;

	QRgb* rgbPtr = reinterpret_cast<QRgb*>(m_originalImage->bits());
	for (int row = 0; row < height; ++row)
	{
		for (int col = 0; col < width; ++col)
		{
			// Swap rgb
			read_idx[0] = col;
			read_idx[1] = row;

			int iOriginalCTValue = readAccess.GetPixelByIndex(read_idx);
			int iPixelValue;
			if (iOriginalCTValue < -F_CT_WIDTH)
				iPixelValue = 0;
			else if (iOriginalCTValue > F_CT_WIDTH)
				iPixelValue = 255;
			else
				iPixelValue = round((double)(iOriginalCTValue + F_CT_WIDTH) / (2 * F_CT_WIDTH) * 255);

			*(rgbPtr++) = QColor(iPixelValue, iPixelValue, iPixelValue).rgb();
			//	colorsPtr += 1;
		}
		//rgbPtr -= width * 2;
	}

	*m_originalImage = m_originalImage->scaled(IMAGE_X_DIM, IMAGE_Y_DIM, Qt::KeepAspectRatio);
	this->setPixmap(QPixmap::fromImage(*m_originalImage));

}

//////////////////////////////////////////////////////////////////////////
void CMyDrawableNodule::setBaseImageFromVtkData(vtkImageData* imageData)
{
	if (!imageData)
		m_originalImage = NULL;

	int width = imageData->GetDimensions()[0];
	int height = imageData->GetDimensions()[1];
	m_originalImage = new QImage(width, height, QImage::Format_RGB32);

	QRgb* rgbPtr = reinterpret_cast<QRgb*>(m_originalImage->bits());
	int* colorsPtr = reinterpret_cast<int*>(imageData->GetScalarPointer());
	// mirror vertically
	for (int row = 0; row < height; ++row)
	{
		for (int col = 0; col < width; ++col)
		{
			// Swap rgb
			int iOriginalCTValue = colorsPtr[0];
			int iPixelValue;
			if (iOriginalCTValue < -F_CT_WIDTH)
				iPixelValue = 0;
			else if (iOriginalCTValue > F_CT_WIDTH)
				iPixelValue = 255;
			else
				iPixelValue = round((iOriginalCTValue + F_CT_WIDTH) / (2 * F_CT_WIDTH) * 255);

			*(rgbPtr++) = QColor(iPixelValue, iPixelValue, iPixelValue).rgb();
			colorsPtr += 1;
		}
		//rgbPtr -= width * 2;
	}
	this->setPixmap(QPixmap::fromImage(*m_originalImage));
}