/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/


// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>
#include "org_mitk_example_gui_regiongrowing_Activator.h"

// Qmitk
#include "QmitkRegionGrowingView.h"
#include "QmitkPointListWidget.h"
#include "QmitkRenderWindow.h"
#include "QmitkDataNodeSelectionProvider.h"
#include "mitkNodePredicateDataType.h"
#include "mitkNodePredicateDimension.h"
#include "mitkNodePredicateAnd.h"
#include "mitkDataNodeObject.h"

// MITK
#include "mitkImageAccessByItk.h"
#include "mitkProperties.h"
#include "mitkColorProperty.h"
#include "mitkIOUtil.h"
#include "mitkImagePixelReadAccessor.h"
#include "mitkImagePixelWriteAccessor.h"
#include "mitkImageGenerator.h"
#include "mitkExtractSliceFilter.h"

// Includes for image casting between ITK and MITK
#include "mitkImageCast.h"
#include "mitkITKImageImport.h"

// ITK includes (general)
#include <itkVectorImage.h>
#include <itkImageFileWriter.h>

// ITK
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkIsolatedConnectedImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageMaskSpatialObject.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkLabelToRGBImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkExtractImageFilter.h"
#include "itkJoinSeriesImageFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"
#include "itkParabolicErodeDilateImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkLabelImageToStatisticsLabelMapFilter.h"
#include "itkRGBToLuminanceImageFilter.h"

#include "itkDecisionRule.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkMinimumDecisionRule.h"
#include "itkEuclideanDistanceMetric.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkSampleClassifierFilter.h"
#include "itkNormalVariateGenerator.h"
#include "itkBilateralImageFilter.h"

// hessian related
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"

#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkVotingBinaryHoleFillFloodingImageFilter.h"
// Smoothing
#include <itkMedianImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>

// Threshold
#include <itkBinaryThresholdImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include "CudaBinaryThresholdImageFilter.h"
#include "CudaGrayscaleDilateImageFilter.h"
#include "CudaMedianImageFilter.h"
#include "CudaVesselnessImageFilter.h"
#include "CudaMultiplyByConstantImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "CudaDivideByConstantImageFilter.h"
#include "CudaGrayscaleErodeImageFilter.h"

// Boolean operations
#include <itkOrImageFilter.h>
#include <itkAndImageFilter.h>
#include <itkXorImageFilter.h>

// Qt
#include <QMessageBox>
#include <QFile>
#include <QFileDialog>
#include <QTreeView>  
#include <QList>  
#include <QStandardItem>  
#include <QColorDialog>

// vtk related
#include "vtkImageWriter.h"
#include "vtkJPEGWriter.h"
#include "vtkPNGWriter.h"
#include "vtkRenderLargeImage.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkTestUtilities.h"


//////////////////////////////////////////////////////////////////////////
const std::string QmitkRegionGrowingView::VIEW_ID = "org.mitk.views.example.regiongrowing";
QString G_ANN_INDEX_FILE_NAME = "annIndexes.txt";
QString G_ANN_REGION_FILE_NAME = "annRegionParams.txt";
QString G_ANN_FORECAST_FILE_NAME = "annForecastData.txt";
QString G_ANN_INDEX_VALUES_PREFIX = "annIndexValues_";

QString STR_ANN_RESULTS_TYPE[QmitkRegionGrowingView::EANNResultsType::EANNRT_COUNT] = {
    "ISO",
    "WALL",
    "VESSEL",
    "GGO",
};

const unsigned int RESULT_PRECISION = 12;
const unsigned int PIXEL_PRECISION = 4;

// typedef itk::BinaryThresholdImageFilter<ImageType, ImageType>			ThresholdFilterType;
typedef itk::AndImageFilter< ImageType, ImageType >						AndImageFilterType;
typedef itk::MultiplyImageFilter< ImageType, ImageType, ImageType >		MultiplyFilterType;



//////////////////////////////////////////////////////////////////////////
template<typename T>
inline void limit_number(T& intput, T minimum, T maximum)
{
    if (intput < minimum)
        intput = minimum;
    if (intput > maximum)
        intput = maximum;
}

//////////////////////////////////////////////////////////////////////////
QmitkRegionGrowingView::QmitkRegionGrowingView()
    : m_PointListWidget(NULL), m_WorkingDir(""), m_BackgroundColor(QColor(0, 0, 0))

{
}

void QmitkRegionGrowingView::SetFocus()
{
    m_Controls.buttonPerformImageProcessing->setFocus();
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::CreateQtPartControl(QWidget *parent)
{
    //-------------------------------------------------------------------------
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.buttonPerformImageProcessing, SIGNAL(clicked()), this, SLOT(DoRegionGrowing()));
    connect(m_Controls.btnHessianAnalysis, SIGNAL(clicked()), this, SLOT(DoHessian()));
    connect(m_Controls.btnDoThreshold, SIGNAL(clicked()), this, SLOT(DoThreshold()));
    // 	connect(m_Controls.btNormalize, SIGNAL(clicked()), this, SLOT(DoNormalization()));
    // 	connect(m_Controls.btnMethodTest, SIGNAL(clicked()), this, SLOT(DoMethodTest()));
    // 	connect(m_Controls.btMethodTest2, SIGNAL(clicked()), this, SLOT(DoMethodTest2()));
    connect(m_Controls.btSelectRegion, SIGNAL(clicked()), this, SLOT(DoSelectRegion()));
    connect(m_Controls.btFillHole, SIGNAL(clicked()), this, SLOT(DoFillHole()));
    connect(m_Controls.btExtractRealData, SIGNAL(clicked()), this, SLOT(DoExtractRealData()));
    connect(m_Controls.btErosion, SIGNAL(clicked()), this, SLOT(DoErosion()));
    connect(m_Controls.btLabelImage, SIGNAL(clicked()), this, SLOT(DoColorLabel()));
    connect(m_Controls.btLabelImage_2, SIGNAL(clicked()), this, SLOT(DoColorLabel_2()));
    connect(m_Controls.btOutputForecast, SIGNAL(clicked()), this, SLOT(DoOutputForecast()));
    // 	connect(m_Controls.btImageCrop, SIGNAL(clicked()), this, SLOT(DoCrop()));
    // 	connect(m_Controls.btLoadIndexes, SIGNAL(clicked()), this, SLOT(DoLoadIndexes()));
    // 	connect(m_Controls.btLoadIndexValues, SIGNAL(clicked()), this, SLOT(DoLoadIndexValues()));
    connect(m_Controls.btShowANNResults, SIGNAL(clicked()), this, SLOT(DoShowAnnResults()));
    // 	connect(m_Controls.btRollingBall, SIGNAL(clicked()), this, SLOT(DoRollingBall()));
    connect(m_Controls.btSelWorkDir, SIGNAL(clicked()), this, SLOT(DoSelectWorkingDir()));
    connect(m_Controls.btDoThinning, SIGNAL(clicked()), this, SLOT(DoThinning()));
    connect(m_Controls.btOstyThreshold, SIGNAL(clicked()), this, SLOT(DoOSTUThreshold()));
    //     connect(m_Controls.btKMeans, SIGNAL(clicked()), this, SLOT(DoKMeans()));
    connect(m_Controls.btSelectLabel, SIGNAL(clicked()), this, SLOT(DoSelectLabel()));
    connect(m_Controls.btDoDilation, SIGNAL(clicked()), this, SLOT(DoDilation()));
    connect(m_Controls.btDoMedian, SIGNAL(clicked()), this, SLOT(DoMedian()));
    connect(m_Controls.btFillBorder, SIGNAL(clicked()), this, SLOT(DoFillBorder()));
    connect(m_Controls.btTransMask, SIGNAL(clicked()), this, SLOT(DoTransferToMask()));
    connect(m_Controls.rbBGPixelRed, SIGNAL(clicked()), this, SLOT(DoChangeBGColor()));
    connect(m_Controls.rbBGPixelBlack, SIGNAL(clicked()), this, SLOT(DoChangeBGColor()));
    connect(m_Controls.btScreenShot, SIGNAL(clicked()), this, SLOT(DoTakeScreenShots()));
    connect(m_Controls.btOutputSelectNodules, SIGNAL(clicked()), this, SLOT(DoOutputScoringResults()));
    connect(m_Controls.btnRedFaces, SIGNAL(clicked()), this, SLOT(DoReduceFaces()));
    connect(m_Controls.btRefineRG, SIGNAL(clicked()), this, SLOT(DoRefineRegionGrowing()));
    connect(m_Controls.btOutputBows, SIGNAL(clicked()), this, SLOT(DoOutputBOWSamples()));
    //-------------------------------------------------------------------------
    // test methods
    connect(m_Controls.btnMethodTest1, SIGNAL(clicked()), this, SLOT(DoMethodTest1()));
    connect(m_Controls.btnMethodTest2, SIGNAL(clicked()), this, SLOT(DoMethodTest2()));
    connect(m_Controls.btnMethodTest3, SIGNAL(clicked()), this, SLOT(DoMethodTest3()));
    connect(m_Controls.btnMethodTest4, SIGNAL(clicked()), this, SLOT(DoMethodTest4()));


    //-------------------------------------------------------------------------
    mitk::NodePredicateDimension::Pointer dimensionPredicate = mitk::NodePredicateDimension::New(3);
    mitk::NodePredicateDataType::Pointer imagePredicate = mitk::NodePredicateDataType::New("Image");
    m_Controls.cbOriginalImage->SetDataStorage(GetDataStorage());
    m_Controls.cbOriginalImage->SetPredicate(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate));
    m_Controls.cbMaskImage->SetDataStorage(GetDataStorage());
    m_Controls.cbMaskImage->SetPredicate(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate));
    m_Controls.cbOriginalImage_2->SetDataStorage(GetDataStorage());
    m_Controls.cbOriginalImage_2->SetPredicate(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate));
    m_Controls.cbMaskImage_2->SetDataStorage(GetDataStorage());
    m_Controls.cbMaskImage_2->SetPredicate(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate));
    m_Controls.cbStatisticImage->SetDataStorage(GetDataStorage());
    m_Controls.cbStatisticImage->SetPredicate(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate));
    m_Controls.cbStatisticFeatureImage->SetDataStorage(GetDataStorage());
    m_Controls.cbStatisticFeatureImage->SetPredicate(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate));

    //-------------------------------------------------------------------------
    /// create a QmitkPointListWidget and add it to the widget created from .ui file
    m_PointListWidget = new QmitkPointListWidget();
    m_Controls.verticalLayout_10->addWidget(m_PointListWidget, 1);
    // retrieve a possibly existing IRenderWindowPart
    if (mitk::IRenderWindowPart* renderWindowPart = GetRenderWindowPart())
    {
        // let the point set widget know about the render window part (crosshair updates)
        RenderWindowPartActivated(renderWindowPart);
    }
    // create a new DataNode containing a PointSet with some interaction
    m_PointSet = mitk::PointSet::New();
    mitk::DataNode::Pointer pointSetNode = mitk::DataNode::New();
    pointSetNode->SetData(m_PointSet);
    pointSetNode->SetName("seed points for region growing");
    pointSetNode->SetProperty("helper object", mitk::BoolProperty::New(true));
    pointSetNode->SetProperty("layer", mitk::IntProperty::New(1024));
    // add the pointset to the data storage (for rendering and access by other modules)
    GetDataStorage()->Add(pointSetNode);
    // tell the GUI widget about the point set
    m_PointListWidget->SetPointSetNode(pointSetNode);

    //-------------------------------------------------------------------------
    m_Controls.treeWidget->setColumnCount(2);
    QStringList header_name;
    header_name << "Index";
    header_name << "Origin";
    header_name << "Size";
    header_name << "Value";
    m_Controls.treeWidget->setHeaderLabels(header_name);

    //-------------------------------------------------------------------------
    QString styleSheet = "background-color:rgb(0,0,0)";
    m_Controls.btBGColor->setStyleSheet(styleSheet);
}


//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/,
    const QList<mitk::DataNode::Pointer>& nodes)
{
    // iterate all selected objects, adjust warning visibility
    foreach(mitk::DataNode::Pointer node, nodes)
    {
        mitk::Image *pSel = dynamic_cast<mitk::Image *>(node->GetData());
        if (NULL != pSel && node.IsNotNull())
        {
            m_Controls.buttonPerformImageProcessing->setEnabled(true);
            return;
        }
    }
    m_Controls.buttonPerformImageProcessing->setEnabled(false);
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::RenderWindowPartActivated(mitk::IRenderWindowPart* renderWindowPart)
{
    // let the point set widget know about the slice navigation controllers
    // in the active render window part (crosshair updates)
    foreach(QmitkRenderWindow* renderWindow, renderWindowPart->GetQmitkRenderWindows().values())
    {
        m_PointListWidget->AddSliceNavigationController(renderWindow->GetSliceNavigationController());
    }

    if (this->m_IRenderWindowPart != renderWindowPart)
    {
        this->m_IRenderWindowPart = renderWindowPart;
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::RenderWindowPartDeactivated(mitk::IRenderWindowPart* renderWindowPart)
{
    foreach(QmitkRenderWindow* renderWindow, renderWindowPart->GetQmitkRenderWindows().values())
    {
        m_PointListWidget->RemoveSliceNavigationController(renderWindow->GetSliceNavigationController());
    }
}

//////////////////////////////////////////////////////////////////////////
mitk::DataNode* QmitkRegionGrowingView::getSelectedNode()
{
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (!nodes.empty())
    {
        mitk::DataNode* node = nodes.front();
        if (node)
        {
            return node;
        }
    }
    return NULL;
}

//////////////////////////////////////////////////////////////////////////
mitk::Image* QmitkRegionGrowingView::getSelectedImage()
{
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (!nodes.empty())
    {
        mitk::DataNode* node = nodes.front();
        if (node)
        {
            // a node itself is not very useful, we need its data item (the image)
            mitk::BaseData*	data = node->GetData();

            if (data)
            {
                // test if this data item is an image or not (could also be a surface or something totally different)
                mitk::Image* image = dynamic_cast<mitk::Image*>(data);
                if (image)
                {
                    return image;
                }
            }
        }
    }
    MITK_WARN << "No Selected Image!!";
    return NULL;
}

//////////////////////////////////////////////////////////////////////////
QString QmitkRegionGrowingView::getLastOpenPath()
{
    berry::IPreferencesService* prefService = mitk::PluginActivator::GetInstance()->GetPreferencesService();

    berry::IPreferences::Pointer  prefs = berry::IPreferences::Pointer(0);
    if (prefService)
        prefs = prefService->GetSystemPreferences()->Node("/General");

    if (prefs.IsNotNull())
    {
        QString qLastPath = prefs->Get("LastFileOpenPath", "");
        int		nStartIndex = qLastPath.lastIndexOf("/");
        qLastPath = qLastPath.mid(0, nStartIndex);
        return qLastPath;
    }
    return QString();
}

////////////////////////////////////////////////////////////////////////////
//QString QmitkRegionGrowingView::getCurrentDataPath()
//{
//    mitk::Image::Pointer mitkImage = getSelectedImage();
//    if (NULL != mitkImage)
//    {
//        mitk::PropertyList* pProList = mitkImage->GetPropertyList();
//        std::string strPath;
//        pProList->GetStringProperty("path", strPath);
//        QString qSTRPath(strPath.c_str());
//        int		nStartIndex = qSTRPath.lastIndexOf("/");
//        QString qSTRImageName = qSTRPath.mid(0, nStartIndex);
//        return qSTRImageName;
//    }
//
//    return "";
//}

////////////////////////////////////////////////////////////////////////////
//void QmitkRegionGrowingView::DoMethodTest2()
//{
//	mitk::Image::Pointer mitkImage = getSelectedImage();
//
//	if (NULL != mitkImage)
//	{
//
//	}
//}

//////////////////////////////////////////////////////////////////////////

//
////////////////////////////////////////////////////////////////////////////
//void QmitkRegionGrowingView::DoMethodTest()
//{
//	QString strFormat;
//	mitk::Image* pSel = getSelectedImage();
//	if (pSel == NULL)
//		return;
//	mitk::Image::Pointer pNew = pSel->Clone();
//	ImageType::Pointer itkImage;
//	mitk::CastToItkImage(pNew, itkImage);
//
//	//-------------------------------------------------------------------------
//	/// cuda median filter test
//#if 0
//	// median filter to smooth the images
//	typedef itk::MedianImageFilter< ImageType, ImageType> MedianFilterType;
//	MedianFilterType::Pointer medianFilter = MedianFilterType::New();
//	MedianFilterType::InputSizeType size;
//	int param1 = 3;
//	size.Fill(param1);
//	medianFilter->SetRadius(size);
//	medianFilter->SetInput(itkImage);
//	PROFILE_FUNCTION_BEGIN;
//	medianFilter->Update();
//	PROFILE_FUNCTION_H_END(CPU_MEDIAN);
//	mitk::Image::Pointer mitkResult = mitk::Image::New();
//	mitk::CastToMitkImage(medianFilter->GetOutput(), mitkResult);
//	strFormat.sprintf("CPU_Median_%d", param1);
//	mitk::DataNode::Pointer newNode = mitk::DataNode::New();
//	newNode->SetData(mitkResult);
//	newNode->SetProperty("name", mitk::StringProperty::New(strFormat.toStdString()));
//	this->GetDataStorage()->Add(newNode);
//	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
//
//	typedef itk::CudaMedianImageFilter< ImageType, ImageType> CudaMedianFilterType;
//	CudaMedianFilterType::Pointer medianFilter2 = CudaMedianFilterType::New();
//	medianFilter2->SetRadius(size);
//	medianFilter2->SetInput(itkImage);
//	PROFILE_FUNCTION_BEGIN;
//	//medianFilter2->UpdateLargestPossibleRegion();
//	medianFilter2->Update();
//	PROFILE_FUNCTION_H_END(GPU_MEDIAN);
//	mitk::Image::Pointer mitkResult2 = mitk::Image::New();
//	mitk::CastToMitkImage(medianFilter2->GetOutput(), mitkResult2);
//	strFormat.sprintf("GPU_Median_%d", param1);
//	mitk::DataNode::Pointer newNode2 = mitk::DataNode::New();
//	newNode2->SetData(mitkResult2);
//	newNode2->SetProperty("name", mitk::StringProperty::New(strFormat.toStdString()));
//	this->GetDataStorage()->Add(newNode2);
//	mitk::RenderingManager::GetInstance()->RequestUpdateAll();
//#endif
//
//#if 0
//	//-------------------------------------------------------------------------
//	// cuda hessian test
//	typedef itk::Image<double> VesselnessOutputType;
//	typedef itk::CudaVesselnessImageFilter< ImageType, VesselnessOutputType> CudaVesselnessFilterType;
//	CudaVesselnessFilterType::InputSizeType size;
//	int param1 = 3;
//	size.Fill(param1);
//	CudaVesselnessFilterType::Pointer vesselnessFilter = CudaVesselnessFilterType::New();
//	vesselnessFilter->SetRadius(size);
//	vesselnessFilter->SetInput(itkImage);
//	PROFILE_FUNCTION_BEGIN;
//	vesselnessFilter->Update();
//	PROFILE_FUNCTION_H_END(GPU_HESSIAN);
//#endif
//
//}


////////////////////////////////////////////////////////////////////////////
//void QmitkRegionGrowingView::DoNormalization()
//{
//	mitk::Image*		mitkImage = getSelectedImage();
//	mitk::DataNode*		mitkNode = getSelectedNode();
//	ImageType::Pointer itkImage;
//	mitk::CastToItkImage(mitkImage, itkImage);
//
//	typedef itk::ResampleImageFilter< ImageType, ImageType >			ResampleImageFilterType;
//	typedef itk::LinearInterpolateImageFunction< ImageType, double >	LinearInterpolatorType;
//
//	if (NULL != mitkImage)
//	{
//		ResampleImageFilterType::Pointer	resampler = ResampleImageFilterType::New();
//		LinearInterpolatorType::Pointer		interpolator = LinearInterpolatorType::New();
//		resampler->SetInterpolator(interpolator);
//		resampler->SetInput(itkImage);
//		resampler->SetOutputOrigin(itkImage->GetOrigin());
//		resampler->SetDefaultPixelValue(-1024); // Hounsfield Units for Air
//
//		ImageType::SizeType			input_size = itkImage->GetLargestPossibleRegion().GetSize();
//		ImageType::SpacingType		input_spacing = itkImage->GetSpacing();
//		ImageType::SpacingValueType input_min_spacing = std::min(input_spacing[0], input_spacing[1]);
//		input_min_spacing = std::min(input_min_spacing, input_spacing[2]);
//
//		ImageType::SizeType		output_size;
//		ImageType::SpacingType	output_spacing;
//
//		output_size[0] = input_size[0] * (input_spacing[0] / input_min_spacing);
//		output_size[1] = input_size[1] * (input_spacing[1] / input_min_spacing);
//		output_size[2] = input_size[2] * (input_spacing[2] / input_min_spacing);
//		output_spacing[0] = input_min_spacing;
//		output_spacing[1] = input_min_spacing;
//		output_spacing[2] = input_min_spacing;
//
//		QString debugSTR;
//		debugSTR.sprintf("Resample...Original[Spacing(%f, %f, %f), Size(%f, %f, %f)], Output[Spacing(%f, %f, %f), Size(%f, %f, %f)]",
//			input_spacing[0], input_spacing[1], input_spacing[2], input_size[0], input_size[1], input_size[2],
//			output_spacing[0], output_spacing[1], output_spacing[2], output_size[0], output_size[1], output_size[2]);
//		MITK_INFO << debugSTR;
//
//		resampler->SetSize(output_size);
//		resampler->SetOutputSpacing(output_spacing);
//		resampler->SetOutputDirection(itkImage->GetDirection());
//		resampler->UpdateLargestPossibleRegion();
//
//		// show results
//		mitk::Image::Pointer pResult = mitk::Image::New();
//		mitk::CastToMitkImage(resampler->GetOutput(), pResult);
//		QString qNodeName;
//		qNodeName.sprintf("%s_resample_%f", mitkNode->GetName().c_str(), input_min_spacing);
//		//		mStrOperation += qNodeName;
//		ShowProcessedDataNode(pResult, qNodeName.toStdString());
//	}
//}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoThreshold()
{
    QString strParam1 = m_Controls.leLowThreshold->text();
    QString strParam2 = m_Controls.leHighThreshold->text();
    float	param1 = strParam1.toInt();
    float	param2 = strParam2.toInt();
    mitk::Image*	pSelected = getSelectedImage();
    mitk::DataNode* pSelNode = getSelectedNode();

    if (pSelNode == NULL || pSelected == NULL)
    {
        MITK_INFO << "NO Image Selected";
        return;
    }

    mitk::Image::Pointer mitkNewImage = pSelected->Clone();
    ImageType::Pointer itkImage  /*= ImageType::New()*/;
    mitk::CastToItkImage(mitkNewImage, itkImage);

    typedef itk::CudaBinaryThresholdImageFilter<ImageType, ImageType>		ThresholdFilterType;
    // binary threshold
    ThresholdFilterType::Pointer thFilter = ThresholdFilterType::New();
    MITK_INFO << "Thresholding with" << param1 << " and " << param2;
    thFilter->SetInput(itkImage);
    thFilter->SetLowerThreshold(param1);
    thFilter->SetUpperThreshold(param2);
    thFilter->SetInsideValue(1);
    thFilter->SetOutsideValue(0);
    thFilter->Update();

    //  cast image to binary
    typedef itk::CastImageFilter< ImageType, BinaryImageType > ImageCastFilter;
    ImageCastFilter::Pointer castFilter = ImageCastFilter::New();
    castFilter->SetInput(thFilter->GetOutput());
    castFilter->Update();

    // show results
    mitk::Image::Pointer pResult = mitk::Image::New();
    mitk::CastToMitkImage(castFilter->GetOutput(), pResult);
    QString qNodeName;
    qNodeName.sprintf("%s_thre_%f_%f", pSelNode->GetName().c_str(), param1, param2);
    ShowProcessedDataNode(pResult, qNodeName.toStdString(), true, pSelNode);
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoHessian()
{
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        QString			strSigStart = m_Controls.leHessian1->text();
        QString			strSigEnd = m_Controls.leHessian2->text();
        QString			strSigInterval = m_Controls.leHessian3->text();
        float			fSigStart = strSigStart.toFloat();
        float			fSigEnd = strSigEnd.toFloat();
        int				nNumOfSigmas = strSigInterval.toInt();

        ImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        float falpha = 0.5f;
        float fbeta = 0.5f;
        float fc = 500.0f;

        typedef itk::SymmetricSecondRankTensor< double, VDIM >	HessianPixelType;
        typedef itk::Image< HessianPixelType, VDIM>				HessianImageType;

        ImageType::SizeType imageSize = mCurrentRegion.GetSize();
        if (imageSize[0] == 0 || imageSize[1] == 0 || imageSize[2] == 0)
        {
            MITK_WARN << "No ROI Set, Using Whole Image..." << mCurrentRegion;
            mCurrentRegion = itkImage->GetLargestPossibleRegion();
        }

#if 0
        typedef itk::HessianToObjectnessMeasureImageFilter< HessianImageType, /*ImageType*/DoubleImageType> ObjectnessFilterType;
        ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
        if (m_Controls.cbHessianRevert->isChecked()){
            objectnessFilter->SetBrightObject(true);
        }
        else{
            objectnessFilter->SetBrightObject(false);
        }
        objectnessFilter->SetScaleObjectnessMeasure(false);
        objectnessFilter->SetAlpha(falpha);
        objectnessFilter->SetBeta(fbeta);
        objectnessFilter->SetGamma(fc);

        typedef itk::MultiScaleHessianBasedMeasureImageFilter< ImageType, HessianImageType, /*ImageType */DoubleImageType>
            MultiScaleEnhancementFilterType;
        MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter =
            MultiScaleEnhancementFilterType::New();
        multiScaleEnhancementFilter->SetInput(itkImage);
        multiScaleEnhancementFilter->SetHessianToMeasureFilter(objectnessFilter);
        multiScaleEnhancementFilter->SetSigmaStepMethodToEquispaced();
        multiScaleEnhancementFilter->SetSigmaMinimum(fSigStart);
        multiScaleEnhancementFilter->SetSigmaMaximum(fSigEnd);
        multiScaleEnhancementFilter->SetNumberOfSigmaSteps(nNumOfSigmas);
        MITK_INFO << "Operating on " << mCurrentRegion;
        multiScaleEnhancementFilter->Update();
#else
        // Declare the type of enhancement filter - use ITK’s 3D vesselness (Sato)
        typedef itk::Hessian3DToVesselnessMeasureImageFilter<double> VesselnessFilterType;
        // Declare the type of multiscale enhancement filter
        typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType, HessianImageType, DoubleImageType> MultiScaleEnhancementFilterType;

        VesselnessFilterType::Pointer vesselnessFilter = VesselnessFilterType::New();
        vesselnessFilter->SetAlpha1(0.25);
        vesselnessFilter->SetAlpha2(0.5);

        // Instantiate the multiscale filter and set the input image
        MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
        multiScaleEnhancementFilter->SetInput(itkImage);
        multiScaleEnhancementFilter->SetHessianToMeasureFilter(vesselnessFilter);
        multiScaleEnhancementFilter->SetSigmaStepMethodToEquispaced();
        multiScaleEnhancementFilter->SetSigmaMinimum(fSigStart);
        multiScaleEnhancementFilter->SetSigmaMaximum(fSigEnd);
        multiScaleEnhancementFilter->SetNumberOfSigmaSteps(nNumOfSigmas);
        multiScaleEnhancementFilter->Update();

#endif
        mitk::Image::Pointer outs, outs_2;

        //// set some properties
        QString qNodeName;
        //qNodeName.sprintf("double_%s_hessian_%f_%f_%d", "image", fSigStart, fSigEnd, nNumOfSigmas);
        //mitk::CastToMitkImage(multiScaleEnhancementFilter->GetOutput(), outs_2);
        //ShowProcessedDataNode(outs_2, qNodeName.toStdString());

        typedef itk::RescaleIntensityImageFilter< DoubleImageType, BinaryImageType >RescaleFilterType;
        RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
        rescaleFilter->SetInput(multiScaleEnhancementFilter->GetOutput());
        rescaleFilter->Update();
        mitk::CastToMitkImage(rescaleFilter->GetOutput(), outs);

        // set some properties
        qNodeName.sprintf("rescale_%s_hessian_%f_%f_%d", "image", fSigStart, fSigEnd, nNumOfSigmas);
        ShowProcessedDataNode(outs, qNodeName.toStdString());
    }
}



//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoRegionGrowing()
{
    mitk::Image* mitkImage = getSelectedImage();
    mitk::DataNode* mitkNode = getSelectedNode();
    QString qSTRDebug;
    if (NULL != mitkImage)
    {
        m_PointSet = m_PointListWidget->GetPointSet();
        mitk::BaseGeometry::Pointer imageGeometry = mitkImage->GetGeometry();

        if (m_PointSet->GetSize() == 0)
        {
            // no points there. Not good for region growing
            QMessageBox::information(NULL,
                "Region growing functionality",
                "Please set some seed points inside the image first.\n"
                "(hold Shift key and click left mouse button inside the image.)");
            return;
        }

        // actually perform region growing. Here we have both an image and some seed points
        FloatImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        // instantiate an ITK region growing filter, set its parameters
        typedef itk::ConnectedThresholdImageFilter<FloatImageType, FloatImageType> RegionGrowingFilterType;
        RegionGrowingFilterType::Pointer regionGrower = RegionGrowingFilterType::New();

        // determine a thresholding interval
        FloatImageType::IndexType seedIndex;
        TPixelType min(std::numeric_limits<TPixelType>::max());
        TPixelType max(std::numeric_limits<TPixelType>::min());

        mitk::PointSet::PointsContainer* points = m_PointSet->GetPointSet()->GetPoints();
        for (mitk::PointSet::PointsConstIterator pointsIterator = points->Begin(); pointsIterator != points->End(); ++pointsIterator)
        {
            // first test if this point is inside the image at all
            if (!imageGeometry->IsInside(pointsIterator.Value())){
                continue;
            }

            // convert world coordinates to image indices
            imageGeometry->WorldToIndex(pointsIterator.Value(), seedIndex);

            // get the pixel value at this point
            TPixelType currentPixelValue = itkImage->GetPixel(seedIndex);

            // adjust minimum and maximum values
            if (currentPixelValue > max)
                max = currentPixelValue;

            if (currentPixelValue < min)
                min = currentPixelValue;

            regionGrower->AddSeed(seedIndex);
        }

        // operation 
        std::stringstream strAddition("");
        // region grower
        QString qStrThreshold = m_Controls.leRGThreshold->text();
        int iThreshold = qStrThreshold.toInt();
        min -= iThreshold;
        max += iThreshold;
        strAddition << "_rg_" << min << "_" << max;

        regionGrower->SetInput(itkImage);
        // set thresholds and execute filter
        regionGrower->SetLower(min);
        regionGrower->SetUpper(max);
        regionGrower->SetReplaceValue(1);
        regionGrower->Update();

        MITK_INFO << "reigon growing done! start to dilating";

        mitk::Image::Pointer resultImage = mitk::Image::New();
        mitk::CastToMitkImage(regionGrower->GetOutput(), resultImage);

        std::string strNode = "image";
        strNode.append(strAddition.str());
        ShowProcessedDataNode(resultImage, strNode, true);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::ShowProcessedDataNode(ImageType::Pointer itkImage, std::string strNodeName, bool bSetDefaultProperties /*= false*/, mitk::DataNode::Pointer nodeParent /*= NULL*/)
{
    mitk::Image::Pointer resultImage = mitk::Image::New();
    mitk::CastToMitkImage(itkImage, resultImage);
    ShowProcessedDataNode(resultImage, strNodeName, bSetDefaultProperties, nodeParent);
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::ShowProcessedDataNode(mitk::Image::Pointer mitkImage, std::string strNodeName, bool bSetDefaultProperties /*= false*/, mitk::DataNode::Pointer nodeParent)
{
    // allocate a new node and show
    mitk::DataNode::Pointer newNode = mitk::DataNode::New();
    newNode->SetData(mitkImage);
    newNode->SetProperty("name", mitk::StringProperty::New(strNodeName));

    if (bSetDefaultProperties){
        // set some properties
        newNode->SetProperty("binary", mitk::BoolProperty::New(true));
        newNode->SetProperty("color", mitk::ColorProperty::New(1.0, 0.0, 0.0));
        newNode->SetProperty("volumerendering", mitk::BoolProperty::New(true));
        newNode->SetProperty("layer", mitk::IntProperty::New(1));
        newNode->SetProperty("opacity", mitk::FloatProperty::New(0.5));
    }

    if (NULL != nodeParent){
        this->GetDataStorage()->Add(newNode, nodeParent);
    }
    else{
        this->GetDataStorage()->Add(newNode);
    }
    // add result to data tree
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();

}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoSelectRegion()
{
    mitk::Image::Pointer pSel = getSelectedImage();
    if (NULL != pSel)
    {
        mitk::Image::Pointer mitkImage = pSel->Clone();
        ImageType::Pointer itkImage /*= ImageType::New()*/;
        mitk::CastToItkImage(mitkImage, itkImage);

        //  cast image to binary
        typedef itk::Image<unsigned char, VDIM> ImageCastOutType;
        typedef itk::CastImageFilter< ImageType, ImageCastOutType > ImageCastFilter;
        ImageCastFilter::Pointer castFilter = ImageCastFilter::New();
        castFilter->SetInput(itkImage);
        castFilter->Update();

        // compute the bounding box and save
        typedef itk::ImageMaskSpatialObject<VDIM> ImageMaskSpatialObjectType;
        ImageMaskSpatialObjectType::Pointer imageMaskSpatialObject = ImageMaskSpatialObjectType::New();
        imageMaskSpatialObject->SetImage(castFilter->GetOutput());
        ImageType::RegionType boundingBoxRegion = imageMaskSpatialObject->GetAxisAlignedBoundingBoxRegion();
        MITK_INFO << "Bounding Box Region: " << boundingBoxRegion;
        mCurrentRegion = boundingBoxRegion;

        QString strFormat;
        ImageType::SizeType imageSize = mCurrentRegion.GetSize();
        ImageType::IndexType imageIdx = mCurrentRegion.GetIndex();
        strFormat.sprintf("RegionStart:[%d %d %d], Size:[%d %d %d]", imageIdx[0], imageIdx[1], imageIdx[2], imageSize[0], imageSize[1], imageSize[2]);
        m_Controls.lbImageCrop->setText(strFormat);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoFillHole()
{
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer mitkNode = getSelectedNode();

    if (NULL != mitkImage)
    {
        mitk::Image::Pointer pNew = mitkImage->Clone();
        ImageType::Pointer itkImage;
        mitk::CastToItkImage(pNew, itkImage);

        // 		typedef itk::Image< BinaryPixelType, VDIM >   BinaryImageType;
        // 		typedef itk::Image< BinaryPixelType, 2 >	Binary2DImageType;

        typedef itk::CastImageFilter< ImageType, BinaryImageType > ImageCastFilter;
        ImageCastFilter::Pointer castFilter = ImageCastFilter::New();
        castFilter->SetInput(itkImage);
        castFilter->Update();

        typedef itk::CudaMultiplyByConstantImageFilter<BinaryImageType, BinaryImageType> MultiplyCastType;
        MultiplyCastType::Pointer multiplyFilter = MultiplyCastType::New();
        multiplyFilter->SetInput(castFilter->GetOutput());
        multiplyFilter->SetConstant(255);
        multiplyFilter->Update();

        typedef itk::SliceBySliceImageFilter< BinaryImageType, BinaryImageType> SbSFilterType;
        SbSFilterType::Pointer SbSFilter = SbSFilterType::New();

        typedef itk::BinaryFillholeImageFilter< BinaryImageType2D > FilterType;
        FilterType::Pointer HoleFillFilter = FilterType::New();
        SbSFilter->SetFilter(HoleFillFilter);
        SbSFilter->SetInput(multiplyFilter->GetOutput());
        SbSFilter->Update();

        typedef itk::CastImageFilter< BinaryImageType, ImageType> ImageCastFilter2;
        ImageCastFilter2::Pointer castFilter2 = ImageCastFilter2::New();
        castFilter2->SetInput(SbSFilter->GetOutput());
        castFilter2->Update();

        QString strNode;
        strNode.sprintf("%s%s", mitkNode->GetName().c_str(), "_fillHole");

        mitk::Image::Pointer mitkMultiplied = mitk::Image::New();
        mitk::CastToMitkImage(castFilter2->GetOutput(), mitkMultiplied);
        ShowProcessedDataNode(mitkMultiplied, strNode.toStdString(), false, mitkNode);


    }

}

////////////////////////////////////////////////////////////////////////////
//void QmitkRegionGrowingView::DoRollingBall()
//{
//	mitk::Image::Pointer mitkImage = getSelectedImage();
//	mitk::DataNode::Pointer mitkNode = getSelectedNode();
//
//	mitk::Image::Pointer mitkNew = mitkImage->Clone();
//	BinaryImageType::Pointer itkImage;
//	mitk::CastToItkImage(mitkNew, itkImage);
//
//	if (NULL != mitkImage)
//	{
//		QString strBallSize = m_Controls.leBallSize->text();
//		float iBallSize = strBallSize.toFloat();
//
//		typedef itk::CudaMultiplyByConstantImageFilter<BinaryImageType, BinaryImageType> MultiplyCastType;
//		MultiplyCastType::Pointer multiplyFilter = MultiplyCastType::New();
//		multiplyFilter->SetInput(itkImage);
//		multiplyFilter->SetConstant(255);
//		multiplyFilter->Update();
//
//		int iRadius = 8;
//		typedef itk::VotingBinaryHoleFillFloodingImageFilter< BinaryImageType2D, BinaryImageType2D> FilterType;
//		FilterType::InputSizeType radius;
//		radius.Fill(iRadius);
//
//		typedef itk::SliceBySliceImageFilter< BinaryImageType, BinaryImageType> SbSFilterType;
//		SbSFilterType::Pointer SbSFilter = SbSFilterType::New();
//
//		FilterType::Pointer filter = FilterType::New();
//		//		filter->SetInput(itkImage);
//		filter->SetRadius(radius);
//		filter->SetBackgroundValue(0);
//		filter->SetForegroundValue(255);
//		filter->SetMaximumNumberOfIterations(20);
//		//		filter->Update();
//		SbSFilter->SetFilter(filter);
//		SbSFilter->SetInput(multiplyFilter->GetOutput());
//
//		PROFILE_FUNCTION_BEGIN;
//		SbSFilter->Update();
//		PROFILE_FUNCTION_H_END(Voting);
//
//		QString strNode;
//		strNode.sprintf("%s%s_%d", mitkNode->GetName().c_str(), "_rolling", iRadius);
//		mitk::Image::Pointer mitkMultipliedRolling = mitk::Image::New();
//		mitk::CastToMitkImage(SbSFilter->GetOutput(), mitkMultipliedRolling);
//		ShowProcessedDataNode(mitkMultipliedRolling, strNode.toStdString(), false, mitkNode);
//	}
//}
//


//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoExtractRealData()
{
    mitk::DataNode::Pointer originalNode = m_Controls.cbOriginalImage->GetSelectedNode();
    mitk::DataNode::Pointer maskNode = m_Controls.cbMaskImage->GetSelectedNode();
    if (NULL == originalNode || NULL == maskNode){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }
    mitk::Image::Pointer	mitkOriginal = dynamic_cast<mitk::Image*> (originalNode->GetData());
    mitk::Image::Pointer	mitkMask = dynamic_cast<mitk::Image*>(maskNode->GetData());
    // check if images are valid
    if ((!mitkOriginal) || (!mitkMask) || (mitkOriginal->IsInitialized() == false) || (mitkMask->IsInitialized() == false)){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    QString qValue = m_Controls.leBackgroundValue->text();
    int iBackground = qValue.toInt();

    mitk::Image::Pointer mitkTarget = mitkOriginal->Clone();
    ImageType::Pointer itkImage;
    mitk::CastToItkImage(mitkTarget, itkImage);
    ImageType::Pointer itkMask;
    mitk::CastToItkImage(mitkMask, itkMask);

    itk::ImageRegionIterator<ImageType> imageRGIterator(itkMask, itkMask->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType> imageITKIterator(itkImage, itkImage->GetLargestPossibleRegion());
    while (!imageRGIterator.IsAtEnd() && !imageITKIterator.IsAtEnd())
    {
        ImageType::PixelType val = imageRGIterator.Get();
        if (val > 0){
        }
        else{
            imageITKIterator.Set(iBackground);
        }
        ++imageRGIterator;
        ++imageITKIterator;
    }

    typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
    FilterType::Pointer roiFilter = FilterType::New();
    roiFilter->SetRegionOfInterest(mCurrentRegion);
    roiFilter->SetInput(itkImage);
    roiFilter->Update();

    ImageType::Pointer itkNewImage = roiFilter->GetOutput();
    mitk::Image::Pointer mitkNewImage = mitk::Image::New();
    mitk::CastToMitkImage(itkNewImage, mitkNewImage);

    QString strNode;
    strNode.sprintf("%s_extract_bg_%d", "image"/*mitkNode->GetName().c_str()*/, iBackground);
    ShowProcessedDataNode(mitkNewImage, strNode.toStdString());
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoErosion()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        FloatImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        typedef itk::BinaryBallStructuringElement<unsigned char, VDIM>              BallType;
        typedef itk::CudaGrayscaleErodeImageFilter<FloatImageType, FloatImageType, BallType> GrayScaleErodeFilter;
        //		typedef itk::GrayscaleErodeImageFilter<ImageType, ImageType, BallType> GrayScaleErodeFilter;
        BallType binaryBall;
        int param2 = m_Controls.spinBox_2->value();
        //         const int G_KERNEL_SIZE_3 = 7;
        //         unsigned int G_KERNEL_3[G_KERNEL_SIZE_3 * G_KERNEL_SIZE_3  *G_KERNEL_SIZE_3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        BinaryImageType::SizeType radius;
        radius.Fill(param2);
        binaryBall.SetRadius(radius);
        binaryBall.CreateStructuringElement();

        GrayScaleErodeFilter::Pointer erodeFilter = GrayScaleErodeFilter::New();
        erodeFilter->SetKernel(binaryBall);
        erodeFilter->SetInput(itkImage);
        erodeFilter->Update();

        QString strNode;
        strNode.sprintf("%s_erode_%d", mitkNode->GetName().c_str(), (int)param2);
        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(erodeFilter->GetOutput(), mitkResult);
        //mStrOperation += strNode;
        ShowProcessedDataNode(mitkResult, strNode.toStdString(), false, mitkNode);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoColorLabel_2()
{
    QString strAddition = "";
    QString strFormat;

    mitk::Image::Pointer mitkImage = getSelectedImage();
    mitk::DataNode::Pointer mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        ImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        //  cast image to binary
        typedef itk::CastImageFilter< ImageType, BinaryImageType > ImageCastFilter;
        ImageCastFilter::Pointer castFilter = ImageCastFilter::New();
        castFilter->SetInput(itkImage);
        castFilter->Update();

        typedef itk::CudaMultiplyByConstantImageFilter<BinaryImageType, BinaryImageType> MultiplyCastType;
        MultiplyCastType::Pointer multiplyFilter = MultiplyCastType::New();
        multiplyFilter->SetInput(castFilter->GetOutput());
        multiplyFilter->SetConstant(255);
        multiplyFilter->Update();

        typedef itk::SliceBySliceImageFilter< BinaryImageType, ImageType> SbSFilterType;
        SbSFilterType::Pointer SbSFilter = SbSFilterType::New();
        typedef itk::ConnectedComponentImageFilter <BinaryImageType2D, ImageType2D> ConnectedComponentImageFilterType;
        ConnectedComponentImageFilterType::Pointer connectedComponentImageFilter = ConnectedComponentImageFilterType::New();
        //connectedComponentImageFilter->FullyConnectedOn();
        SbSFilter->SetFilter(connectedComponentImageFilter);
        SbSFilter->SetInput(multiplyFilter->GetOutput());
        SbSFilter->Update();

        ImageType::Pointer connected_itk_image = SbSFilter->GetOutput();
        ImageType::RegionType connected_region = connected_itk_image->GetRequestedRegion();
        ImageType::IndexType connected_index = connected_region.GetIndex();
        ImageType::SizeType connected_size = connected_region.GetSize();

        typedef itk::JoinSeriesImageFilter<ImageType2D, ImageType> JoinSeriesImageFilterType;
        JoinSeriesImageFilterType::Pointer joinFilter = JoinSeriesImageFilterType::New();
        joinFilter->SetSpacing(connected_itk_image->GetSpacing()[2]);
        joinFilter->SetOrigin(connected_itk_image->GetOrigin()[2]);

        MITK_INFO << "Size" << connected_size;
        MITK_INFO << "Index" << connected_index;

        QString strMin = m_Controls.leLabelMinVoxel->text();
        QString strMax = m_Controls.leLabelMaxVoxel->text();
        QString strRound = m_Controls.leLabelMinRoundness->text();
        float	N_MIN_R = strMin.toFloat();
        float	N_MAX_R = strMax.toFloat();
        float	N_ROUND = strRound.toFloat();
        strAddition.append(QString("_maxpixel_%1").arg(N_MAX_R));
        strAddition.append(QString("_minpixel_%1").arg(N_MIN_R));
        strAddition.append(QString("_minroundness_%1").arg(N_ROUND));

        for (int iSlice = 0; iSlice < connected_size[2]; iSlice++)
        {
            typedef itk::ExtractImageFilter<ImageType, ImageType2D> ExtractFilter;
            ExtractFilter::Pointer	extract_filter = ExtractFilter::New();
            ImageType::RegionType	extract_region;
            ImageType::IndexType	extract_index = connected_index;
            ImageType::SizeType		extract_size = connected_size;
            extract_index[2] = iSlice;
            extract_size[2] = 0;
            extract_region.SetIndex(extract_index);
            extract_region.SetSize(extract_size);
            extract_filter->SetInput(connected_itk_image);
            extract_filter->SetExtractionRegion(extract_region);
            extract_filter->SetDirectionCollapseToIdentity();
            extract_filter->Update();

            typedef itk::LabelImageToShapeLabelMapFilter<ImageType2D> ShapeLabelMapFilter_2D;
            ShapeLabelMapFilter_2D::Pointer shape_filter = ShapeLabelMapFilter_2D::New();
            shape_filter->SetInput(extract_filter->GetOutput());
            shape_filter->Update();

            //  By default, objects with an attribute value smaller than Lamba are removed. 
            typedef itk::ShapeOpeningLabelMapFilter< ShapeLabelMapFilter_2D::OutputImageType > ShapeOpeningLabelMapFilterType_2D;
            ShapeOpeningLabelMapFilterType_2D::Pointer shapeOpeningLabelMapFilter0 = ShapeOpeningLabelMapFilterType_2D::New();
            shapeOpeningLabelMapFilter0->SetInput(shape_filter->GetOutput());
            shapeOpeningLabelMapFilter0->SetLambda(N_MIN_R);
            //shapeOpeningLabelMapFilter0->ReverseOrderingOn();
            shapeOpeningLabelMapFilter0->SetAttribute(ShapeOpeningLabelMapFilterType_2D::LabelObjectType::NUMBER_OF_PIXELS);
            shapeOpeningLabelMapFilter0->Update();

            ShapeOpeningLabelMapFilterType_2D::Pointer shapeOpeningLabelMapFilter1 = ShapeOpeningLabelMapFilterType_2D::New();
            shapeOpeningLabelMapFilter1->SetInput(shapeOpeningLabelMapFilter0->GetOutput());
            shapeOpeningLabelMapFilter1->SetLambda(N_MAX_R);
            shapeOpeningLabelMapFilter1->ReverseOrderingOn();
            shapeOpeningLabelMapFilter1->SetAttribute(ShapeOpeningLabelMapFilterType_2D::LabelObjectType::NUMBER_OF_PIXELS);
            shapeOpeningLabelMapFilter1->Update();

            ShapeOpeningLabelMapFilterType_2D::Pointer shapeOpeningLabelMapFilter2 = ShapeOpeningLabelMapFilterType_2D::New();
            shapeOpeningLabelMapFilter2->SetInput(shapeOpeningLabelMapFilter1->GetOutput());
            shapeOpeningLabelMapFilter2->SetLambda(N_ROUND);
            shapeOpeningLabelMapFilter2->SetAttribute(ShapeOpeningLabelMapFilterType_2D::LabelObjectType::ROUNDNESS);
            shapeOpeningLabelMapFilter2->Update();

            // Create a label image
            typedef itk::LabelMapToLabelImageFilter<ShapeOpeningLabelMapFilterType_2D::OutputImageType, ImageType2D> LabelMapToLabelImageFilterType;
            LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
            labelMapToLabelImageFilter->SetInput(shapeOpeningLabelMapFilter2->GetOutput());
            labelMapToLabelImageFilter->Update();

            joinFilter->SetInput(iSlice, labelMapToLabelImageFilter->GetOutput());
        }
        joinFilter->Update();

        //-------------------------------------------------------------------------
        // 2d filtered color image
        typedef itk::LabelToRGBImageFilter<ImageType, RGBImageType> RGBFilterType;
        RGBFilterType::Pointer rgbFilter = RGBFilterType::New();
        rgbFilter->SetInput(joinFilter->GetOutput());
        rgbFilter->Update();

        strFormat.sprintf("2d_connected_%s", strAddition.toStdString().c_str());
        mitk::Image::Pointer mitkConnected;
        mitk::CastToMitkImage(rgbFilter->GetOutput(), mitkConnected);
        ShowProcessedDataNode(mitkConnected, strFormat.toStdString(), false, mitkNode);

        //-------------------------------------------------------------------------
        /// 2d filtered mask
        typedef itk::BinaryThresholdImageFilter<ImageType, BinaryImageType> BinaryThresholdFilter;
        BinaryThresholdFilter::Pointer binary_threshold = BinaryThresholdFilter::New();
        binary_threshold->SetInput(joinFilter->GetOutput());
        binary_threshold->SetOutsideValue(0);
        binary_threshold->SetInsideValue(1);
        binary_threshold->SetLowerThreshold(0);
        binary_threshold->Update();

        strFormat.sprintf("2d_mask_connected_%s", strAddition.toStdString().c_str());
        mitk::Image::Pointer mitkConnectedMask;
        mitk::CastToMitkImage(binary_threshold->GetOutput(), mitkConnectedMask);
        ShowProcessedDataNode(mitkConnectedMask, strFormat.toStdString(), false, mitkNode);

        // transform to [0-255]
        typedef itk::MultiplyImageFilter<BinaryImageType, BinaryImageType> MultiplyFilter_3D;
        MultiplyFilter_3D::Pointer mulfilter_3d = MultiplyFilter_3D::New();
        mulfilter_3d->SetInput(binary_threshold->GetOutput());
        mulfilter_3d->SetConstant(255);
        mulfilter_3d->Update();

        QString strMin_3D = m_Controls.leMinVoxel->text();
        QString strMax_3D = m_Controls.leMaxVoxel->text();
        QString strMinSphere_3D = m_Controls.leMinSphereness->text();
        float	N_MIN_R_3D = strMin_3D.toFloat();
        float	N_MAX_R_3D = strMax_3D.toFloat();
        float	N_MIN_SPH_3D = strMinSphere_3D.toFloat();
        strAddition.clear();
        strAddition.append(QString("_maxvoxel_%1").arg(N_MAX_R_3D));
        strAddition.append(QString("_minvoxel_%1").arg(N_MIN_R_3D));
        strAddition.append(QString("_minsphereness_%1").arg(N_MIN_SPH_3D));

        typedef itk::BinaryImageToShapeLabelMapFilter<BinaryImageType> ShapeLabelMapFilter_3D;
        ShapeLabelMapFilter_3D::Pointer shape_filter_3d = ShapeLabelMapFilter_3D::New();
        shape_filter_3d->SetInput(mulfilter_3d->GetOutput());
        shape_filter_3d->Update();

        // By default, objects with an attribute value smaller than Lamba are removed. 
        typedef itk::ShapeOpeningLabelMapFilter< ShapeLabelMapFilter_3D::OutputImageType > ShapeOpeningLabelMapFilterType_3D;
        ShapeOpeningLabelMapFilterType_3D::Pointer shapeOpeningLabelMapFilter0_3D = ShapeOpeningLabelMapFilterType_3D::New();
        shapeOpeningLabelMapFilter0_3D->SetInput(shape_filter_3d->GetOutput());
        shapeOpeningLabelMapFilter0_3D->SetLambda(N_MIN_R_3D);
        shapeOpeningLabelMapFilter0_3D->SetAttribute(ShapeOpeningLabelMapFilterType_3D::LabelObjectType::NUMBER_OF_PIXELS);
        shapeOpeningLabelMapFilter0_3D->Update();

        ShapeOpeningLabelMapFilterType_3D::Pointer shapeOpeningLabelMapFilter1_3D = ShapeOpeningLabelMapFilterType_3D::New();
        shapeOpeningLabelMapFilter1_3D->SetInput(shapeOpeningLabelMapFilter0_3D->GetOutput());
        shapeOpeningLabelMapFilter1_3D->SetLambda(N_MAX_R_3D);
        shapeOpeningLabelMapFilter1_3D->ReverseOrderingOn();
        shapeOpeningLabelMapFilter1_3D->SetAttribute(ShapeOpeningLabelMapFilterType_3D::LabelObjectType::NUMBER_OF_PIXELS);
        shapeOpeningLabelMapFilter1_3D->Update();

        ShapeOpeningLabelMapFilterType_3D::Pointer shapeOpeningLabelMapFilter2_3D = ShapeOpeningLabelMapFilterType_3D::New();
        shapeOpeningLabelMapFilter2_3D->SetInput(shapeOpeningLabelMapFilter1_3D->GetOutput());
        shapeOpeningLabelMapFilter2_3D->SetLambda(N_MIN_SPH_3D);
        shapeOpeningLabelMapFilter2_3D->SetAttribute(ShapeOpeningLabelMapFilterType_3D::LabelObjectType::ROUNDNESS);
        shapeOpeningLabelMapFilter2_3D->Update();

        typedef itk::LabelMapToLabelImageFilter<ShapeLabelMapFilter_3D::OutputImageType, ImageType> LabelMapToLabelImageFilterType_3D;
        LabelMapToLabelImageFilterType_3D::Pointer labelMapToLabelImageFilter_3d = LabelMapToLabelImageFilterType_3D::New();
        labelMapToLabelImageFilter_3d->SetInput(shapeOpeningLabelMapFilter2_3D->GetOutput());
        labelMapToLabelImageFilter_3d->Update();

        typedef itk::ShapeLabelObject<itk::SizeValueType, VDIM>	ShapeLabelObjectType;
        typedef itk::LabelMap< ShapeLabelObjectType >	ShapeLabelMapType;
        ShapeLabelMapType::Pointer shape_labelMap = shapeOpeningLabelMapFilter2_3D->GetOutput();

        // 	int szLabelMap0 = statistics_labelMap->GetNumberOfLabelObjects();
        int szLabelMap = shape_labelMap->GetNumberOfLabelObjects();
        // 	MITK_INFO << "There are " << szLabelMap0 << " labels in statistics label.";
        MITK_INFO << "There are " << szLabelMap << " labels in shape label.";

        //-------------------------------------------------------------------------
        // 3d filtered color image
        RGBFilterType::Pointer rgbFilter_3d = RGBFilterType::New();
        rgbFilter_3d->SetInput(labelMapToLabelImageFilter_3d->GetOutput());
        rgbFilter_3d->Update();
        strFormat.sprintf("3d_connected_%s", strAddition.toStdString().c_str());
        mitk::Image::Pointer mitkConnected_3D;
        mitk::CastToMitkImage(rgbFilter_3d->GetOutput(), mitkConnected_3D);
        ShowProcessedDataNode(mitkConnected_3D, strFormat.toStdString(), false, mitkNode);

        //-------------------------------------------------------------------------
        // 3d filtered mask image
        typedef itk::LabelMapToBinaryImageFilter<ShapeLabelMapFilter_3D::OutputImageType, BinaryImageType> LabelMapToBinaryFilter;
        LabelMapToBinaryFilter::Pointer labeltobinaryFilter = LabelMapToBinaryFilter::New();
        labeltobinaryFilter->SetInput(shapeOpeningLabelMapFilter2_3D->GetOutput());
        labeltobinaryFilter->Update();

        strFormat.sprintf("3d_mask_connected_%s", strAddition.toStdString().c_str());
        mitk::Image::Pointer mitkConnectedMask_3D;
        mitk::CastToMitkImage(labeltobinaryFilter->GetOutput(), mitkConnectedMask_3D);
        ShowProcessedDataNode(mitkConnectedMask_3D, strFormat.toStdString(), false, mitkNode);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoColorLabel()
{
    QString strAddition = "";
    QString strFormat;

    //mitk::Image::Pointer mitkImage = getSelectedImage();
    //mitk::DataNode::Pointer mitkNode = getSelectedNode();

    // validate data 
    mitk::DataNode::Pointer mitkNode = m_Controls.cbStatisticImage->GetSelectedNode();
    mitk::DataNode::Pointer featureNode = m_Controls.cbStatisticFeatureImage->GetSelectedNode();

    if (NULL == mitkNode || NULL == featureNode){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    mitk::Image::Pointer	mitkImage = dynamic_cast<mitk::Image*> (mitkNode->GetData());
    mitk::Image::Pointer	mitkFeatureImage = dynamic_cast<mitk::Image*>(featureNode->GetData());
    // check if images are valid
    if ((!mitkImage) || (!mitkFeatureImage) || (mitkImage->IsInitialized() == false) || (mitkFeatureImage->IsInitialized() == false)){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    // image conversion
    BinaryImageType::Pointer itkImage;
    mitk::CastToItkImage(mitkImage, itkImage);

    ImageType::Pointer itkFeatureImage;
    mitk::CastToItkImage(mitkFeatureImage, itkFeatureImage);

    typedef itk::SliceBySliceImageFilter< BinaryImageType, ImageType> SbSFilterType;
    SbSFilterType::Pointer SbSFilter = SbSFilterType::New();
    // using connectness to label image
    typedef itk::ConnectedComponentImageFilter <BinaryImageType2D, ImageType2D> ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connectedComponentImageFilter = ConnectedComponentImageFilterType::New();
    //     connectedComponentImageFilter->FullyConnectedOn();
    SbSFilter->SetFilter(connectedComponentImageFilter);
    SbSFilter->SetInput(/*multiplyFilter->GetOutput()*/itkImage);
    SbSFilter->Update();

    ImageType::Pointer connected_itk_image = SbSFilter->GetOutput();
    ImageType::RegionType connected_region = connected_itk_image->GetRequestedRegion();
    ImageType::IndexType connected_index = connected_region.GetIndex();
    ImageType::SizeType connected_size = connected_region.GetSize();

    typedef itk::JoinSeriesImageFilter<ImageType2D, ImageType> JoinSeriesImageFilterType;
    JoinSeriesImageFilterType::Pointer joinFilter = JoinSeriesImageFilterType::New();
    joinFilter->SetSpacing(connected_itk_image->GetSpacing()[2]);
    joinFilter->SetOrigin(connected_itk_image->GetOrigin()[2]);

    MITK_INFO << "Size" << connected_size;
    MITK_INFO << "Index" << connected_index;

    QString strMin = m_Controls.leLabelMinVoxel->text();
    QString strMax = m_Controls.leLabelMaxVoxel->text();
    QString strRound = m_Controls.leLabelMinRoundness->text();
    float	N_MIN_R = strMin.toFloat();
    float	N_MAX_R = strMax.toFloat();
    float	N_ROUND = strRound.toFloat();
    strAddition.append(QString("_maxpixel_%1").arg(N_MAX_R));
    strAddition.append(QString("_minpixel_%1").arg(N_MIN_R));
    strAddition.append(QString("_minroundness_%1").arg(N_ROUND));

    // reduce image slice by slice
    for (int iSlice = 0; iSlice < connected_size[2]; iSlice++)
    {
        typedef itk::ExtractImageFilter<ImageType, ImageType2D> ExtractFilter;
        ExtractFilter::Pointer	extract_filter = ExtractFilter::New();
        ImageType::RegionType	extract_region;
        ImageType::IndexType	extract_index = connected_index;
        ImageType::SizeType		extract_size = connected_size;
        extract_index[2] = iSlice;
        extract_size[2] = 0;
        extract_region.SetIndex(extract_index);
        extract_region.SetSize(extract_size);
        extract_filter->SetInput(connected_itk_image);
        extract_filter->SetExtractionRegion(extract_region);
        extract_filter->SetDirectionCollapseToIdentity();
        extract_filter->Update();

        // convert label to shape image to evaluate the shape properties
        typedef itk::LabelImageToShapeLabelMapFilter<ImageType2D> ShapeLabelMapFilter_2D;
        ShapeLabelMapFilter_2D::Pointer shape_filter = ShapeLabelMapFilter_2D::New();
        shape_filter->SetInput(extract_filter->GetOutput());
        shape_filter->Update();

        //  By default, objects with an attribute value smaller than Lamba are removed. 
        typedef itk::ShapeOpeningLabelMapFilter< ShapeLabelMapFilter_2D::OutputImageType > ShapeOpeningLabelMapFilterType_2D;
        ShapeOpeningLabelMapFilterType_2D::Pointer shapeOpeningLabelMapFilter0 = ShapeOpeningLabelMapFilterType_2D::New();
        shapeOpeningLabelMapFilter0->SetInput(shape_filter->GetOutput());
        shapeOpeningLabelMapFilter0->SetLambda(N_MIN_R);
        //shapeOpeningLabelMapFilter0->ReverseOrderingOn();
        shapeOpeningLabelMapFilter0->SetAttribute(ShapeOpeningLabelMapFilterType_2D::LabelObjectType::NUMBER_OF_PIXELS);
        shapeOpeningLabelMapFilter0->Update();

        ShapeOpeningLabelMapFilterType_2D::Pointer shapeOpeningLabelMapFilter1 = ShapeOpeningLabelMapFilterType_2D::New();
        shapeOpeningLabelMapFilter1->SetInput(shapeOpeningLabelMapFilter0->GetOutput());
        shapeOpeningLabelMapFilter1->SetLambda(N_MAX_R);
        shapeOpeningLabelMapFilter1->ReverseOrderingOn();
        shapeOpeningLabelMapFilter1->SetAttribute(ShapeOpeningLabelMapFilterType_2D::LabelObjectType::NUMBER_OF_PIXELS);
        shapeOpeningLabelMapFilter1->Update();

        ShapeOpeningLabelMapFilterType_2D::Pointer shapeOpeningLabelMapFilter2 = ShapeOpeningLabelMapFilterType_2D::New();
        shapeOpeningLabelMapFilter2->SetInput(shapeOpeningLabelMapFilter1->GetOutput());
        shapeOpeningLabelMapFilter2->SetLambda(N_ROUND);
        shapeOpeningLabelMapFilter2->SetAttribute(ShapeOpeningLabelMapFilterType_2D::LabelObjectType::ROUNDNESS);
        shapeOpeningLabelMapFilter2->Update();

        // convert back to label image
        typedef itk::LabelMapToLabelImageFilter<ShapeOpeningLabelMapFilterType_2D::OutputImageType, ImageType2D> LabelMapToLabelImageFilterType;
        LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
        labelMapToLabelImageFilter->SetInput(shapeOpeningLabelMapFilter2->GetOutput());
        labelMapToLabelImageFilter->Update();

        joinFilter->SetInput(iSlice, labelMapToLabelImageFilter->GetOutput());
    }
    joinFilter->Update();

    //-------------------------------------------------------------------------
    // show 2d filtered color image
    typedef itk::LabelToRGBImageFilter<ImageType, RGBImageType> RGBFilterType;
    RGBFilterType::Pointer rgbFilter = RGBFilterType::New();
    rgbFilter->SetInput(joinFilter->GetOutput());
    rgbFilter->Update();

    strFormat.sprintf("sort_2d_%s", strAddition.toStdString().c_str());
    mitk::Image::Pointer mitkConnected;
    mitk::CastToMitkImage(rgbFilter->GetOutput(), mitkConnected);
    ShowProcessedDataNode(mitkConnected, strFormat.toStdString(), false, mitkNode);

    /// show 2d filtered mask
    typedef itk::BinaryThresholdImageFilter<ImageType, BinaryImageType> BinaryThresholdFilter;
    BinaryThresholdFilter::Pointer binary_threshold = BinaryThresholdFilter::New();
    binary_threshold->SetInput(joinFilter->GetOutput());
    binary_threshold->SetOutsideValue(0);
    binary_threshold->SetInsideValue(1);
    binary_threshold->SetLowerThreshold(0);
    binary_threshold->Update();

    strFormat.sprintf("sort_2d_mask_%s", strAddition.toStdString().c_str());
    mitk::Image::Pointer mitkConnectedMask;
    mitk::CastToMitkImage(binary_threshold->GetOutput(), mitkConnectedMask);
    ShowProcessedDataNode(mitkConnectedMask, strFormat.toStdString(), false, mitkNode);

    //-------------------------------------------------------------------------
    // do 3d processing
    // convert the image to a label image using connect
    typedef itk::ConnectedComponentImageFilter <BinaryImageType, ImageType> ConnectedComponentImageFilterType_3D;
    ConnectedComponentImageFilterType_3D::Pointer connectedComponentImageFilter_3D = ConnectedComponentImageFilterType_3D::New();
    //     connectedComponentImageFilter_3D->FullyConnectedOn();
    connectedComponentImageFilter_3D->SetInput(/*mulfilter_3d*/binary_threshold->GetOutput());
    connectedComponentImageFilter_3D->Update();

    // statistical analysis
    typedef itk::LabelImageToStatisticsLabelMapFilter<ImageType, ImageType> StatisticsLabelMapFilter_3D;
    StatisticsLabelMapFilter_3D::Pointer statistics_filter_3d = StatisticsLabelMapFilter_3D::New();
    statistics_filter_3d->SetInput(connectedComponentImageFilter_3D->GetOutput());
    statistics_filter_3d->SetFeatureImage(itkFeatureImage);
    statistics_filter_3d->Update();

    // shape analysis
    typedef itk::LabelImageToShapeLabelMapFilter<ImageType> ShapeLabelMapFilter_3D;
    ShapeLabelMapFilter_3D::Pointer shape_filter_3d = ShapeLabelMapFilter_3D::New();
    shape_filter_3d->SetInput(connectedComponentImageFilter_3D->GetOutput());
    shape_filter_3d->Update();

    // analysis according to shape & statistic properties
    typedef itk::StatisticsLabelObject<TPixelType, VDIM>StatisticsLabelObjectType;
    typedef itk::LabelMap< StatisticsLabelObjectType >	StatisticsLabelMapType;
    StatisticsLabelMapType::Pointer statistics_labelMap = statistics_filter_3d->GetOutput();

    typedef itk::ShapeLabelObject<TPixelType, VDIM>	ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType >	ShapeLabelMapType;
    ShapeLabelMapType::Pointer shape_labelMap = shape_filter_3d->GetOutput();

    int szLabelMap0 = statistics_labelMap->GetNumberOfLabelObjects();
    int szLabelMap = shape_labelMap->GetNumberOfLabelObjects();

    m_vLabeledRegions.clear();
    //     m_vLabeledRegions.resize(szLabelMap);

    for (unsigned int n = 0; n < szLabelMap; ++n){
        ShapeLabelObjectType*	labelObject = shape_labelMap->GetNthLabelObject(n);
        ImageType::RegionType	bound = labelObject->GetBoundingBox();
        m_vLabeledRegions.push_back(bound);
    }

    // 	MITK_INFO << "There are " << szLabelMap0 << " labels in statistics label.";
    MITK_INFO << "There are " << szLabelMap << " labels in shape label.";

    QString strMin_3D = m_Controls.leMinVoxel->text();
    QString strMax_3D = m_Controls.leMaxVoxel->text();
    QString strSelect = m_Controls.leLabelCount->text();
    QString strMinSphere_3D = m_Controls.leMinSphereness->text();
    float	N_MIN_R_3D = strMin_3D.toFloat();
    float	N_MAX_R_3D = strMax_3D.toFloat();
    float	N_MIN_SPH_3D = strMinSphere_3D.toFloat();
    int		nMinR3D = N_MIN_R_3D;
    int		nMaxR3D = N_MAX_R_3D;
    int		nSelect = strSelect.toInt();

    double* dwVotingScore = new double[szLabelMap];

#define DECLARE_VOTING_VAR(name) double*dwVotingScore##name=new double[szLabelMap];\
for (int iVotingScore##name=0; iVotingScore##name<szLabelMap; iVotingScore##name++){\
dwVotingScore##name[iVotingScore##name] = 0.0;}\
double dwMax##name = -999999999.0;\
double dwMin##name = 999999999.0;\
for (unsigned int n = 0; n < szLabelMap; ++n){\
ShapeLabelObjectType*		labelObject = shape_labelMap->GetNthLabelObject(n);\
StatisticsLabelObjectType*		statisticObject = statistics_labelMap->GetNthLabelObject(n);\
ImageType::RegionType	bound = labelObject->GetBoundingBox();\
ImageType::SizeType bound_size = bound.GetSize();\
ImageType::IndexType bound_index = bound.GetIndex();

#define NORMALIZE_VOTING_VAR(name)\
if(dwVotingScore##name[n] > dwMax##name)\
dwMax##name = dwVotingScore##name[n];\
if (dwVotingScore##name[n] < dwMin##name)\
dwMin##name = dwVotingScore##name[n];\
}\
for(unsigned int n= 0;n<szLabelMap;++n){\
dwVotingScore##name[n] = (dwVotingScore##name[n]-dwMin##name)/(dwMax##name-dwMin##name);}

#if 1
    QFile qTEST_Original("statistics_orig.txt");
    qTEST_Original.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate);
    QTextStream tsTEST_Original(&qTEST_Original);
    tsTEST_Original.setRealNumberNotation(QTextStream::FixedNotation);
    tsTEST_Original.setRealNumberPrecision(6);
    tsTEST_Original << "Index\t\t" << "Round\t\t" << "NumPix\t\t" << "Skew\t\t" << "Elong\t\t" << "Start\t\t" << "Size\t\t" << "\n";
    for (int i = 0; i < szLabelMap; i++)
    {
        ShapeLabelObjectType*		    labelObject = shape_labelMap->GetNthLabelObject(i);
        StatisticsLabelObjectType*		statisticObject = statistics_labelMap->GetNthLabelObject(i);
        ImageType::RegionType	bound = labelObject->GetBoundingBox();
        ImageType::SizeType     bound_size = bound.GetSize();
        ImageType::IndexType    bound_index = bound.GetIndex();
        tsTEST_Original << i << "\t\t" << labelObject->GetRoundness() << "\t\t" << labelObject->GetNumberOfPixels() << "\t\t"
            << statisticObject->GetSkewness() << "\t\t" << statisticObject->GetElongation() << "\t\t"
            << bound_index[0] + bound_size[0] / 2 << " " << bound_index[1] + bound_size[1] / 2 << " " << bound_index[2] + bound_size[2] / 2 << "\t\t"
            << bound_size[0] << " " << bound_size[1] << " " << bound_size[2] << "\n";
    }
    qTEST_Original.close();
#endif


    DECLARE_VOTING_VAR(Roundness);
    // 1. roundness
    double dwRoundness = labelObject->GetRoundness();
    dwVotingScoreRoundness[n] = labelObject->GetRoundness();
    if (dwRoundness < N_MIN_SPH_3D) {
        dwVotingScoreRoundness[n] = 0.0;
    }
    NORMALIZE_VOTING_VAR(Roundness);

    // 2. size
    DECLARE_VOTING_VAR(NumPixels);
    unsigned int nPixels = labelObject->GetNumberOfPixels();
    if (nPixels > nMinR3D && nPixels < nMaxR3D){
        dwVotingScoreNumPixels[n] = (1.0 + ((double)nPixels - nMinR3D) / (nMaxR3D - nMinR3D));
    }
    else if (nPixels >= nMaxR3D){
        // dwVotingScore[n] *= (1.0 / (exp((nPixels - nMaxR3D) / (float)nMaxR3D)));
        //dwVotingScoreNumPixels[n] = -1.0;
        if (nPixels < 2 * nMaxR3D)
            dwVotingScoreNumPixels[n] = (1.0 - (nPixels - nMaxR3D) / (double)nMaxR3D);
        else
            dwVotingScoreNumPixels[n] = 0.0;
    }
    else if (nPixels <= nMinR3D){
        //dwVotingScore[n] *= (1.0 / (exp((nMinR3D - nPixels) / (float)nMinR3D)));
        dwVotingScoreNumPixels[n] = (1.0 - (nMinR3D - nPixels) / (double)nMinR3D);
        //dwVotingScoreNumPixels[n] = -1.0;
    }
    NORMALIZE_VOTING_VAR(NumPixels);

    //
    DECLARE_VOTING_VAR(Slices);
    double N_MIN_SLICES_LOW = 3.0;
    double N_MIN_SLICES_HIGH = 15.0;
    int szSlices = bound_size[2];
    if (szSlices < N_MIN_SLICES_LOW){
        dwVotingScoreSlices[n] = 0.0;
    }
    else if ((szSlices >= N_MIN_SLICES_LOW) && (szSlices <= N_MIN_SLICES_HIGH)){
        dwVotingScoreSlices[n] = (1.0 + ((double)szSlices - N_MIN_SLICES_LOW) / (N_MIN_SLICES_HIGH - N_MIN_SLICES_LOW));
    }
    else if (szSlices > N_MIN_SLICES_HIGH){
        dwVotingScoreSlices[n] = (1.0 + 1.0 / exp(((double)szSlices - N_MIN_SLICES_HIGH) / (double)N_MIN_SLICES_HIGH));
    }
    NORMALIZE_VOTING_VAR(Slices);

    DECLARE_VOTING_VAR(Elongation);
    // 3. Elongation 
    double DW_ACCEPT_ELONGATION = 3.0;
    double dwDeltaElongation = abs(labelObject->GetElongation() - 1.0);
    if (dwDeltaElongation < DW_ACCEPT_ELONGATION){
        dwVotingScoreElongation[n] = (1.0 + (DW_ACCEPT_ELONGATION - (double)dwDeltaElongation) / DW_ACCEPT_ELONGATION);
    }
    else{
        dwVotingScoreElongation[n] = 0.1;
    }
    NORMALIZE_VOTING_VAR(Elongation);


    // 4. distribution
    DECLARE_VOTING_VAR(Skewness);
    double DW_ACCEPT_KURTOSIS = 1.0;
    double DW_ACCEPT_SKEWNESS = 3.0;
    double dwKurtosis = statisticObject->GetKurtosis();
    double dwSkewness = statisticObject->GetSkewness();
    double dwDeltaKurtosis = abs(dwKurtosis);
    double dwDeltaSkewness = abs(dwSkewness);

    //if (dwDeltaKurtosis < DW_ACCEPT_KURTOSIS){
    //	dwVotingScore[n] *= (1.0 + 1.0 - dwDeltaKurtosis/DW_ACCEPT_KURTOSIS);
    //}
    //else{
    //	dwVotingScore[n] *= (1.0 / exp(dwDeltaKurtosis - DW_ACCEPT_KURTOSIS));
    //}

    if (dwDeltaSkewness < DW_ACCEPT_SKEWNESS){
        dwVotingScoreSkewness[n] = (1.0 + (DW_ACCEPT_SKEWNESS - (double)dwDeltaSkewness) / DW_ACCEPT_SKEWNESS);
    }
    else{
        dwVotingScoreSkewness[n] = 0.1;
    }
    NORMALIZE_VOTING_VAR(Skewness);


    for (int i = 0; i < szLabelMap; i++){
        StatisticsLabelObjectType*		statisticObject = statistics_labelMap->GetNthLabelObject(i);
        double dwPixelMean = statisticObject->GetMean();
        TPixelType pixMaxmium = statisticObject->GetMaximum();
        dwVotingScore[i] =
            0.15 * dwVotingScoreNumPixels[i] +
            0.25* dwVotingScoreRoundness[i] +
            0.2 * dwVotingScoreSlices[i] +
            0.2 * dwVotingScoreSkewness[i] +
            0.2 * dwVotingScoreElongation[i];

        if (dwPixelMean > -200.0 && pixMaxmium > 800.0){
            if (dwVotingScoreRoundness[i] > 0.3){
                //dwVotingScore[i] = 10.0;
                dwVotingScore[i] *= 2.0;
                MITK_WARN << "Reweighting label " << i << ", dwPixelMean: " << dwPixelMean << ", pixMax: " << pixMaxmium;
            }
        }
    }
#if 1
    QFile qTEST("statistics_norm.txt");
    qTEST.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate);
    QTextStream tsTEST(&qTEST);
    tsTEST.setRealNumberNotation(QTextStream::FixedNotation);
    tsTEST.setRealNumberPrecision(6);
    tsTEST << "Max: " << dwMaxRoundness << "\t" << dwMaxNumPixels << "\t" << dwMaxSlices << "\t" << dwMaxSkewness << "\t" << dwMaxElongation << "\n";
    tsTEST << "Min: " << dwMinRoundness << "\t" << dwMinNumPixels << "\t" << dwMinSlices << "\t" << dwMinSkewness << "\t" << dwMinElongation << "\n";

    tsTEST << "Index\t\t" << "Round\t\t" << "NumPix\t\t" << "Slices\t\t" << "Skew\t\t" << "Elong\t\t" << "Total" << "\n";
    for (int i = 0; i < szLabelMap; i++){
        tsTEST << i << "\t\t" << dwVotingScoreRoundness[i] << "\t\t" << dwVotingScoreNumPixels[i] << "\t\t" << dwVotingScoreSlices[i] << "\t\t"
            << dwVotingScoreSkewness[i] << "\t\t" << dwVotingScoreElongation[i] << "\t\t" << dwVotingScore[i] << "\n";
    }
    qTEST.close();
#endif

    typedef std::multimap<double, int/*, std::greater<double>*/> SORT_MAP;
    SORT_MAP sort_map;
    for (unsigned int n = 0; n < szLabelMap; ++n){
        ShapeLabelObjectType*	labelObject = shape_labelMap->GetNthLabelObject(n);
        sort_map.insert(std::make_pair(dwVotingScore[n], labelObject->GetLabel()));
    }

    SORT_MAP::iterator sort_map_it = sort_map.begin(), sort_map_end = sort_map.end();
    if (szLabelMap > nSelect)
    {
        int nCounter = szLabelMap - nSelect;
        for (int i = nCounter; i > 0 && sort_map_it != sort_map_end; i--, sort_map_it++){
            shape_labelMap->RemoveLabel(sort_map_it->second);
        }
    }
    int szSortedLabelMap = shape_labelMap->GetNumberOfLabelObjects();
    MITK_INFO << "LabelMap after Sort " << szSortedLabelMap;

    // transform to label image
    typedef itk::LabelMapToLabelImageFilter<ShapeLabelMapFilter_3D::OutputImageType, ImageType> LabelMapToLabelImageFilterType_3D;
    LabelMapToLabelImageFilterType_3D::Pointer labelMapToLabelImageFilter_3d = LabelMapToLabelImageFilterType_3D::New();
    labelMapToLabelImageFilter_3d->SetInput(shape_labelMap);
    labelMapToLabelImageFilter_3d->Update();

    // show rgb image
    RGBFilterType::Pointer rgbFilter_3d = RGBFilterType::New();
    rgbFilter_3d->SetInput(labelMapToLabelImageFilter_3d->GetOutput());
    rgbFilter_3d->Update();
    strFormat.sprintf("sort_%d_pix_%d_%d_round_%f", nSelect, nMinR3D, nMaxR3D, N_MIN_SPH_3D);
    mitk::Image::Pointer mitkConnected_3D;
    mitk::CastToMitkImage(rgbFilter_3d->GetOutput(), mitkConnected_3D);
    ShowProcessedDataNode(mitkConnected_3D, strFormat.toStdString(), false, mitkNode);

    // show mask
    BinaryThresholdFilter::Pointer binary_threshold_3d = BinaryThresholdFilter::New();
    binary_threshold_3d->SetInput(labelMapToLabelImageFilter_3d->GetOutput());
    binary_threshold_3d->SetOutsideValue(0);
    binary_threshold_3d->SetInsideValue(1);
    binary_threshold_3d->SetLowerThreshold(0);
    binary_threshold_3d->Update();
    strFormat.sprintf("mask_sort_%d_pix_%d_%d_round_%f", nSelect, nMinR3D, nMaxR3D, N_MIN_SPH_3D);
    mitk::Image::Pointer mitkConnectedMask3D;
    mitk::CastToMitkImage(binary_threshold_3d->GetOutput(), mitkConnectedMask3D);
    ShowProcessedDataNode(mitkConnectedMask3D, strFormat.toStdString(), false, mitkNode);

    SAFE_DELETE(dwVotingScore);
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::collectRegions(BinaryImageType::Pointer itkBinaryImage, VVisualizeNodules& vOutput)
{
    typedef itk::BinaryImageToShapeLabelMapFilter<BinaryImageType> ShapeLabelMapFilter_3D;
    ShapeLabelMapFilter_3D::Pointer shape_filter_3d = ShapeLabelMapFilter_3D::New();
    shape_filter_3d->SetInput(/*mulfilter_3d->GetOutput()*/itkBinaryImage);
    shape_filter_3d->Update();

    typedef itk::ShapeLabelObject<itk::SizeValueType, VDIM>	ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType >			LabelMapType;
    LabelMapType::Pointer labelMap = shape_filter_3d->GetOutput();
    unsigned int number_of_objects = labelMap->GetNumberOfLabelObjects();
    MITK_INFO << "There are " << number_of_objects << " labels.";
    // Retrieve all attributes
    vOutput.resize(number_of_objects);
    for (unsigned int n = 0; n < number_of_objects; ++n)
    {
        ShapeLabelObjectType*	labelObject = labelMap->GetNthLabelObject(n);
        ImageType::RegionType	bound = labelObject->GetBoundingBox();
        for (unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++){
            vOutput[n].nodule_indexes.push_back(labelObject->GetIndex(pixelId));
        }
        vOutput[n].nodule_size = vOutput[n].nodule_indexes.size();
        vOutput[n].nodule_region = bound;
    }
}


// #define USE_NODULE_CENTERS
//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoOutputForecast()
{
    // for debugging formation
    QString strQDebug;

    // validate data 
    mitk::DataNode::Pointer textureNode = m_Controls.cbOriginalImage_2->GetSelectedNode();
    mitk::DataNode::Pointer maskNode = m_Controls.cbMaskImage_2->GetSelectedNode();

    if (NULL == textureNode || NULL == maskNode){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    mitk::Image::Pointer	mitkTexture = dynamic_cast<mitk::Image*> (textureNode->GetData());
    mitk::Image::Pointer	mitkMask = dynamic_cast<mitk::Image*>(maskNode->GetData());
    // check if images are valid
    if ((!mitkTexture) || (!mitkMask) || (mitkTexture->IsInitialized() == false) || (mitkMask->IsInitialized() == false)){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    //-------------------------------------------------------------------------
    // path and image preparation
    QString strSavePath = m_WorkingDir;

    mitk::Image::Pointer mitkForecastImage = mitkTexture->Clone();
    ImageType::Pointer itkForecast;
    mitk::CastToItkImage(mitkForecastImage, itkForecast);
    itkForecast->FillBuffer(0);

    // transform the mitk texture
    ImageType::Pointer itkTexture;
    mitk::CastToItkImage(mitkTexture, itkTexture);
    BinaryImageType::Pointer itkMask;
    mitk::CastToItkImage(mitkMask, itkMask);

    //-------------------------------------------------------------------------
    // start
    int iWidth = mitkTexture->GetDimension(0);
    int iHeight = mitkTexture->GetDimension(1);
    int iSlice = mitkTexture->GetDimension(2);

    // file header: sz_data szInput szOutput
    // 	int szData = N_NODULE_D * N_NODULE_D * N_NODULE_D;
    int szInput = N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D;
    int szOutput = 1;
    int szAllData = 0;
    ImageType::SizeType image_size;
    image_size[0] = iWidth;
    image_size[1] = iHeight;
    image_size[2] = iSlice;

    // collect regions
    VVisualizeNodules vNoduleRegions;
    collectRegions(itkMask, vNoduleRegions);

    int szData = 0;
    szAllData = vNoduleRegions.size();
    for (int iVCenter = 0; iVCenter < szAllData; iVCenter++)
    {
        SVisualizeNodule region_infor = vNoduleRegions[iVCenter];
        ImageType::RegionType nodule_region = region_infor.nodule_region;
        ImageType::IndexType nodule_index = nodule_region.GetIndex();
        ImageType::SizeType nodule_size = nodule_region.GetSize();
#ifdef USE_NODULE_CENTERS
        int iRegionX = std::min((int)nodule_size[0], N_NODULE_D);
        int iRegionY = std::min((int)nodule_size[1], N_NODULE_D);
        int iRegionZ = std::min((int)nodule_size[2], N_NODULE_D);
        szData += iRegionX*iRegionY*iRegionZ;
#else
        szData += region_infor.nodule_indexes.size();
#endif
    }
    //     // collect some information
    //     itk::ImageRegionIterator<BinaryImageType>  iteratorMask(itkMask, itkMask->GetRequestedRegion());
    //     iteratorMask.GoToBegin();
    //     while (!iteratorMask.IsAtEnd()){
    //         BinaryImageType::PixelType val = iteratorMask.Get();
    //         if (val > 0){
    //             szAllData++;
    //         }
    //         ++iteratorMask;
    //     }

    MITK_INFO << "Forecast Data Count: " << szAllData;

    // texture 
    // open file writer pointers
    QFile qForecast(strSavePath + G_ANN_FORECAST_FILE_NAME);
    MITK_INFO << "Saving Forecast to " << strSavePath + G_ANN_FORECAST_FILE_NAME;
    if (qForecast.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
    {
        QTextStream tsForecast(&qForecast);
        tsForecast.setRealNumberNotation(QTextStream::FixedNotation);
        tsForecast.setRealNumberPrecision(PIXEL_PRECISION);

        // container to store indexes 
        VIndexes vIndexes;

        tsForecast << szData << "\t" << szInput << "\t" << szOutput << "\n";
        for (int iVCenter = 0; iVCenter < szAllData; iVCenter++)
        {
            strQDebug.sprintf("Outputting Forecast Nodule %d", iVCenter);
            MITK_INFO << strQDebug.toStdString();

            SVisualizeNodule region_infor = vNoduleRegions[iVCenter];
            ImageType::RegionType noduleRegion = region_infor.nodule_region;
            ImageType::IndexType noduleIndex = noduleRegion.GetIndex();
            ImageType::SizeType noduleSize = noduleRegion.GetSize();

#ifdef USE_NODULE_CENTERS
            ImageType::IndexType centerIndex;
            centerIndex[0] = noduleIndex[0] + noduleSize[0] / 2;
            centerIndex[1] = noduleIndex[1] + noduleSize[1] / 2;
            centerIndex[2] = noduleIndex[2] + noduleSize[2] / 2;

            int iRegionX = std::min((int)noduleSize[0], N_NODULE_D);
            int iRegionY = std::min((int)noduleSize[1], N_NODULE_D);
            int iRegionZ = std::min((int)noduleSize[2], N_NODULE_D);

            ImageType::SizeType regionSize;
            regionSize[0] = iRegionX;
            regionSize[1] = iRegionY;
            regionSize[2] = iRegionZ;

            ImageType::IndexType regionIndex;
            regionIndex[0] = centerIndex[0] - iRegionX / 2 - 1;
            regionIndex[1] = centerIndex[1] - iRegionY / 2 - 1;
            regionIndex[2] = centerIndex[2] - iRegionZ / 2 - 1;

            // start at index and look for neighbors at each point inside of index+size;
            ImageType::RegionType region;
            // region size
            region.SetSize(regionSize);
            // region start
            region.SetIndex(regionIndex);

            ImageType::SizeType radius;
            radius[0] = N_LOCAL_WINDOW_R;
            radius[1] = N_LOCAL_WINDOW_R;
            radius[2] = N_LOCAL_WINDOW_R;

            itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
            // one piece of training data
            while (!iterator.IsAtEnd())
            {
                int nTotalIn = 0;
                ImageType::IndexType center_index = iterator.GetIndex();
                vIndexes.push_back(center_index);

                for (unsigned int i = 0; i < iterator.Size(); i++)
                {
                    ImageType::IndexType local_index = iterator.GetIndex(i);

                    bool IsInBounds;
                    int neighborValue = iterator.GetPixel(i, IsInBounds);

                    if (!IsInBounds){
                        neighborValue = -(int)F_CT_WIDTH;
                    }
                    if (neighborValue > F_CT_WIDTH){
                        neighborValue = F_CT_WIDTH;
                    }
                    if (neighborValue < -F_CT_WIDTH){
                        neighborValue = -F_CT_WIDTH;
                    }
                    if (IsInBounds){
                        itkForecast->SetPixel(local_index, neighborValue);
                    }
                    double dwNormalized = ((double)(neighborValue + F_CT_WIDTH)) / ((double)(F_CT_WIDTH * 2));
                    tsForecast << dwNormalized/*neighborValue */ << "\t";
                }
                tsForecast << "\n";
                tsForecast << 0;
                tsForecast << "\n";

                ++iterator;
            }
#else
            int szIndexes = region_infor.nodule_indexes.size();
            for (int i = 0; i < szIndexes; i++)
            {
                ImageType::SizeType regionSize;
                regionSize[0] = 1;
                regionSize[1] = 1;
                regionSize[2] = 1;

                ImageType::IndexType regionIndex;
                regionIndex = region_infor.nodule_indexes[i];

                // start at index and look for neighbors at each point inside of index+size;
                ImageType::RegionType region;
                // region size
                region.SetSize(regionSize);
                // region start
                region.SetIndex(regionIndex);

                ImageType::SizeType radius;
                radius[0] = N_LOCAL_WINDOW_R;
                radius[1] = N_LOCAL_WINDOW_R;
                radius[2] = N_LOCAL_WINDOW_R;

                itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
                // one piece of training data
                while (!iterator.IsAtEnd())
                {
                    int nTotalIn = 0;
                    ImageType::IndexType center_index = iterator.GetIndex();
                    vIndexes.push_back(center_index);

                    for (unsigned int i = 0; i < iterator.Size(); i++)
                    {
                        ImageType::IndexType local_index = iterator.GetIndex(i);

                        bool IsInBounds;
                        int neighborValue = iterator.GetPixel(i, IsInBounds);

                        if (!IsInBounds){
                            neighborValue = -(int)F_CT_WIDTH;
                        }
                        if (neighborValue > F_CT_WIDTH){
                            neighborValue = F_CT_WIDTH;
                        }
                        if (neighborValue < -F_CT_WIDTH){
                            neighborValue = -F_CT_WIDTH;
                        }
                        if (IsInBounds){
                            itkForecast->SetPixel(local_index, neighborValue);
                        }
                        double dwNormalized = ((double)(neighborValue + F_CT_WIDTH)) / ((double)(F_CT_WIDTH * 2));
                        tsForecast << dwNormalized/*neighborValue */ << "\t";
                    }
                    tsForecast << "\n";
                    tsForecast << 0;
                    tsForecast << "\n";

                    ++iterator;
                }
            }
#endif
        }
        // #else
        //         // output params
        //         tsForecast << szAllData << "\t" << szInput << "\t" << szOutput << "\n";
        //         iteratorMask.GoToBegin();
        //         unsigned int nOutputPixelCount = 0;
        //         while (!iteratorMask.IsAtEnd())
        //         {
        // //             if (nOutputPixelCount % 100 == 0)
        // //             {
        // //                 strQDebug.sprintf("Outputting Nodule Pixel %d", nOutputPixelCount);
        // //                 MITK_INFO << strQDebug.toStdString();
        // //             }
        // // 
        // //             nOutputPixelCount++;
        //             BinaryImageType::PixelType val = iteratorMask.Get();
        //             if (val > 0)
        //             {
        //                 BinaryImageType::IndexType index = iteratorMask.GetIndex();
        // 
        //                 // remember output indexes
        //                 vIndexes.push_back(index);
        // 
        //                 ImageType::SizeType regionSize;
        //                 regionSize[0] = 1;
        //                 regionSize[1] = 1;
        //                 regionSize[2] = 1;
        // 
        //                 ImageType::IndexType regionIndex;
        //                 regionIndex[0] = index[0];
        //                 regionIndex[1] = index[1];
        //                 regionIndex[2] = index[2];
        // 
        //                 // start at index and look for neighbors at each point inside of index+size;
        //                 ImageType::RegionType region;
        //                 // region size
        //                 region.SetSize(regionSize);
        //                 // region start
        //                 region.SetIndex(regionIndex);
        // 
        //                 ImageType::SizeType radius;
        //                 radius[0] = N_LOCAL_WINDOW_R;
        //                 radius[1] = N_LOCAL_WINDOW_R;
        //                 radius[2] = N_LOCAL_WINDOW_R;
        // 
        //                 itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
        //                 // one piece of training data
        //                 while (!iterator.IsAtEnd())
        //                 {
        //                     int nTotalIn = 0;
        //                     for (unsigned int i = 0; i < iterator.Size(); i++)
        //                     {
        //                         ImageType::IndexType local_index = iterator.GetIndex(i);
        // 
        //                         bool IsInBounds;
        //                         int neighborValue = iterator.GetPixel(i, IsInBounds);
        // 
        //                         if (!IsInBounds)
        //                         {
        //                             neighborValue = -(int)F_CT_WIDTH;
        //                         }
        //                         if (!IsInBounds){
        //                             neighborValue = -(int)F_CT_WIDTH;
        //                         }
        //                         if (neighborValue > F_CT_WIDTH){
        //                             neighborValue = F_CT_WIDTH;
        //                         }
        //                         if (neighborValue < -F_CT_WIDTH){
        //                             neighborValue = -F_CT_WIDTH;
        //                         }
        //                         if (IsInBounds){
        //                             itkForecast->SetPixel(local_index, neighborValue);
        //                         }
        // 
        //                         double dwNormalized = ((double)(neighborValue + F_CT_WIDTH)) / ((double)(F_CT_WIDTH * 2));
        //                         tsForecast << dwNormalized/*neighborValue */ << "\t";
        // 
        //                         nTotalIn++;
        //                     }
        //                     tsForecast << "\n";
        //                     tsForecast << 0;
        //                     tsForecast << "\n";
        // 
        //                     if (nTotalIn != szInput)
        //                     {
        //                         strQDebug.sprintf("number of output training data not equal with sample size, Needed: %d, Actual: %d", szInput, nTotalIn);
        //                         MITK_WARN << strQDebug;
        //                     }
        //                     ++iterator;
        //                 }
        //             }
        //             ++iteratorMask;
        //         }
        // #endif //USE_NODULE_CENTERS

        // output remembered indexes
        SaveIndexes(strSavePath + G_ANN_INDEX_FILE_NAME, image_size, vNoduleRegions);
        SaveRegionParams(strSavePath + G_ANN_REGION_FILE_NAME, vNoduleRegions);
    }
    else
    {
        MITK_WARN << "Open " << strSavePath + G_ANN_FORECAST_FILE_NAME << " Failed...";
    }

}

////////////////////////////////////////////////////////////////////////////
//void QmitkRegionGrowingView::DoCrop()
//{
//	ImageType::SizeType imageSize = mCurrentRegion.GetSize();
//	ImageType::IndexType imageIdx = mCurrentRegion.GetIndex();
//	if (imageSize[0] == 0 || imageSize[1] == 0 || imageSize[2] == 0)
//	{
//		MITK_WARN << "No ROI Set, Select ROI first!!";
//		return;
//	}
//
//	mitk::Image::Pointer	pSel = getSelectedImage();
//	mitk::DataNode::Pointer pSelNode = getSelectedNode();
//	mitk::Image::Pointer	mitkImage = pSel->Clone();
//	ImageType::Pointer itkImage;
//	mitk::CastToItkImage(mitkImage, itkImage);
//
//	typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
//	FilterType::Pointer roiFilter = FilterType::New();
//	roiFilter->SetRegionOfInterest(mCurrentRegion);
//	roiFilter->SetInput(itkImage);
//	roiFilter->Update();
//
//	std::string str_nodename = pSelNode->GetName();
//	QString qstr_name(str_nodename.c_str());
//	if (qstr_name.length() > 64){
//		qstr_name = "image";
//	}
//
//	QString strFormat;
//	strFormat.sprintf("%s_crop_d_%d_%d", qstr_name.toStdString().c_str(), imageSize[0], imageSize[1], imageSize[2]);
//	mitk::Image::Pointer mitkNewColorImage = mitk::Image::New();
//	mitk::CastToMitkImage(roiFilter->GetOutput(), mitkNewColorImage);
//	ShowProcessedDataNode(mitkNewColorImage, strFormat.toStdString());
//}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::buildNormalizedGaussian(double***&dwRET, int nGaussianRX, int nGaussianRY, int nGaussianRZ, bool bSigmaAdjustable)
{
    const double	DW_SIGMA_STD = 4.0/*0.39894228040141*/;
    const double    DW_RADIUS_STD = 4.0;

    int nGaussianR = std::max(nGaussianRX, nGaussianRY);
    nGaussianR = std::max(nGaussianR, nGaussianRZ);

    double dwThisSigma = 0;
    if (bSigmaAdjustable)
        dwThisSigma = (nGaussianR / DW_RADIUS_STD)*DW_SIGMA_STD;
    else
        dwThisSigma = DW_SIGMA_STD;

    double dwThisSigma2 = dwThisSigma*dwThisSigma;

    //     MITK_INFO << "building gaussian with radius " << nGaussianR;
    //     int nGaussianD = 2 * nGaussianR + 1;
    int nGaussianDX = 2 * nGaussianRX + 1;
    int nGaussianDY = 2 * nGaussianRY + 1;
    int nGaussianDZ = 2 * nGaussianRZ + 1;

    double dwGaussianSigma = 0.0;
    double dwGaussianMax = 0.0;
    double dwGaussianMin = 1000000.0;
    for (int local_x_idx = 0; local_x_idx < nGaussianDX; local_x_idx++){
        for (int local_y_idx = 0; local_y_idx < nGaussianDY; local_y_idx++){
            for (int local_z_idx = 0; local_z_idx < nGaussianDZ; local_z_idx++){
                int local_x = (local_x_idx - nGaussianRX);
                int local_y = (local_y_idx - nGaussianRY);
                int local_z = (local_z_idx - nGaussianRZ);
                double dwThisVal = 1.0 / sqrt(2 * DW_MATH_PI* dwThisSigma2) * exp(-(local_x*local_x + local_y*local_y + local_z*local_z) / (2 * dwThisSigma2));
                dwRET[local_x_idx][local_y_idx][local_z_idx] = dwThisVal;
                if (dwThisVal > dwGaussianMax){
                    dwGaussianMax = dwThisVal;
                }
                if (dwThisVal < dwGaussianMin){
                    dwGaussianMin = dwThisVal;
                }

                dwGaussianSigma += dwThisVal;
            }
        }
    }

    // normalize 
    for (int local_x_idx = 0; local_x_idx < nGaussianDX; local_x_idx++){
        for (int local_y_idx = 0; local_y_idx < nGaussianDY; local_y_idx++){
            for (int local_z_idx = 0; local_z_idx < nGaussianDZ; local_z_idx++){
                //dwRET[local_x_idx][local_y_idx][local_z_idx] = (dwRET[local_x_idx][local_y_idx][local_z_idx] - dwGaussianMin)/(dwGaussianMax - dwGaussianMin);
                dwRET[local_x_idx][local_y_idx][local_z_idx] /= dwGaussianSigma;
            }
        }
    }
    // 	dwGaussianMax /= dwGaussianSigma;
    // 	return dwGaussianMax;
}

//////////////////////////////////////////////////////////////////////////
bool QmitkRegionGrowingView::normalizeRegionValues(VIndexValue& vIndexValues, VVisualizeNodules& vOutVisNodules)
{
    QString qSTRDBG;

    int nIndexCounter = 0;
    int szRegions = vOutVisNodules.size();

    double* pRegionLocalMaxIndexValue = new double[szRegions];
    double* pRegionLocalMinIndexValue = new double[szRegions];
    for (int i = 0; i < szRegions; i++){
        pRegionLocalMaxIndexValue[i] = -999999999;
        pRegionLocalMinIndexValue[i] = 999999999;
    }

    for (int i = 0; i < szRegions; i++)
    {
#ifdef USE_NODULE_CENTER
        ImageType::RegionType& region = vOutVisNodules[i].nodule_region;
        ImageType::IndexType nodule_index = region.GetIndex();
        ImageType::SizeType nodule_size = region.GetSize();

        int iRegionX = std::min((int)nodule_size[0], N_NODULE_D);
        int iRegionY = std::min((int)nodule_size[1], N_NODULE_D);
        int iRegionZ = std::min((int)nodule_size[2], N_NODULE_D);

        ImageType::SizeType out_nodule_size;
        out_nodule_size[0] = iRegionX;
        out_nodule_size[1] = iRegionY;
        out_nodule_size[2] = iRegionZ;
        vOutVisNodules[i].nodule_region.SetIndex(nodule_index);
        vOutVisNodules[i].nodule_region.SetSize(out_nodule_size);

        int szThisNodule = iRegionX * iRegionY * iRegionZ;
#else
        SVisualizeNodule& region_infor = vOutVisNodules[i];
        ImageType::RegionType& region = region_infor.nodule_region;
        int szThisNodule = region_infor.nodule_indexes.size();
        if (szThisNodule != region_infor.nodule_size)
        {
            qSTRDBG.sprintf("Nodule %d Size Info not Match should %d, actual %d", i, region_infor.nodule_size, szThisNodule);
            MITK_WARN << qSTRDBG.toStdString();
        }

        ImageType::IndexType nodule_index = region.GetIndex();
        ImageType::SizeType nodule_size = region.GetSize();
#endif

        for (int j = 0; j < szThisNodule; j++){
            if (vIndexValues[nIndexCounter + j] > pRegionLocalMaxIndexValue[i]){
                pRegionLocalMaxIndexValue[i] = vIndexValues[nIndexCounter + j];
            }
            if (vIndexValues[nIndexCounter + j] < pRegionLocalMinIndexValue[i]){
                pRegionLocalMinIndexValue[i] = vIndexValues[nIndexCounter + j];
            }
        }
        nIndexCounter += szThisNodule;
    }

    if (nIndexCounter != vIndexValues.size()){
        qSTRDBG.sprintf("Index Counter %d and Size of IndexValues %d not match", nIndexCounter, vIndexValues.size());
        MITK_WARN << qSTRDBG.toStdString();
    }

    double dwMaxInMax = -99999999;
    double dwMinInMin = 99999999;
    nIndexCounter = 0;

    for (int i = 0; i < szRegions; i++)
    {
        if (pRegionLocalMaxIndexValue[i] > dwMaxInMax)
            dwMaxInMax = pRegionLocalMaxIndexValue[i];
        if (pRegionLocalMinIndexValue[i] < dwMinInMin)
            dwMinInMin = pRegionLocalMinIndexValue[i];
    }

    qSTRDBG.sprintf("Normalizing Index Values with Max %f, Min %f", dwMaxInMax, dwMinInMin);
    MITK_INFO << qSTRDBG.toStdString();

    for (int i = 0; i < szRegions; i++)
    {
#ifdef USE_NODULE_CENTER
        ImageType::SizeType out_nodule_size = vOutVisNodules[i].nodule_region.GetSize();
        int szThisNodule = out_nodule_size[0] * out_nodule_size[1] * out_nodule_size[2];
#else
        SVisualizeNodule& region_infor = vOutVisNodules[i];
        ImageType::RegionType& region = region_infor.nodule_region;
        int valid_index_this_region = region_infor.nodule_size;
        int szThisNodule = valid_index_this_region;
#endif

        for (int j = 0; j < szThisNodule; j++)
        {
            double dwScaledValue = (vIndexValues[nIndexCounter + j] - dwMinInMin) / (dwMaxInMax - dwMinInMin);
            vOutVisNodules[i].nodule_index_values.push_back(dwScaledValue);
        }
        nIndexCounter += szThisNodule;
    }

    SAFE_DELETE(pRegionLocalMaxIndexValue);
    SAFE_DELETE(pRegionLocalMinIndexValue);

    MITK_INFO << "normalize VisNodules values: " << vOutVisNodules.size() << " Done";

    return true;
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::SaveRegionParams(QString strFile, VVisualizeNodules& regions)
{
    MITK_INFO << "Saving Regions to: " << strFile;

    QFile qValues(strFile);
    if (qValues.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
    {
        QTextStream tsValues(&qValues);
        int szRegions = regions.size();

        tsValues << szRegions << "\n";

        for (int i = 0; i < szRegions; i++)
        {
            SVisualizeNodule& thisRegionInfo = regions[i];
            ImageType::RegionType& thisRegion = thisRegionInfo.nodule_region;
            ImageType::IndexType index = thisRegion.GetIndex();
            ImageType::SizeType size = thisRegion.GetSize();
            tsValues << index[0] << "\t" << index[1] << "\t" << index[2] << "\t" << size[0] << "\t" << size[1] << "\t" << size[2] << "\t" << thisRegionInfo.nodule_size << "\n";
        }
        qValues.close();
    }
}

bool QmitkRegionGrowingView::LoadRegionParams(QString strFile, VVisualizeNodules& regions_out)
{
    MITK_INFO << "Loading Region Params: " << strFile;
    QString strDBG;

    QFile qValues(strFile);

    if (qValues.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QTextStream tsValues(&qValues);
        QString strLine;
        QStringList strList;

        ImageType::IndexType index;
        ImageType::SizeType size;
        int valid_index_count;
        SVisualizeNodule region_infor;

        // size of data set
        QString strSizeAll = tsValues.readLine();
        int szAll = strSizeAll.toInt();

        while (!tsValues.atEnd())
        {
            strLine = tsValues.readLine();
            strList = strLine.split(QRegExp("\\s+"));
            index[0] = strList[0].toInt();
            index[1] = strList[1].toInt();
            index[2] = strList[2].toInt();
            size[0] = strList[3].toInt();
            size[1] = strList[4].toInt();
            size[2] = strList[5].toInt();
            valid_index_count = strList[6].toInt();
            region_infor.nodule_region.SetSize(size);
            region_infor.nodule_region.SetIndex(index);
            region_infor.nodule_size = valid_index_count;
            regions_out.push_back(region_infor);

            //             strDBG.sprintf("Region Param Index(%d %d %d), Size(%d %d %d), NoduleSize %d",
            //                 index[0], index[1], index[2], size[0], size[1], size[2], valid_index_count);
            //             MITK_INFO << strDBG.toStdString();
        }

        qValues.close();

        if (szAll != regions_out.size()){
            MITK_WARN << "Region Size not Match ! " << "Want: " << szAll << ", Actual: " << regions_out.size();
        }

        return true;
    }
    return false;

}


#define ANN_SHOW_GAUSSIAN
//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoShowAnnResults()
{
    mitk::DataStorage::Pointer pDataStorage = GetDataStorage();
    mitk::DataNode::Pointer baseimagenode = pDataStorage->GetNamedNode("image_extract_bg_-2048");

    QString str_node_prefix = "ANNGaussianResults_";
    for (int i = 0; i < EANNRT_COUNT; i++)
    {

        mitk::DataNode::Pointer ann_reuslt_node = pDataStorage->GetNamedNode((str_node_prefix + STR_ANN_RESULTS_TYPE[i]).toStdString());
        if (ann_reuslt_node != nullptr)
            pDataStorage->Remove(ann_reuslt_node);
    }

    str_node_prefix = "ANNResults_";
    for (int i = 0; i < EANNRT_COUNT; i++)
    {

        mitk::DataNode::Pointer ann_reuslt_node = pDataStorage->GetNamedNode((str_node_prefix + STR_ANN_RESULTS_TYPE[i]).toStdString());
        if (ann_reuslt_node != nullptr)
            pDataStorage->Remove(ann_reuslt_node);
    }

    QString     strFile = m_WorkingDir;
    bool        bShowGaussian = m_Controls.cbShowGaussian->isChecked();
    bool        bShowGround = m_Controls.cbGroundMask->isChecked();
    QString     strShownResults = m_Controls.leShownPriorNum->text();
    int         nShownResults = strShownResults.toInt();

    m_Controls.treeWidget->clear();
    m_vvTotalNoduleIndex.clear();
    m_vvTotalNoduleIndex.resize(EANNRT_COUNT);

    if (m_Controls.cbShowVessel->isChecked())
        showOneAnnResult(strFile, EANNRT_VESSEL, baseimagenode, m_vvTotalNoduleIndex[EANNRT_VESSEL], bShowGaussian, bShowGround, nShownResults);
    if (m_Controls.cbShowGGO->isChecked())
        showOneAnnResult(strFile, EANNRT_GGO, baseimagenode, m_vvTotalNoduleIndex[EANNRT_GGO], bShowGaussian, bShowGround, nShownResults);
    if (m_Controls.cbShowWall->isChecked())
        showOneAnnResult(strFile, EANNRT_WALL, baseimagenode, m_vvTotalNoduleIndex[EANNRT_WALL], bShowGaussian, bShowGround, nShownResults);
    if (m_Controls.cbShowISO->isChecked())
        showOneAnnResult(strFile, EANNRT_ISO, baseimagenode, m_vvTotalNoduleIndex[EANNRT_ISO], bShowGaussian, bShowGround, nShownResults);
}

//////////////////////////////////////////////////////////////////////////
bool QmitkRegionGrowingView::LoadIndexValues(QString strFile, VVisualizeNodules& regions)
{
    MITK_INFO << "Loading Index Values: " << strFile;

    QFile qValues(strFile);

    VIndexValue vIndexValues;
    if (qValues.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QTextStream tsValues(&qValues);
        // skip the comment
        tsValues.readLine();
        while (!tsValues.atEnd())
        {
            double dwVal;
            dwVal = tsValues.readLine().toDouble();
            vIndexValues.push_back(dwVal);
        }
        qValues.close();
        normalizeRegionValues(vIndexValues, regions);
        return true;
    }

    return false;
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::SaveIndexes(QString strFile, ImageType::SizeType& image_size, VVisualizeNodules& regions)
{
    MITK_INFO << "Saving Indexes to: " << strFile;

    QFile qValues(strFile);
    if (qValues.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
    {
        QTextStream tsValues(&qValues);
        tsValues << image_size[0] << "\t" << image_size[1] << "\t" << image_size[2] << "\n";

        for (int i = 0; i < regions.size(); i++)
        {
            VIndexes& vIndexes = regions[i].nodule_indexes;
            int szIndexes = vIndexes.size();
            for (int j = 0; j < szIndexes; j++)
            {
                ImageType::IndexType& index = vIndexes[j];
                tsValues << index[0] << "\t" << index[1] << "\t" << index[2] << "\n";
            }
        }
        qValues.close();
    }
}

//////////////////////////////////////////////////////////////////////////
bool QmitkRegionGrowingView::LoadIndexes(QString strFile, ImageType::SizeType& image_size, VVisualizeNodules& regions)
{
    MITK_INFO << "Loading Indexes: " << strFile;
    QString qSTRDBG;

    QFile qValues(strFile);
    if (qValues.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QTextStream tsValues(&qValues);
        ImageType::IndexType image_index;

        QString		strLine = tsValues.readLine();
        QStringList strList = strLine.split(QRegExp("\\s+"));
        image_size[0] = strList[0].toInt();
        image_size[1] = strList[1].toInt();
        image_size[2] = strList[2].toInt();

        VIndexes vIndexes;
        while (!tsValues.atEnd())
        {
            strLine = tsValues.readLine();
            strList = strLine.split(QRegExp("\\s+"));
            image_index[0] = strList[0].toInt();
            image_index[1] = strList[1].toInt();
            image_index[2] = strList[2].toInt();
            vIndexes.push_back(image_index);
        }
        qValues.close();

        int nIndexCounter = 0;
        VIndexes::iterator index_begin = vIndexes.begin();
        for (int i = 0; i < regions.size(); i++)
        {
            SVisualizeNodule& region_infor = regions[i];
            ImageType::RegionType& region = region_infor.nodule_region;
            int valid_index_this_region = region_infor.nodule_size;

            region_infor.nodule_indexes.insert(region_infor.nodule_indexes.begin(),
                index_begin + nIndexCounter, index_begin + nIndexCounter + valid_index_this_region);
            nIndexCounter += valid_index_this_region;
        }
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////
//void QmitkRegionGrowingView::DoLoadIndexValues()
//{
//	QString strFile = getCurrentDataPath();
//	strFile += "/";
//	strFile += G_ANN_INDEX_VALUES_FILE_NAME;
//
//	VIndexValue vIndexValues;
//	// 	double dwMax, dwSigma;
//	LoadIndexValues(strFile, vIndexValues/*, dwMax, dwSigma*/);
//}

////////////////////////////////////////////////////////////////////////////
//void QmitkRegionGrowingView::DoSaveIndexes()
//{
//	QString strFile = getCurrentDataPath();
//	strFile += "/";
//	strFile += G_ANN_INDEX_FILE_NAME;
//
//	SaveIndexes(strFile);
//}

////////////////////////////////////////////////////////////////////////////
//void QmitkRegionGrowingView::DoLoadIndexes()
//{
//	QString strFile = getCurrentDataPath();
//	strFile += "/";
//	strFile += G_ANN_INDEX_FILE_NAME;
//
//	// 	VVIndexes vvIndexes;
//	// 	ImageType::SizeType image_size;
//	// 	LoadIndexes(strFile, image_size, vvIndexes);
//}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::NodeAdded(const mitk::DataNode* node)
{
    mitk::PropertyList* pProList = node->GetPropertyList();
    std::string strPath;
    pProList->GetStringProperty("path", strPath);
    QString qPath(strPath.c_str());
    if (qPath != "")
    {
        int	idx = qPath.lastIndexOf("/");
        qPath = qPath.mid(0, idx);
        qPath += "/";
        m_WorkingDir = qPath;
        m_Controls.lbWorkingDir->setText(m_WorkingDir);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoSelectWorkingDir()
{
    //     QString tmp = QFileDialog::getExistingDirectory(
    //         NULL, "Select Working Directory", "E:\\LIDC\\LIDC-IDRI", QFileDialog::ShowDirsOnly);
    // 
    //     if (tmp != ""){
    //         m_WorkingDir = tmp;
    //     }
    //     m_WorkingDir += "/";
    // 
    //     m_Controls.lbWorkingDir->setText(m_WorkingDir);

    QString tmp = getLastOpenPath();
    if (tmp != ""){
        m_WorkingDir = tmp;
    }
    m_WorkingDir += "/";

    QMessageBox msgbox;
    msgbox.setText(m_WorkingDir);
    msgbox.exec();

    m_Controls.lbWorkingDir->setText(m_WorkingDir);
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoThinning()
{
    mitk::Image::Pointer mitkImage = getSelectedImage();
    mitk::DataNode::Pointer mitkNode = getSelectedNode();

    if (NULL != mitkImage)
    {
        BinaryImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage->Clone(), itkImage);

        typedef itk::BinaryThinningImageFilter <BinaryImageType, BinaryImageType> BinaryThinningImageFilterType;
        BinaryThinningImageFilterType::Pointer binaryThinningImageFilter = BinaryThinningImageFilterType::New();
        binaryThinningImageFilter->SetInput(itkImage);
        PROFILE_FUNCTION_BEGIN;
        binaryThinningImageFilter->Update();
        PROFILE_FUNCTION_H_END(thinning);
        // Rescale the image so that it can be seen (the output is 0 and 1, we want 0 and 255)
        typedef itk::RescaleIntensityImageFilter< BinaryImageType, BinaryImageType > RescaleType;
        RescaleType::Pointer rescaler = RescaleType::New();
        rescaler->SetInput(binaryThinningImageFilter->GetOutput());
        rescaler->SetOutputMinimum(0);
        rescaler->SetOutputMaximum(255);
        rescaler->Update();

        QString strDBG;
        strDBG.sprintf("%s_thinning", mitkNode->GetName().c_str());
        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(rescaler->GetOutput(), mitkResult);
        ShowProcessedDataNode(mitkResult, strDBG.toStdString(), false, mitkNode);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoOSTUThreshold()
{
    // 	QString strParam1 = m_Controls.leLowThreshold->text();
    // 	QString strParam2 = m_Controls.leHighThreshold->text();
    // 	float	param1 = strParam1.toInt();
    // 	float	param2 = strParam2.toInt();
    mitk::Image*	pSelected = getSelectedImage();
    mitk::DataNode* pSelNode = getSelectedNode();

    mitk::Image::Pointer mitkNewImage = pSelected->Clone();
    DoubleImageType::Pointer itkImage  /*= ImageType::New()*/;
    mitk::CastToItkImage(mitkNewImage, itkImage);

    typedef itk::OtsuThresholdImageFilter< DoubleImageType, BinaryImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(itkImage);
    filter->Update();

    MITK_INFO << "OTSU Threshold: " << filter->GetThreshold();

    ////  cast image to binary
    //typedef itk::CastImageFilter< DoubleImageType, BinaryImageType > ImageCastFilter;
    //ImageCastFilter::Pointer castFilter = ImageCastFilter::New();
    //castFilter->SetInput(filter->GetOutput());
    //castFilter->Update();

    // show results
    mitk::Image::Pointer pResult = mitk::Image::New();
    mitk::CastToMitkImage(filter->GetOutput(), pResult);
    QString qNodeName;
    qNodeName.sprintf("%s_otsu", pSelNode->GetName().c_str());
    ShowProcessedDataNode(pResult, qNodeName.toStdString(), true, pSelNode);

}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoKMeans()
{
    // validate data 
    mitk::DataNode::Pointer mitkNode = m_Controls.cbStatisticImage->GetSelectedNode();
    mitk::DataNode::Pointer featureNode = m_Controls.cbStatisticFeatureImage->GetSelectedNode();

    if (NULL == mitkNode || NULL == featureNode){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    mitk::Image::Pointer	mitkImage = dynamic_cast<mitk::Image*> (mitkNode->GetData());
    mitk::Image::Pointer	mitkFeatureImage = dynamic_cast<mitk::Image*>(featureNode->GetData());
    // check if images are valid
    if ((!mitkImage) || (!mitkFeatureImage) || (mitkImage->IsInitialized() == false) || (mitkFeatureImage->IsInitialized() == false)){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    // image conversion
    BinaryImageType::Pointer itkImage;
    mitk::CastToItkImage(mitkImage, itkImage);

    ImageType::Pointer itkFeatureImage;
    mitk::CastToItkImage(mitkFeatureImage, itkFeatureImage);

    // do 3d processing
    // convert the image to a label image using connect
    typedef itk::ConnectedComponentImageFilter <BinaryImageType, ImageType> ConnectedComponentImageFilterType_3D;
    ConnectedComponentImageFilterType_3D::Pointer connectedComponentImageFilter_3D = ConnectedComponentImageFilterType_3D::New();
    //     connectedComponentImageFilter_3D->FullyConnectedOn();
    connectedComponentImageFilter_3D->SetInput(itkImage);
    connectedComponentImageFilter_3D->Update();

    // statistical analysis
    typedef itk::LabelImageToStatisticsLabelMapFilter<ImageType, ImageType> StatisticsLabelMapFilter_3D;
    StatisticsLabelMapFilter_3D::Pointer statistics_filter_3d = StatisticsLabelMapFilter_3D::New();
    statistics_filter_3d->SetInput(connectedComponentImageFilter_3D->GetOutput());
    statistics_filter_3d->SetFeatureImage(itkFeatureImage);
    PROFILE_FUNCTION_BEGIN;
    statistics_filter_3d->Update();
    PROFILE_FUNCTION_H_END(statistics_filter_3d);

    // shape analysis
    typedef itk::LabelImageToShapeLabelMapFilter<ImageType> ShapeLabelMapFilter_3D;
    ShapeLabelMapFilter_3D::Pointer shape_filter_3d = ShapeLabelMapFilter_3D::New();
    shape_filter_3d->SetInput(connectedComponentImageFilter_3D->GetOutput());
    shape_filter_3d->Update();

    // analysis according to shape & statistic properties
    typedef itk::StatisticsLabelObject<TPixelType, VDIM>StatisticsLabelObjectType;
    typedef itk::LabelMap< StatisticsLabelObjectType >	StatisticsLabelMapType;
    StatisticsLabelMapType::Pointer statistics_labelMap = statistics_filter_3d->GetOutput();

    typedef itk::ShapeLabelObject<TPixelType, VDIM>	ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType >	ShapeLabelMapType;
    ShapeLabelMapType::Pointer shape_labelMap = shape_filter_3d->GetOutput();

    // roundness, nPixels, nSlices, mean, median, variance, Kurtosis, SKEWNESS, 
    typedef itk::Vector< double, 8 > MeasurementVectorType;
    typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
    SampleType::Pointer sample = SampleType::New();
    MeasurementVectorType mv;

    int szLabelMap0 = statistics_labelMap->GetNumberOfLabelObjects();
    int szLabelMap = shape_labelMap->GetNumberOfLabelObjects();

    MITK_INFO << "There are " << szLabelMap << " Label Maps";

    for (unsigned int n = 0; n < szLabelMap; ++n)
    {
        ShapeLabelObjectType*		labelObject = shape_labelMap->GetNthLabelObject(n);
        ImageType::RegionType bound = labelObject->GetBoundingBox();
        ImageType::SizeType bound_size = bound.GetSize();
        StatisticsLabelObjectType* statisticObject = statistics_labelMap->GetNthLabelObject(n);

        mv[0] = labelObject->GetRoundness();
        mv[1] = labelObject->GetNumberOfPixels();
        mv[2] = bound_size[2];
        mv[3] = statisticObject->GetMean();
        mv[4] = statisticObject->GetMedian();
        mv[5] = statisticObject->GetVariance();
        mv[6] = statisticObject->GetKurtosis();
        mv[7] = statisticObject->GetSkewness();
        sample->PushBack(mv);
    }

    typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
    TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

    treeGenerator->SetSample(sample);
    treeGenerator->SetBucketSize(16);
    treeGenerator->Update();

    typedef TreeGeneratorType::KdTreeType TreeType;
    typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
    EstimatorType::Pointer estimator = EstimatorType::New();

    EstimatorType::ParametersType initialMeans(16);
    initialMeans[0] = 0.9; // Cluster 1, mean[0]
    initialMeans[1] = 300.0; // Cluster 1, mean[1]
    initialMeans[2] = 20.0; // Cluster 1, mean[2]
    initialMeans[3] = 0.0; // Cluster 2, mean[0]
    initialMeans[4] = 0.0; // Cluster 2, mean[1]
    initialMeans[5] = 0.0; // Cluster 2, mean[2]
    initialMeans[6] = 0.0; // Cluster 2, mean[2]
    initialMeans[7] = 0.0; // Cluster 2, mean[2]

    initialMeans[8] = 0.5; // Cluster 1, mean[0]
    initialMeans[9] = 16.0; // Cluster 1, mean[1]
    initialMeans[10] = 1.0; // Cluster 1, mean[2]
    initialMeans[11] = 0.0; // Cluster 2, mean[1]
    initialMeans[12] = 0.0; // Cluster 2, mean[2]
    initialMeans[13] = 0.0; // Cluster 2, mean[2]
    initialMeans[14] = 0.0; // Cluster 2, mean[2]
    initialMeans[15] = 0.0; // Cluster 2, mean[2]

    estimator->SetParameters(initialMeans);
    estimator->SetKdTree(treeGenerator->GetOutput());
    estimator->SetMaximumIteration(200);
    estimator->SetCentroidPositionChangesThreshold(0.0);
    estimator->StartOptimization();

    EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();

    for (unsigned int i = 0; i < 16; i += 2)
    {
        std::cout << "cluster[" << i << "] " << std::endl;
        std::cout << "    estimated mean : " << estimatedMeans[i] << " , " << estimatedMeans[i + 1] << std::endl;
    }

    typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType > MembershipFunctionType;
    typedef MembershipFunctionType::Pointer MembershipFunctionPointer;

    typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
    DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();

    typedef itk::Statistics::SampleClassifierFilter< SampleType > ClassifierType;
    ClassifierType::Pointer classifier = ClassifierType::New();
    classifier->SetDecisionRule(decisionRule);
    classifier->SetInput(sample);
    classifier->SetNumberOfClasses(2);

    typedef ClassifierType::ClassLabelVectorObjectType               ClassLabelVectorObjectType;
    typedef ClassifierType::ClassLabelVectorType                     ClassLabelVectorType;
    typedef ClassifierType::MembershipFunctionVectorObjectType       MembershipFunctionVectorObjectType;
    typedef ClassifierType::MembershipFunctionVectorType             MembershipFunctionVectorType;

    ClassLabelVectorObjectType::Pointer  classLabelsObject = ClassLabelVectorObjectType::New();
    classifier->SetClassLabels(classLabelsObject);

    ClassLabelVectorType &  classLabelsVector = classLabelsObject->Get();
    classLabelsVector.push_back(100);
    classLabelsVector.push_back(200);

    MembershipFunctionVectorObjectType::Pointer membershipFunctionsObject = MembershipFunctionVectorObjectType::New();
    classifier->SetMembershipFunctions(membershipFunctionsObject);

    MembershipFunctionVectorType &  membershipFunctionsVector = membershipFunctionsObject->Get();
    MembershipFunctionType::CentroidType origin(sample->GetMeasurementVectorSize());
    int index = 0;
    for (unsigned int i = 0; i < 2; i++)
    {
        MembershipFunctionPointer membershipFunction = MembershipFunctionType::New();
        for (unsigned int j = 0; j < sample->GetMeasurementVectorSize(); j++)
        {
            origin[j] = estimatedMeans[index++];
        }
        membershipFunction->SetCentroid(origin);
        membershipFunctionsVector.push_back(membershipFunction.GetPointer());
    }
    classifier->Update();

    const ClassifierType::MembershipSampleType* membershipSample = classifier->GetOutput();
    ClassifierType::MembershipSampleType::ConstIterator iter = membershipSample->Begin();

    while (iter != membershipSample->End())
    {
        MITK_INFO << "measurement vector = " << iter.GetMeasurementVector()
            << "class label = " << iter.GetClassLabel();
        ++iter;
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoSelectLabel()
{
    mitk::DataNode::Pointer featureNode = m_Controls.cbStatisticFeatureImage->GetSelectedNode();
    mitk::Image::Pointer	mitkFeatureImage = dynamic_cast<mitk::Image*>(featureNode->GetData());
    mitk::BaseGeometry::Pointer geometry = mitkFeatureImage->GetGeometry();

    QString strDBG;
    mitk::Point3D crossPositionInWorldCoordinates = m_IRenderWindowPart->GetSelectedPosition();
    strDBG.sprintf("Mouse Corr (%f %f %f)", crossPositionInWorldCoordinates[0], crossPositionInWorldCoordinates[1], crossPositionInWorldCoordinates[2]);
    MITK_INFO << strDBG.toStdString();

    int szStoredLabels = m_vLabeledRegions.size();
    for (int i = 1; i < szStoredLabels; i++)
    {
        ImageType::RegionType& thisRegion = m_vLabeledRegions[i];
        ImageType::IndexType index = thisRegion.GetIndex();
        ImageType::SizeType size = thisRegion.GetSize();

        mitk::Point3D cornerPoint1InWorldCoordinates;
        geometry->IndexToWorld(index, cornerPoint1InWorldCoordinates);

        if (crossPositionInWorldCoordinates[0] >= cornerPoint1InWorldCoordinates[0] && crossPositionInWorldCoordinates[0] <= cornerPoint1InWorldCoordinates[0] + size[0] &&
            crossPositionInWorldCoordinates[1] >= cornerPoint1InWorldCoordinates[1] && crossPositionInWorldCoordinates[1] <= cornerPoint1InWorldCoordinates[1] + size[1] &&
            crossPositionInWorldCoordinates[2] >= cornerPoint1InWorldCoordinates[2] && crossPositionInWorldCoordinates[2] <= cornerPoint1InWorldCoordinates[2] + size[2])
        {
            strDBG.sprintf("%d", i);
            QMessageBox msgBox;
            msgBox.setText(strDBG);
            msgBox.exec();
        }

    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoMedian()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        FloatImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

#if 0
        typedef itk::MedianImageFilter< FloatImageType, FloatImageType> MedianFilterType;
#else
        typedef itk::CudaMedianImageFilter< FloatImageType, FloatImageType> MedianFilterType;
#endif
        MedianFilterType::Pointer medianFilter = MedianFilterType::New();
        MedianFilterType::InputSizeType size;
        int param1 = m_Controls.spinBox_1->value();
        size.Fill(param1);
        medianFilter->SetRadius(size);
        medianFilter->SetInput(itkImage);
        medianFilter->Update();

        QString strNode;
        strNode.sprintf("%s_median_%d", mitkNode->GetName().c_str(), (int)param1);
        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(medianFilter->GetOutput(), mitkResult);
        ShowProcessedDataNode(mitkResult, strNode.toStdString(), false, mitkNode);
    }
}
//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoDilation()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        FloatImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        // dilation to include the contour
        typedef itk::BinaryBallStructuringElement<unsigned char, VDIM>              BallType;

#if 0
        typedef itk::GrayscaleDilateImageFilter<FloatImageType, FloatImageType, BallType>    DilationFilterType;
#else
        typedef itk::CudaGrayscaleDilateImageFilter<FloatImageType, FloatImageType, BallType>    DilationFilterType;
#endif

        BallType binaryBall;
        int param2 = m_Controls.spinBox_2->value();
        binaryBall.SetRadius(param2);
        binaryBall.CreateStructuringElement();

        DilationFilterType::Pointer dilationFilter = DilationFilterType::New();
        dilationFilter->SetInput(itkImage);
        dilationFilter->SetKernel(binaryBall);
        dilationFilter->Update();

        QString strNode;
        strNode.sprintf("%s_dilate_%d", mitkNode->GetName().c_str(), (int)param2);
        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(dilationFilter->GetOutput(), mitkResult);
        ShowProcessedDataNode(mitkResult, strNode.toStdString(), false, mitkNode);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoFillBorder()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        FloatImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        // dilation to include the contour
        typedef itk::FlatStructuringElement<VDIM> StructuringElementType;
        StructuringElementType::RadiusType elementRadius;
        int radius = m_Controls.spinBox_2->value();
        elementRadius.Fill(radius);
        StructuringElementType structuringElement = StructuringElementType::Box(elementRadius);

        typedef itk::GrayscaleDilateImageFilter<FloatImageType, FloatImageType, StructuringElementType>    DilationFilterType;
        DilationFilterType::Pointer dilationFilter = DilationFilterType::New();
        dilationFilter->SetInput(itkImage);
        dilationFilter->SetKernel(structuringElement);
        dilationFilter->Update();

        QString strNode;
        strNode.sprintf("%s_dilate_%d", mitkNode->GetName().c_str(), (int)radius);
        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(dilationFilter->GetOutput(), mitkResult);
        ShowProcessedDataNode(mitkResult, strNode.toStdString(), false, mitkNode);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoTransferToMask()
{
    const unsigned char BG_TOLARANCE = 1;

    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        RGBImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        QString qstr_pixelR = m_Controls.leBGPixelR->text();
        QString qstr_pixelG = m_Controls.leBGPixelG->text();
        QString qstr_pixelB = m_Controls.leBGPixelB->text();

        unsigned char BG_R = qstr_pixelR.toUInt();
        unsigned char BG_G = qstr_pixelG.toUInt();
        unsigned char BG_B = qstr_pixelB.toUInt();

        BinaryImageType::Pointer itkBinaryImage = BinaryImageType::New();
        itkBinaryImage->SetRegions(itkImage->GetLargestPossibleRegion());
        itkBinaryImage->CopyInformation(itkImage);
        itkBinaryImage->Allocate();
        itkBinaryImage->FillBuffer(0);


        itk::ImageRegionIterator<RGBImageType> imageITKIterator(itkImage, itkImage->GetLargestPossibleRegion());
        int nCounter = 0;
        while (!imageITKIterator.IsAtEnd())
        {
            RGBImageType::PixelType val = imageITKIterator.Get();
            RGBImageType::IndexType index = imageITKIterator.GetIndex();

            unsigned char delta_red = abs(val[0] - BG_R);
            unsigned char delta_green = abs(val[1] - BG_G);
            unsigned char delta_blue = abs(val[2] - BG_B);

            // bg color
            if (!(delta_red <= BG_TOLARANCE && delta_green <= BG_TOLARANCE && delta_blue <= BG_TOLARANCE))
                itkBinaryImage->SetPixel(index, 1);
            ++imageITKIterator;
        }

        typedef itk::MultiplyImageFilter<BinaryImageType, BinaryImageType> MultiplyFilter_3D;
        MultiplyFilter_3D::Pointer mulfilter_3d = MultiplyFilter_3D::New();
        mulfilter_3d->SetInput(itkBinaryImage);
        mulfilter_3d->SetConstant(255);
        mulfilter_3d->Update();

        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(mulfilter_3d->GetOutput(), mitkResult);

        qStrDBG.sprintf("%s", "mask_rgb2binary");
        ShowProcessedDataNode(mitkResult, qStrDBG.toStdString(), false, mitkNode);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoChangeBGColor()
{
    bool bChecked1 = m_Controls.rbBGPixelRed->isChecked();
    bool bChecked2 = m_Controls.rbBGPixelBlack->isChecked();

    if (bChecked1)
    {
        m_Controls.leBGPixelR->setText("139");
        m_Controls.leBGPixelG->setText("35");
        m_Controls.leBGPixelB->setText("35");
    }
    else if (bChecked2)
    {
        m_Controls.leBGPixelR->setText("0");
        m_Controls.leBGPixelG->setText("0");
        m_Controls.leBGPixelB->setText("0");
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::showOneAnnResult(QString basepath, EANNResultsType result_type, mitk::DataNode::Pointer baseNode, VecInt& vOutTotalNoduleIndex,
    bool bOutputGaussian /*= false*/, bool bOutputGGroundtruth /*= false*/, int nShownResults /*= 4*/)
{
    QString strDBG;

    mitk::BaseData*	data = baseNode->GetData();
    mitk::Image::Pointer mitkImage = NULL;
    if (data)
        mitkImage = dynamic_cast<mitk::Image*>(data);

    mitk::Image::Pointer mitkNewImage1, mitkNewImage2, mitkNewImage3;
    if (NULL != mitkImage){
        mitkNewImage1 = mitkImage->Clone();
        mitkNewImage2 = mitkImage->Clone();
        mitkNewImage3 = mitkImage->Clone();
    }
    else{
        MITK_ERROR << "Select Base Evaluation Image";
        return;
    }
    MITK_INFO << "Image Creation Done!";

    VVisualizeNodules vVisualizeNodules;
    ImageType::SizeType image_size;

    // region params should be loaded first
    bool bLoad3 = LoadRegionParams(basepath + G_ANN_REGION_FILE_NAME, vVisualizeNodules);
    bool bLoad2 = LoadIndexes(basepath + G_ANN_INDEX_FILE_NAME, image_size, vVisualizeNodules);
    bool bLoad1 = LoadIndexValues(basepath + G_ANN_INDEX_VALUES_PREFIX + STR_ANN_RESULTS_TYPE[result_type] + ".txt", vVisualizeNodules);

    // check data
    if (!bLoad1){
        MITK_WARN << "Load Index Values Failed";
        return;
    }
    if (!bLoad2){
        MITK_WARN << "Load Indexes Failed";
        return;
    }
    if (!bLoad3){
        MITK_WARN << "Load Regions Failed";
        return;
    }

    // work image and copy image
    DoubleImageType::Pointer itkImage;
    DoubleImageType::Pointer itkDupImage;
    DoubleImageType::Pointer itkDupImage2;
    DoubleImageType::IndexType itkImageIndex;

    mitk::CastToItkImage(mitkNewImage1, itkImage);
    itkImage->FillBuffer(0);

    mitk::CastToItkImage(mitkNewImage2, itkDupImage);
    itkDupImage->FillBuffer(0);

    mitk::CastToItkImage(mitkNewImage3, itkDupImage2);
    itkDupImage2->FillBuffer(0);

    for_each(vVisualizeNodules.begin(), vVisualizeNodules.end(),
        [&itkImageIndex, &itkImage/*, &dwMaxIndexValue*/](SVisualizeNodule& this_nodule)
    {
        int szThisIndex = this_nodule.nodule_indexes.size();
        for (int i = 0; i < szThisIndex; i++)
        {
            DoubleImageType::IndexType& thisIndex = this_nodule.nodule_indexes[i];

            itkImageIndex[0] = thisIndex[0];
            itkImageIndex[1] = thisIndex[1];
            itkImageIndex[2] = thisIndex[2];
            // normalize and set 
            itkImage->SetPixel(itkImageIndex, this_nodule.nodule_index_values[i]);
        }
    });


    int szNodules = vVisualizeNodules.size();
    double* dwNoduleScore = new double[szNodules];
    if (bOutputGaussian)
    {
        // process
        for (int iNodule = 0; iNodule < szNodules; iNodule++)
        {
            dwNoduleScore[iNodule] = 0.0;
            SVisualizeNodule& this_nodule = vVisualizeNodules[iNodule];
            VIndexValue& this_index_values = this_nodule.nodule_index_values;
            VIndexes& this_indexes = this_nodule.nodule_indexes;

            ImageType::RegionType& this_reigon = this_nodule.nodule_region;
            ImageType::SizeType nodule_size = this_reigon.GetSize();
            ImageType::IndexType nodule_index = this_reigon.GetIndex();

            int iRegionX = std::min((int)nodule_size[0], N_NODULE_D);
            int iRegionY = std::min((int)nodule_size[1], N_NODULE_D);
            int iRegionZ = std::min((int)nodule_size[2], N_NODULE_D);

            int iRegionXR = std::ceil((iRegionX + 0.5) / 2);
            int iRegionYR = std::ceil((iRegionY + 0.5) / 2);
            int iRegionZR = std::ceil((iRegionZ + 0.5) / 2);

            //             int xCoorIndexSigma = 0;
            //             int yCoorIndexSigma = 0;
            //             int zCoorIndexSigma = 0;
            // 
            //             int nValidSigmaCount = 0;
            //             for (int i = 0; i < this_indexes.size(); i++)
            //             {
            //                 ImageType::IndexType& index_coord = this_indexes[i];
            //                 double index_value = this_index_values[i];
            //                 if (index_value > 0.4)
            //                 {
            //                     xCoorIndexSigma += index_coord[0];
            //                     yCoorIndexSigma += index_coord[1];
            //                     zCoorIndexSigma += index_coord[2];
            //                     nValidSigmaCount++;
            //                 }
            //             }
            // 
            //             double max_pixel_value = -9999999;
            //             int max_pixel_value_index = -1;
            //             for (unsigned int iIdxValue = 0; iIdxValue < this_index_values.size(); iIdxValue++)
            //             {
            //                 double dwThisValue = this_index_values[iIdxValue];
            //                 if (dwThisValue > max_pixel_value)
            //                 {
            //                     max_pixel_value = dwThisValue;
            //                     max_pixel_value_index = iIdxValue;
            //                 }
            //             }
            DoubleImageType::IndexType center_index;
            //             if (max_pixel_value_index >= 0)
            //             {
            //                 center_index[0] = this_indexes[max_pixel_value_index][0];
            //                 center_index[1] = this_indexes[max_pixel_value_index][1];
            //                 center_index[2] = this_indexes[max_pixel_value_index][2];
            //             }
            //             else
            //             {
            //                 center_index[0] = nodule_index[0] + nodule_size[0] / 2;
            //                 center_index[1] = nodule_index[1] + nodule_size[1] / 2;
            //                 center_index[2] = nodule_index[2] + nodule_size[2] / 2;
            //             }
            //             center_index[0] = xCoorIndexSigma / this_nodule.nodule_indexes.size();
            //             center_index[1] = yCoorIndexSigma / this_nodule.nodule_indexes.size();
            //             center_index[2] = zCoorIndexSigma / this_nodule.nodule_indexes.size();
            center_index[0] = nodule_index[0] + nodule_size[0] / 2;
            center_index[1] = nodule_index[1] + nodule_size[1] / 2;
            center_index[2] = nodule_index[2] + nodule_size[2] / 2;

            DoubleImageType::SizeType gaussian_size;
            gaussian_size[0] = iRegionX;
            gaussian_size[1] = iRegionY;
            gaussian_size[2] = iRegionZ;

            DoubleImageType::IndexType gaussian_index;
            gaussian_index[0] = center_index[0] - iRegionXR;
            gaussian_index[1] = center_index[1] - iRegionYR;
            gaussian_index[2] = center_index[2] - iRegionZR;

            if (gaussian_index[0] < 0)
                gaussian_index[0] = 0;
            if (gaussian_index[1] < 0)
                gaussian_index[1] = 0;
            if (gaussian_index[2] < 0)
                gaussian_index[2] = 0;

            // build 3D gaussian kernel
            //             int iGaussianR = std::max(iRegionXR, iRegionYR);
            //             iGaussianR = std::max(iGaussianR, iRegionZR);
            //        int iGaussianR = N_LOCAL_WINDOW_R;
            //int iGaussianD = iGaussianR * 2 + 1;
            int iGaussianDX = iRegionXR * 2 + 1;
            int iGaussianDY = iRegionYR * 2 + 1;
            int iGaussianDZ = iRegionZR * 2 + 1;

            double ***dwGaussian = (double***)malloc(iGaussianDX * sizeof(double**));
            for (int i = 0; i < iGaussianDX; i++){
                dwGaussian[i] = (double**)malloc(iGaussianDY * sizeof(double*));
                for (int j = 0; j < iGaussianDY; j++){
                    dwGaussian[i][j] = (double*)malloc(iGaussianDZ * sizeof(double));
                }
            }

            bool bAdjustSigma = nShownResults > 5 ? true : false;
            buildNormalizedGaussian(dwGaussian, iRegionXR, iRegionYR, iRegionZR, bAdjustSigma);

            //             strDBG.sprintf("Operating on region %d with R:%d, Center:(%d %d %d), RegionCenter:%s",
            //                 iNodule, iGaussianR, center_index[0], center_index[1], center_index[2], max_pixel_value_index >= 0 ? "NO" : "YES");
            //             MITK_INFO << strDBG.toStdString();
#if 1
            assert(this_index_values.size() == this_indexes.size());
            for (int i = 0; i < this_indexes.size(); i++)
            {
                ImageType::IndexType& local_this_index = this_indexes[i];
                int local_x = local_this_index[0] - center_index[0] + /*iGaussianR*/iRegionXR;
                int local_y = local_this_index[1] - center_index[1] + /*iGaussianR*/iRegionYR;
                int local_z = local_this_index[2] - center_index[2] + /*iGaussianR*/iRegionZR;

                if (local_x < 0)
                    local_x = 0;
                if (local_y < 0)
                    local_y = 0;
                if (local_z < 0)
                    local_z = 0;
                if (local_x >= iGaussianDX)
                    local_x = iGaussianDX - 1;
                if (local_y >= iGaussianDY)
                    local_y = iGaussianDY - 1;
                if (local_z >= iGaussianDZ)
                    local_z = iGaussianDZ - 1;
                //                 strDBG.sprintf("(%d %d %d)", local_x, local_y, local_z);
                //                 MITK_INFO << strDBG.toStdString();

                double dwThisGaussian = dwGaussian[local_x][local_y][local_z];
                double dwThisValue = this_index_values[i];

                itkDupImage2->SetPixel(local_this_index, dwThisGaussian);

                dwNoduleScore[iNodule] += dwThisGaussian * dwThisValue;
            }

#else

            // start at index and look for neighbors at each point inside of index+size;
            DoubleImageType::RegionType process_region;
            // region size
            process_region.SetSize(gaussian_size);
            // region start
            process_region.SetIndex(gaussian_index);

            DoubleImageType::SizeType gaussian_radius;
            //         gaussian_radius[0] = iRegionXR;
            //         gaussian_radius[1] = iRegionYR;
            //         gaussian_radius[2] = iRegionZR;
            gaussian_radius[0] = iGaussianR;
            gaussian_radius[1] = iGaussianR;
            gaussian_radius[2] = iGaussianR;

            itk::NeighborhoodIterator<DoubleImageType> iterator(gaussian_radius, itkImage, process_region/*process_region*/);
            // convolve
            while (!iterator.IsAtEnd())
            {
                DoubleImageType::PixelType dwPixelValue = 0.0;
                DoubleImageType::IndexType local_center_index = iterator.GetIndex();

                for (unsigned int iITER = 0; iITER < iterator.Size(); iITER++)
                {
                    DoubleImageType::IndexType current_index = iterator.GetIndex(iITER);
                    bool IsInBounds = false;
                    DoubleImageType::PixelType current_value = iterator.GetPixel(iITER, IsInBounds);
                    if (IsInBounds)
                    {
                        int local_x = current_index[0] - local_center_index[0];
                        int local_y = current_index[1] - local_center_index[1];
                        int local_z = current_index[2] - local_center_index[2];
                        int local_index_x = local_x + iGaussianR;
                        int local_index_y = local_y + iGaussianR;
                        int local_index_z = local_z + iGaussianR;

                        //                     strDBG.sprintf("(%d, %d, %d), GaussianD: %d", iRegionX, iRegionY, iRegionZ, iGaussianD);
                        //                     MITK_INFO << strDBG.toStdString();
                        double dwThisGaussian = dwGaussian[local_index_x][local_index_y][local_index_z];
                        dwPixelValue += current_value * dwThisGaussian;
                        //                     strDBG.sprintf("current value (%f, %f)", current_value, dwThisGaussian);
                        //                     MITK_WARN << strDBG.toStdString();
                    }
                    else
                    {
                        dwPixelValue += 0.0;
                    }
                }
                strDBG.sprintf("final value %f", dwPixelValue);
                MITK_WARN << strDBG.toStdString();
                itkDupImage->SetPixel(local_center_index, dwPixelValue);
                dwNoduleScore[iNodule] += dwPixelValue;
                ++iterator;
            }

            for (int i = 0; i < iGaussianD; i++){
                for (int j = 0; j < iGaussianD; j++){
                    SAFE_DELETE(dwGaussian[i][j]);
                }
                SAFE_DELETE(dwGaussian[i]);
            }
            SAFE_DELETE(dwGaussian);
#endif
        }
    }

    // rescale for better vision
    typedef itk::MultiplyImageFilter<DoubleImageType, ImageType> MultiplyCastType;
    typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > ROIFilter;

    // original result image
    MultiplyCastType::Pointer multiplyFilter = MultiplyCastType::New();
    multiplyFilter->SetInput(itkImage);
    multiplyFilter->SetConstant(255);
    multiplyFilter->Update();

    mitk::Image::Pointer mitkResults;
    mitk::CastToMitkImage(multiplyFilter->GetOutput(), mitkResults);

    QString strFormat;
    strFormat.sprintf("%s", "ANNResults_");
    strFormat += STR_ANN_RESULTS_TYPE[result_type];
    // allocate a new node and show
    mitk::DataNode::Pointer newNode = mitk::DataNode::New();
    newNode->SetData(mitkResults);
    newNode->SetProperty("name", mitk::StringProperty::New(strFormat.toStdString()));
    newNode->SetProperty("visible", mitk::BoolProperty::New(false));
    this->GetDataStorage()->Add(newNode, baseNode);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();

    // show gaussianed results
    if (bOutputGaussian)
    {
        typedef std::multimap<double, int, std::greater<double>> SORT_MAP;
        SORT_MAP sort_map;
        for (unsigned int n = 0; n < szNodules; ++n){
            sort_map.insert(std::make_pair(dwNoduleScore[n], n));
        }

        SORT_MAP::iterator sort_map_it = sort_map.begin(), sort_map_end = sort_map.end();
        for (int i = 0; i < nShownResults; i++)
        {
            SVisualizeNodule& this_nodule = vVisualizeNodules[sort_map_it->second];
            VIndexValue& this_index_values = this_nodule.nodule_index_values;
            VIndexes& this_indexes = this_nodule.nodule_indexes;

            for_each(this_indexes.begin(), this_indexes.end(),
                [&itkDupImage](ImageType::IndexType indexes)
            {
                itkDupImage->SetPixel(indexes, 255);
            });

            vOutTotalNoduleIndex.push_back(sort_map_it->second);
            sort_map_it++;
        }

        mitk::Image::Pointer mitkGaussianProcessed;
        mitk::CastToMitkImage(itkDupImage, mitkGaussianProcessed);
        strFormat.sprintf("%s", "ANNGaussianResults_");
        strFormat += STR_ANN_RESULTS_TYPE[result_type];
        ShowProcessedDataNode(mitkGaussianProcessed, strFormat.toStdString(), false, baseNode);
    }

    // show built gaussian kernel
    if (bOutputGGroundtruth)
    {
        MultiplyCastType::Pointer multiplyFilter2 = MultiplyCastType::New();
        multiplyFilter2->SetInput(itkDupImage2);
        multiplyFilter2->SetConstant(255);
        multiplyFilter2->Update();

        mitk::Image::Pointer mitkGaussianProcessed;
        mitk::CastToMitkImage(multiplyFilter2->GetOutput(), mitkGaussianProcessed);
        strFormat.sprintf("%s", "ANNGaussians_");
        strFormat += STR_ANN_RESULTS_TYPE[result_type];
        ShowProcessedDataNode(mitkGaussianProcessed, strFormat.toStdString(), false, baseNode);
    }

    // showing in tree control
    if (bOutputGaussian)
    {
        QTreeWidgetItem *imageItem1 = new QTreeWidgetItem(m_Controls.treeWidget, QStringList(STR_ANN_RESULTS_TYPE[result_type]));
        for (int i = 0; i < szNodules; i++)
        {
            QString str_data;

            SVisualizeNodule& this_nodule = vVisualizeNodules[i];
            ImageType::IndexType this_index = this_nodule.nodule_region.GetIndex();
            ImageType::SizeType this_size = this_nodule.nodule_region.GetSize();

            QTreeWidgetItem *imageItem1_1 = new QTreeWidgetItem(imageItem1, QStringList(QString("Band1")));
            str_data.sprintf("%d %d %d", this_index[0], this_index[1], this_index[2]);
            imageItem1_1->setText(0, QString::number(i));
            imageItem1_1->setText(1, str_data);
            str_data.sprintf("%d %d %d", this_size[0], this_size[1], this_size[2]);
            imageItem1_1->setText(2, str_data);
            str_data.sprintf("%f", dwNoduleScore[i]);
            imageItem1_1->setText(3, str_data);
            imageItem1->addChild(imageItem1_1);
        }

    }
    SAFE_DELETE(dwNoduleScore);
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoTakeScreenShots()
{
    vtkRenderer* renderer = m_IRenderWindowPart->GetActiveQmitkRenderWindow()->GetRenderer()->GetVtkRenderer();
    QString fileName = "test";
    if (renderer != nullptr)
        this->TakeScreenshot(renderer, 1, m_WorkingDir + fileName);
}


//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::TakeScreenshot(vtkRenderer* renderer, unsigned int magnificationFactor, QString fileName)
{
    if ((renderer == nullptr) || (magnificationFactor < 1) || fileName.isEmpty())
        return;

    bool doubleBuffering(renderer->GetRenderWindow()->GetDoubleBuffer());
    renderer->GetRenderWindow()->DoubleBufferOff();

    vtkImageWriter* fileWriter;

    QFileInfo fi(fileName);
    QString suffix = fi.suffix();
    if (suffix.compare("jpg", Qt::CaseInsensitive) == 0)
    {
        vtkJPEGWriter* w = vtkJPEGWriter::New();
        w->SetQuality(100);
        w->ProgressiveOff();
        fileWriter = w;

    }
    else  // default is png
    {
        fileWriter = vtkPNGWriter::New();
    }
    vtkRenderLargeImage* magnifier = vtkRenderLargeImage::New();
    magnifier->SetInput(renderer);
    magnifier->SetMagnification(magnificationFactor);
    //magnifier->Update();
    fileWriter->SetInputConnection(magnifier->GetOutputPort());
    fileWriter->SetFileName(fileName.toLatin1());

    // vtkRenderLargeImage has problems with different layers, therefore we have to
    // temporarily deactivate all other layers.
    // we set the background to white, because it is nicer than black...
    double oldBackground[3];
    renderer->GetBackground(oldBackground);

    //  QColor color = QColorDialog::getColor();
    double bgcolor[] = { m_BackgroundColor.red() / 255.0, m_BackgroundColor.green() / 255.0, m_BackgroundColor.blue() / 255.0 };
    renderer->SetBackground(bgcolor);

    fileWriter->Write();
    fileWriter->Delete();
    renderer->SetBackground(oldBackground);

    renderer->GetRenderWindow()->SetDoubleBuffer(doubleBuffering);
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::SelectBackgroundColor()
{
    m_BackgroundColor = QColorDialog::getColor();

    m_Controls.btBGColor->setAutoFillBackground(true);

    QString styleSheet = "background-color:rgb(";
    styleSheet.append(QString::number(m_BackgroundColor.red()));
    styleSheet.append(",");
    styleSheet.append(QString::number(m_BackgroundColor.green()));
    styleSheet.append(",");
    styleSheet.append(QString::number(m_BackgroundColor.blue()));
    styleSheet.append(")");
    m_Controls.btBGColor->setStyleSheet(styleSheet);
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoOutputScoringResults()
{
    if (m_vvTotalNoduleIndex.size() != EANNRT_COUNT)
    {
        MITK_ERROR << "m_vvTotalNoduleIndex Size NOT Equal With EANNRT_COUNT";
        return;
    }
    OutputScoringResultImages(m_vvTotalNoduleIndex);
}


//////////////////////////////////////////////////////////////////////////
bool QmitkRegionGrowingView::ClearDirectory(QString folderDir)
{
    QDir dir(folderDir);
    QFileInfoList fileList;
    QFileInfo curFile;
    if (!dir.exists())  { return false; }//文件不存，则返回false
    fileList = dir.entryInfoList(QDir::Dirs | QDir::Files
        | QDir::Readable | QDir::Writable
        | QDir::Hidden | QDir::NoDotAndDotDot
        , QDir::Name);
    while (fileList.size() > 0)//跳出条件
    {
        int infoNum = fileList.size();
        for (int i = infoNum - 1; i >= 0; i--)
        {
            curFile = fileList[i];
            if (curFile.isFile())//如果是文件，删除文件
            {
                QFile fileTemp(curFile.filePath());
                fileTemp.remove();
                fileList.removeAt(i);
            }
            if (curFile.isDir())//如果是文件夹
            {
                QDir dirTemp(curFile.filePath());
                QFileInfoList fileList1 = dirTemp.entryInfoList(QDir::Dirs | QDir::Files
                    | QDir::Readable | QDir::Writable
                    | QDir::Hidden | QDir::NoDotAndDotDot
                    , QDir::Name);
                if (fileList1.size() == 0)//下层没有文件或文件夹
                {
                    dirTemp.rmdir(".");
                    fileList.removeAt(i);
                }
                else//下层有文件夹或文件
                {
                    for (int j = 0; j < fileList1.size(); j++)
                    {
                        if (!(fileList.contains(fileList1[j])))
                            fileList.append(fileList1[j]);
                    }
                }
            }
        }
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::OutputScoringResultImages(VVecInt& vVNodules)
{
    QString qstr_debug;

    if (m_vvTotalNoduleIndex.empty())
    {
        MITK_ERROR << "Data container empty";
        return;
    }

    if (m_WorkingDir == "")
    {
        MITK_ERROR << "Select Working Directory First!!";
        return;
    }

    // load init files
    VVisualizeNodules vVisualizeNodules;
    ImageType::SizeType image_size;

    bool bLoad3 = LoadRegionParams(m_WorkingDir + G_ANN_REGION_FILE_NAME, vVisualizeNodules);
    bool bLoad2 = LoadIndexes(m_WorkingDir + G_ANN_INDEX_FILE_NAME, image_size, vVisualizeNodules);

    if (!bLoad2){
        MITK_WARN << "Load Indexes Failed";
        return;
    }
    if (!bLoad3){
        MITK_WARN << "Load Regions Failed";
        return;
    }

    QDir q_dir;
    QString q_working_dir;

    // get texture node and calculate offset 
    mitk::DataStorage::Pointer pDataStorage = GetDataStorage();
    mitk::DataNode::Pointer texture_image_node = pDataStorage->GetNamedNode("image_extract_bg_-2048");
    if (texture_image_node == nullptr){
        MITK_ERROR << "No image_extract_bg_-2048 node found!!";
        return;
    }
    mitk::BaseData*	data = texture_image_node->GetData();
    mitk::Image::Pointer mitkTexture = NULL;
    mitk::BaseGeometry::Pointer geometry = NULL;
    if (data)
        mitkTexture = dynamic_cast<mitk::Image*>(data);
    else
        return;

    // transformed image offset
    geometry = mitkTexture->GetGeometry();
    mitk::AffineTransform3D::Pointer transform = geometry->GetIndexToWorldTransform();
    mitk::AffineTransform3D::OutputVectorType offset3D = transform->GetOffset();
    //mitk::Point3D offset3D = geometry->GetOrigin();

    float nDim0_0 = offset3D[0];
    float nDim0_1 = offset3D[1];
    float nDim0_2 = offset3D[2];
    MITK_INFO << "From Origin: " << nDim0_0 << " " << nDim0_1 << " " << nDim0_2;

    texture_image_node = pDataStorage->GetNamedNode("normalized_texture");
    if (texture_image_node == nullptr){
        MITK_ERROR << "No nodule_texture node found!!";
        return;
    }
    data = texture_image_node->GetData();
    mitkTexture = NULL;
    if (data)
        mitkTexture = dynamic_cast<mitk::Image*>(data);
    else
        return;

    // original image offset
    geometry = mitkTexture->GetGeometry();
    transform = geometry->GetIndexToWorldTransform();
    offset3D = transform->GetOffset();
    //offset3D = geometry->GetOrigin();

    float nDim1_0 = offset3D[0];
    float nDim1_1 = offset3D[1];
    float nDim1_2 = offset3D[2];
    MITK_INFO << "To Origin: " << nDim1_0 << " " << nDim1_1 << " " << nDim1_2;

    // relative substract
    ImageType::IndexType offset_index;
    offset_index[0] = std::round(nDim0_0 - nDim1_0);
    offset_index[1] = std::round(nDim0_1 - nDim1_1);
    offset_index[2] = std::round(nDim0_2 - nDim1_2);

    MITK_INFO << "Offset: " << offset_index[0] << " " << offset_index[1] << " " << offset_index[2];

    /// arrange output dirs
    q_working_dir = m_WorkingDir + "annResults";
    if (!q_dir.exists(q_working_dir)){
        q_dir.mkdir(q_working_dir);
    }
    else{
        ClearDirectory(q_working_dir);
    }

    for (int i = 0; i < EANNResultsType::EANNRT_COUNT; i++)
    {
        // for each ann type
        VecInt& vResultIndex = vVNodules[i];

        q_working_dir = m_WorkingDir + "annResults" + "/" + STR_ANN_RESULTS_TYPE[i];

        // no need for checking anymore since all data are deleted in last step.
        //         if (q_dir.exists(q_working_dir)){
        //             q_dir.rmdir(q_working_dir);
        //         }
        q_dir.mkdir(q_working_dir);

        if (vResultIndex.empty())
        {
            MITK_WARN << "No Output Indexes for " << STR_ANN_RESULTS_TYPE[i] << " Node";
            continue;
        }

        // get result node
        mitk::DataNode::Pointer gaussian_image_node = pDataStorage->GetNamedNode(QString("ANNResults_" + STR_ANN_RESULTS_TYPE[i]).toStdString());
        if (texture_image_node == nullptr){
            MITK_WARN << "No ANNResults_" << STR_ANN_RESULTS_TYPE[i] << " Node found !!";
            continue;
        }

        mitk::BaseData*	gaussian_data = gaussian_image_node->GetData();
        mitk::Image::Pointer mitkGaussian = NULL;
        if (data)
            mitkGaussian = dynamic_cast<mitk::Image*>(gaussian_data);
        else
            continue;


        for (unsigned int j = 0; j < vResultIndex.size(); j++)
        {
            int this_output_index = vResultIndex[j];

            if (this_output_index >= vVisualizeNodules.size())
            {
                qstr_debug.sprintf("output_index(%d) out of vVisualizeNodules size(%d)", this_output_index, vVisualizeNodules.size());
                MITK_ERROR << qstr_debug.toStdString();
                continue;
            }

            SVisualizeNodule& this_output_nodule = vVisualizeNodules[this_output_index];
            QString this_working_dir;
            this_working_dir.sprintf("%s/%d", q_working_dir.toStdString().c_str(), j);
            q_dir.mkdir(this_working_dir);
            this_working_dir += "/";

            ExportROIImage(this_working_dir, mitkTexture, mitkGaussian, this_output_nodule.nodule_region/*export_roi_reigon*/, offset_index);

            qstr_debug.sprintf("Nodule %d in %s ExportROIImage Done !!!", j, STR_ANN_RESULTS_TYPE[i]);
            MITK_INFO << qstr_debug.toStdString();
        }
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::ExportROIImage(QString outputFileName,
    mitk::Image::Pointer& texture_image, mitk::Image::Pointer& result_image,
    ImageType::RegionType& region_roi, ImageType::IndexType& offset)
{
    QString str_DBG;

    typedef itk::ExtractImageFilter< ImageType, ImageType2D >               ExtractFilterType;
    typedef itk::RegionOfInterestImageFilter< ImageType2D, ImageType2D >    ROIFilterType;
    typedef itk::CastImageFilter< ImageType2D, BinaryImageType2D >          CastFilterType;
    typedef itk::ImageFileWriter< BinaryImageType2D >                       WriterType;

    ImageType::Pointer itkTexture = NULL;
    ImageType::Pointer itkGaussian = NULL;

    mitk::CastToItkImage(texture_image, itkTexture);
    mitk::CastToItkImage(result_image, itkGaussian);

    // roi region  
    ImageType::IndexType   index_value = region_roi.GetIndex();
    ImageType::SizeType    size_value = region_roi.GetSize();

    // find iterator center
    ImageType::SizeType neighbour_iterator_size;
    neighbour_iterator_size.Fill(1);
    ImageType::IndexType neighbour_iterator_index;
    ImageType::RegionType neighbour_iterator_region;
    neighbour_iterator_region.SetSize(neighbour_iterator_size);

    // create the local working image
    BinaryImageType2D::Pointer itkLocalImage = BinaryImageType2D::New();
    BinaryImageType2D::IndexType local_start;
    local_start[0] = 0;
    local_start[1] = 0;
    BinaryImageType2D::SizeType  local_size;
    local_size[0] = N_EXPORT_DIM;
    local_size[1] = N_EXPORT_DIM;
    BinaryImageType2D::RegionType local_region;
    local_region.SetSize(local_size);
    local_region.SetIndex(local_start);
    itkLocalImage->SetRegions(local_region);
    itkLocalImage->Allocate();

    ImageType::SizeType radius;
    radius[0] = N_EXPORT_DIM / 2;
    radius[1] = N_EXPORT_DIM / 2;
    radius[2] = 0;

    WriterType::Pointer writer1 = WriterType::New();

    //slice by slice
    QString q_tmp_path;
    int szSlices = size_value[2];
    for (int iSlice = 0; iSlice < szSlices; iSlice++)
    {
        //-------------------------------------------------------------------------
        // texture
        neighbour_iterator_index[0] = index_value[0] + offset[0] + size_value[0] / 2;
        neighbour_iterator_index[1] = index_value[1] + offset[1] + size_value[1] / 2;
        neighbour_iterator_index[2] = index_value[2] + offset[2] + iSlice;
        neighbour_iterator_region.SetIndex(neighbour_iterator_index);

        // clear working image
        itkLocalImage->FillBuffer(0);

        itk::NeighborhoodIterator<ImageType> texture_iterator(radius, itkTexture, neighbour_iterator_region);
        while (!texture_iterator.IsAtEnd())
        {
            // Get the value of the current pixel
            ImageType::IndexType center_index = texture_iterator.GetIndex();
            for (unsigned int i = 0; i < texture_iterator.Size(); i++)
            {
                ImageType::IndexType local_index = texture_iterator.GetIndex(i);
                bool IsInBounds;
                int neighborValue = texture_iterator.GetPixel(i, IsInBounds);
                if (!IsInBounds)
                    neighborValue = -(int)F_CT_WIDTH;

                limit_number<int>(neighborValue, -F_CT_WIDTH, F_CT_WIDTH);
                unsigned char chNormalized = std::round((double(neighborValue + F_CT_WIDTH)) / (F_CT_WIDTH * 2) * 255);

                BinaryImageType2D::IndexType local_image_index;
                local_image_index[0] = local_index[0] - neighbour_iterator_index[0] + radius[0];
                local_image_index[1] = local_index[1] - neighbour_iterator_index[1] + radius[1];
                limit_number<BinaryImageType2D::IndexType::IndexValueType>(local_image_index[0], 0, N_EXPORT_DIM - 1);
                limit_number<BinaryImageType2D::IndexType::IndexValueType>(local_image_index[1], 0, N_EXPORT_DIM - 1);

                itkLocalImage->SetPixel(local_image_index, chNormalized);
            }
            ++texture_iterator;
        }
        q_tmp_path.sprintf("%s%s%d.png", outputFileName.toStdString().c_str(), "texture_", index_value[2] + offset[2] + iSlice);
        writer1->SetFileName(q_tmp_path.toStdString());
        writer1->SetInput(itkLocalImage);
        writer1->Update();

        //-------------------------------------------------------------------------
        // gaussian result
        neighbour_iterator_index[0] = index_value[0] + size_value[0] / 2;
        neighbour_iterator_index[1] = index_value[1] + size_value[1] / 2;
        neighbour_iterator_index[2] = index_value[2] + iSlice;
        neighbour_iterator_region.SetIndex(neighbour_iterator_index);

        // clear working image
        itkLocalImage->FillBuffer(0);

        itk::NeighborhoodIterator<ImageType> gaussian_iterator(radius, itkGaussian, neighbour_iterator_region);
        while (!gaussian_iterator.IsAtEnd())
        {
            // Get the value of the current pixel
            ImageType::IndexType center_index = gaussian_iterator.GetIndex();
            for (unsigned int i = 0; i < gaussian_iterator.Size(); i++)
            {
                ImageType::IndexType local_index = gaussian_iterator.GetIndex(i);
                bool IsInBounds;
                int neighborValue = gaussian_iterator.GetPixel(i, IsInBounds);
                if (!IsInBounds)
                    neighborValue = 0;

                BinaryImageType2D::IndexType local_image_index;
                local_image_index[0] = local_index[0] - neighbour_iterator_index[0] + radius[0];
                local_image_index[1] = local_index[1] - neighbour_iterator_index[1] + radius[1];

                limit_number<BinaryImageType2D::IndexType::IndexValueType>(local_image_index[0], 0, N_EXPORT_DIM - 1);
                limit_number<BinaryImageType2D::IndexType::IndexValueType>(local_image_index[1], 0, N_EXPORT_DIM - 1);
                itkLocalImage->SetPixel(local_image_index, neighborValue);
            }
            ++gaussian_iterator;
        }
        q_tmp_path.sprintf("%s%s%d.png", outputFileName.toStdString().c_str(), "gaussian_", index_value[2] + iSlice);
        writer1->SetFileName(q_tmp_path.toStdString());
        writer1->SetInput(itkLocalImage);
        writer1->Update();
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::OutputBOWSampleData()
{
    // for debugging formation
    QString strQDebug;

    // validate data 
    mitk::DataNode::Pointer textureNode = m_Controls.cbOriginalImage_2->GetSelectedNode();
    mitk::DataNode::Pointer maskNode = m_Controls.cbMaskImage_2->GetSelectedNode();

    if (NULL == textureNode || NULL == maskNode){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    //-------------------------------------------------------------------------
    // image validation
    mitk::Image::Pointer	mitkTexture = dynamic_cast<mitk::Image*> (textureNode->GetData());
    mitk::Image::Pointer	mitkMask = dynamic_cast<mitk::Image*>(maskNode->GetData());
    // check if images are valid
    if ((!mitkTexture) || (!mitkMask) || (mitkTexture->IsInitialized() == false) || (mitkMask->IsInitialized() == false)){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    //-------------------------------------------------------------------------
    // path and image preparation
    QString strSavePath = m_WorkingDir;

//     mitk::Image::Pointer mitkForecastImage = mitkTexture->Clone();
//     ImageType::Pointer itkForecast;
//     mitk::CastToItkImage(mitkForecastImage, itkForecast);
//     itkForecast->FillBuffer(0);

    // transform the mitk texture
    ImageType::Pointer itkTexture;
    mitk::CastToItkImage(mitkTexture, itkTexture);
    BinaryImageType::Pointer itkMask;
    mitk::CastToItkImage(mitkMask, itkMask);

    //-------------------------------------------------------------------------
    // start
    int iWidth = mitkTexture->GetDimension(0);
    int iHeight = mitkTexture->GetDimension(1);
    int iSlice = mitkTexture->GetDimension(2);

    // file header: sz_data szInput szOutput
    // 	int szData = N_NODULE_D * N_NODULE_D * N_NODULE_D;
    //     int szInput = N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D;
    //     int szOutput = 1;
    //     int szAllData = 0;
    //     ImageType::SizeType image_size;
    //     image_size[0] = iWidth;
    //     image_size[1] = iHeight;
    //     image_size[2] = iSlice;

    // collect regions
    VVisualizeNodules vNoduleRegions;
    collectRegions(itkMask, vNoduleRegions);

    int szAllData = vNoduleRegions.size();
    //     for (int iVCenter = 0; iVCenter < szAllData; iVCenter++)
    //     {
    //         SVisualizeNodule region_infor = vNoduleRegions[iVCenter];
    //         ImageType::RegionType nodule_region = region_infor.nodule_region;
    //         ImageType::IndexType nodule_index = nodule_region.GetIndex();
    //         ImageType::SizeType nodule_size = nodule_region.GetSize();
    // 
    //         int iRegionX = std::min((int)nodule_size[0], N_NODULE_D);
    //         int iRegionY = std::min((int)nodule_size[1], N_NODULE_D);
    //         int iRegionZ = std::min((int)nodule_size[2], N_NODULE_D);
    //         strQDebug.sprintf("RegionX %d, RegionY %d, RegionZ %d", iRegionX, iRegionY, iRegionZ);
    //         MITK_INFO << << strQDebug.toStdString();
    // //         szData += iRegionX*iRegionY*iRegionZ;
    //     }



    QDir q_dir;
    QString qSTRTestDir = "";

    qSTRTestDir = m_WorkingDir + "bow";
    if (!q_dir.exists(qSTRTestDir)){
        q_dir.mkdir(qSTRTestDir);
    }
    else{
        ClearDirectory(qSTRTestDir);
    }

    qSTRTestDir = m_WorkingDir + "bow\\testing";
    if (!q_dir.exists(qSTRTestDir)){
        q_dir.mkdir(qSTRTestDir);
    }
    else{
        ClearDirectory(qSTRTestDir);
    }

    for (int iVCenter = 0; iVCenter < szAllData; iVCenter++)
    {
        QString file_tmp;
        file_tmp.sprintf("%d.txt", iVCenter);
        qSTRTestDir = m_WorkingDir + "bow\\testing\\" + file_tmp;
        MITK_INFO << "Saving bow sample to " << qSTRTestDir.toStdString();

        QFile qBOW(qSTRTestDir);
        if (qBOW.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
        {
            QTextStream tsBOW(&qBOW);
            tsBOW.setRealNumberNotation(QTextStream::FixedNotation);
            tsBOW.setRealNumberPrecision(PIXEL_PRECISION);

            SVisualizeNodule region_infor = vNoduleRegions[iVCenter];
            ImageType::RegionType noduleRegion = region_infor.nodule_region;
            ImageType::IndexType noduleIndex = noduleRegion.GetIndex();
            ImageType::SizeType noduleSize = noduleRegion.GetSize();

            ImageType::IndexType centerIndex;
            centerIndex[0] = noduleIndex[0] + noduleSize[0] / 2;
            centerIndex[1] = noduleIndex[1] + noduleSize[1] / 2;
            centerIndex[2] = noduleIndex[2] + noduleSize[2] / 2;

            //             int iRegionX = std::min((int)noduleSize[0], N_NODULE_D);
            //             int iRegionY = std::min((int)noduleSize[1], N_NODULE_D);
            //             int iRegionZ = std::min((int)noduleSize[2], N_NODULE_D);

            ImageType::SizeType regionSize;
            regionSize[0] = N_NODULE_D;
            regionSize[1] = N_NODULE_D;
            regionSize[2] = N_NODULE_D;

            ImageType::IndexType regionIndex;
            regionIndex[0] = centerIndex[0] - N_NODULE_D / 2 - 1;
            regionIndex[1] = centerIndex[1] - N_NODULE_D / 2 - 1;
            regionIndex[2] = centerIndex[2] - N_NODULE_D / 2 - 1;

            // start at index and look for neighbors at each point inside of index+size;
            ImageType::RegionType region;
            // region size
            region.SetSize(regionSize);
            // region start
            region.SetIndex(regionIndex);

            ImageType::SizeType radius;
            radius[0] = 0;
            radius[1] = 0;
            radius[2] = 0;

            itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
            // one piece of training data
            while (!iterator.IsAtEnd())
            {
                int nTotalIn = 0;
                ImageType::IndexType center_index = iterator.GetIndex();

                for (unsigned int i = 0; i < iterator.Size(); i++)
                {
                    ImageType::IndexType local_index = iterator.GetIndex(i);

                    bool IsInBounds;
                    int neighborValue = iterator.GetPixel(i, IsInBounds);

                    if (!IsInBounds){
                        neighborValue = -(int)F_CT_WIDTH;
                    }
                    if (neighborValue > F_CT_WIDTH){
                        neighborValue = F_CT_WIDTH;
                    }
                    if (neighborValue < -F_CT_WIDTH){
                        neighborValue = -F_CT_WIDTH;
                    }
//                     if (IsInBounds){
//                         itkForecast->SetPixel(local_index, neighborValue);
//                     }
                    double dwNormalized = ((double)(neighborValue + F_CT_WIDTH)) / ((double)(F_CT_WIDTH * 2));
                    tsBOW << dwNormalized/*neighborValue */ << "\n";
                }

                ++iterator;
            }
            qBOW.close();
        }
        else
        {
            MITK_WARN << "Open File Failed...";
        }
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoOutputBOWSamples()
{
    OutputBOWSampleData();
}


//////////////////////////////////////////////////////////////////////////
#include "itkBinaryMask3DMeshSource.h"
#include "itkQuadEdgeMesh.h"
#include "itkMeshFileWriter.h"
#include "itkQuadEdgeMeshDecimationCriteria.h"
#include "itkSquaredEdgeLengthDecimationQuadEdgeMeshFilter.h"

void QmitkRegionGrowingView::DoReduceFaces()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        ImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        typedef itk::QuadEdgeMesh< double, 3 > MeshType;

        typedef itk::BinaryMask3DMeshSource< ImageType, MeshType > FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(itkImage);
        filter->SetObjectValue(255);
        filter->Update();

        MeshType::Pointer mesh = filter->GetOutput();

        typedef itk::NumberOfFacesCriterion< MeshType > CriterionType;
        typedef itk::SquaredEdgeLengthDecimationQuadEdgeMeshFilter<MeshType, MeshType, CriterionType > DecimationType;

        QString qValue = m_Controls.leMeshFaceCount->text();
        long int N = qValue.toInt();

        CriterionType::Pointer criterion = CriterionType::New();
        criterion->SetTopologicalChange(true);
        criterion->SetNumberOfElements(N);

        DecimationType::Pointer decimate = DecimationType::New();
        decimate->SetInput(mesh);
        decimate->SetCriterion(criterion);
        decimate->Update();

        typedef itk::MeshFileWriter< MeshType > WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName("vtkmesh_data.vtk");
        writer->SetInput(decimate->GetOutput());
        try
        {
            writer->Update();
        }
        catch (itk::ExceptionObject & error)
        {
            MITK_ERROR(error.GetDescription());
        }
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoMethodTest1()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        ImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        typedef itk::SliceBySliceImageFilter< ImageType, ImageType> SbSFilterType;
        SbSFilterType::Pointer SbSFilter = SbSFilterType::New();

        typedef itk::CurvatureFlowImageFilter< ImageType2D, ImageType2D> SmoothingFilter;
        SmoothingFilter::Pointer smoothingFilter = SmoothingFilter::New();
        smoothingFilter->SetTimeStep(0.002);
        smoothingFilter->SetNumberOfIterations(200);

        SbSFilter->SetFilter(smoothingFilter);
        SbSFilter->SetInput(itkImage);
        SbSFilter->Update();

        QString strNode;
        strNode.sprintf("%s_curv", mitkNode->GetName().c_str());
        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(SbSFilter->GetOutput(), mitkResult);
        ShowProcessedDataNode(mitkResult, strNode.toStdString(), false, mitkNode);
    }

}

//////////////////////////////////////////////////////////////////////////
#include "itkBinaryContourImageFilter.h"
void QmitkRegionGrowingView::DoMethodTest2()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        ImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        typedef itk::BinaryContourImageFilter <ImageType, ImageType > Filter;

        Filter::Pointer filter = Filter::New();
        filter->SetInput(itkImage);
        filter->SetBackgroundValue(0);
        filter->SetForegroundValue(255);
        filter->Update();

        QString strNode;
        strNode.sprintf("%s_contour", mitkNode->GetName().c_str());
        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(filter->GetOutput(), mitkResult);
        ShowProcessedDataNode(mitkResult, strNode.toStdString(), false, mitkNode);
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoMethodTest3()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();

    int iWidth = mitkImage->GetDimension(0);
    int iHeight = mitkImage->GetDimension(1);
    int iSlice = mitkImage->GetDimension(2);

    QString qValue = m_Controls.leSliceNumber->text();
    long int N = qValue.toInt();

    if (NULL != mitkImage)
    {
        ImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        ImageType::IndexType desiredStart;
        desiredStart[0] = 0;
        desiredStart[1] = 0;
        desiredStart[2] = N;

        ImageType::SizeType desiredSize;
        desiredSize[0] = iWidth;
        desiredSize[1] = iHeight;
        desiredSize[2] = 0;

        ImageType::RegionType desiredRegion;
        desiredRegion.SetIndex(desiredStart);
        desiredRegion.SetSize(desiredSize);

        // declare an extract filter to extract slices
        typedef itk::ExtractImageFilter< ImageType, ImageType2D > ExtractSliceFilterType;
        ExtractSliceFilterType::Pointer  filter = ExtractSliceFilterType::New();
        filter->SetInput(itkImage);
        filter->SetExtractionRegion(desiredRegion);
        filter->SetDirectionCollapseToIdentity();
        filter->Update();

        QString strNode;
        strNode.sprintf("%s_extracted", mitkNode->GetName().c_str());
        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(filter->GetOutput(), mitkResult);
        ShowProcessedDataNode(mitkResult, strNode.toStdString(), false, mitkNode);
    }
}

//////////////////////////////////////////////////////////////////////////
#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkContourExtractor2DImageFilter.h"
void QmitkRegionGrowingView::DoMethodTest4()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();
    if (NULL != mitkImage)
    {
        ImageType2D::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        typedef  itk::ApproximateSignedDistanceMapImageFilter< ImageType2D, FloatImageType2D > ApproximateSignedDistanceMapImageFilterType;
        ApproximateSignedDistanceMapImageFilterType::Pointer approximateSignedDistanceMapImageFilter = ApproximateSignedDistanceMapImageFilterType::New();
        approximateSignedDistanceMapImageFilter->SetInput(itkImage);
        approximateSignedDistanceMapImageFilter->SetInsideValue(255);
        approximateSignedDistanceMapImageFilter->SetOutsideValue(0);
        approximateSignedDistanceMapImageFilter->Update();

        typedef itk::ContourExtractor2DImageFilter <FloatImageType2D> ContourExtractor2DImageFilterType;
        ContourExtractor2DImageFilterType::Pointer contourExtractor2DImageFilter = ContourExtractor2DImageFilterType::New();
        contourExtractor2DImageFilter->SetInput(approximateSignedDistanceMapImageFilter->GetOutput());
        contourExtractor2DImageFilter->SetContourValue(0);
        contourExtractor2DImageFilter->Update();


        char subfix[MAX_PATH] = "E:\\shape_analysis\\";
        QString qValue = m_Controls.leSliceNumber->text();
        long int N = qValue.toInt();
        QString tmpPath;
        tmpPath.sprintf("%sdata_%d.txt", subfix, N);

        qValue = m_Controls.leContourMinPointCount->text();
        int nMinContourCount = qValue.toInt();

        QFile qTraining(tmpPath);
        QTextStream tsTraining(&qTraining);
        if (qTraining.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
        {
            tsTraining.setRealNumberNotation(QTextStream::FixedNotation);
            tsTraining.setRealNumberPrecision(6);
        }
        else
            return;

        MITK_INFO << "Writing to " << tmpPath;
        MITK_INFO << "There are " << contourExtractor2DImageFilter->GetNumberOfOutputs() << " contours" << std::endl;
        for (unsigned int i = 0; i < contourExtractor2DImageFilter->GetNumberOfOutputs(); ++i)
        {
            int size_this_contour = contourExtractor2DImageFilter->GetOutput(i)->GetVertexList()->Size();
            if (size_this_contour < nMinContourCount)
            {
                MITK_INFO << "Contour " << i << ": " << size_this_contour << " skipped";
                continue;
            }

            ContourExtractor2DImageFilterType::VertexListType::ConstIterator vertexIterator =
                contourExtractor2DImageFilter->GetOutput(i)->GetVertexList()->Begin();
            while (vertexIterator != contourExtractor2DImageFilter->GetOutput(i)->GetVertexList()->End())
            {
                tsTraining << vertexIterator->Value()[0] << " " << vertexIterator->Value()[1] << "\n";
                ++vertexIterator;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////
void QmitkRegionGrowingView::DoRefineRegionGrowing()
{
    QString qStrDBG;
    mitk::Image::Pointer	mitkImage = getSelectedImage();
    mitk::DataNode::Pointer	mitkNode = getSelectedNode();

    if (NULL != mitkImage)
    {
        BinaryImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage, itkImage);

        typedef itk::VotingBinaryIterativeHoleFillingImageFilter< BinaryImageType > FilterType;
        FilterType::InputSizeType radius;
        radius.Fill(15);
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(itkImage);
        filter->SetRadius(radius);
        filter->SetMajorityThreshold(2);
        filter->SetBackgroundValue(0);
        filter->SetForegroundValue(1);
        filter->SetMaximumNumberOfIterations(15);
        filter->Update();

        QString strNode;
        strNode.sprintf("%s_bin_fill", mitkNode->GetName().c_str());
        mitk::Image::Pointer mitkResult;
        mitk::CastToMitkImage(filter->GetOutput(), mitkResult);
        ShowProcessedDataNode(mitkResult, strNode.toStdString(), false, mitkNode);
    }
}