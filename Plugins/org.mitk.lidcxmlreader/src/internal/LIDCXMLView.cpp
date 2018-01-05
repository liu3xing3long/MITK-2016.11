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

// Qmitk
#include "LIDCXMLView.h"
#include "org_mitk_lidcxmlreader_Activator.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QStandardItemModel>
#include <QXmlStreamReader>
#include <QPixmap>
#include <QPainter>
#include "MyDrawableNodule.h"

//mitk image
#include "mitkImage.h"
#include "mitkImageCast.h"
#include "mitkImagePixelReadAccessor.h"
#include "mitkImagePixelWriteAccessor.h"
#include "mitkIOUtil.h"
#include "mitkNodePredicateDataType.h"
#include "mitkNodePredicateDimension.h"
#include "mitkNodePredicateAnd.h"
#include "mitkDataNodeObject.h"

//itk&&vtk headers
#include <mitkImageVtkReadAccessor.h>
#include <mitkImageVtkWriteAccessor.h>

// itk headers
#include <itkNeighborhoodIterator.h>
#include <itkImageRegionIterator.h>
#include <itkCastImageFilter.h>
#include <itkResampleImageFilter.h>
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkShapeOpeningLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkLabelToRGBImageFilter.h"

// cuda headers
// #include "CudaBinaryThresholdImageFilter.h"
// #include "CudaGrayscaleDilateImageFilter.h"
// #include "CudaMedianImageFilter.h"
// #include "CudaVesselnessImageFilter.h"
// #include "CudaMultiplyByConstantImageFilter.h"
// #include "CudaDivideByConstantImageFilter.h"
// #include "CudaGrayscaleErodeImageFilter.h"

// ITK
#include "itkMultiplyImageFilter.h"
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

#include "MyMath.h"
#include <mitkIDataStorageService.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkProperties.h>

#include <berryIEditorReference.h>
#include <berryIWorkbenchPage.h>
#include <berryIWorkbenchWindow.h>
#include <mitkWorkbenchUtil.h>
#include <QmitkIOUtil.h>

const std::string LIDCXMLView::VIEW_ID = "org.mitk.views.lidcxmlview";
const unsigned int RESULT_PRECISION = 12;
const unsigned int PIXEL_PRECISION = 4;

//////////////////////////////////////////////////////////////////////////
enum EOutputNoduleType
{
    EOutType_ISOLATED = 0,
    EOutType_WALL,
    EOutType_VESSEL,
    EOutType_GGO,
    EOutType_TRAIL,
    //EOutType_NONNODULE,
    EOutType_COUNT,
};

QString STR_OUTPUT_TYPE[EOutputNoduleType::EOutType_COUNT] = {
    "ISO",
    "WALL",
    "VESSEL",
    "GGO",
    "TRAIL",
    // "NONNODULE",
};


//////////////////////////////////////////////////////////////////////////
LIDCXMLView::LIDCXMLView()
{
    iCurrengShownImage = -1;
    m_WorkingDir = "";
    m_outputNoduleTypes = NULL;
    m_bgImage = NULL;
}

//////////////////////////////////////////////////////////////////////////
LIDCXMLView::~LIDCXMLView()
{
    destroyContainers();
}


//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::SetFocus()
{
    m_Controls.btLoadAnnotation->setFocus();

}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::CreateQtPartControl(QWidget *parent)
{
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.btLoadAnnotation, SIGNAL(clicked()), this, SLOT(DoLoadAnnotations()));
    connect(m_Controls.btnPriorImage, SIGNAL(clicked()), this, SLOT(showPriorImage()));
    connect(m_Controls.btnNextImage, SIGNAL(clicked()), this, SLOT(showNextImage()));
    // 	connect(m_Controls.btnSaveSpecificImage, SIGNAL(clicked()), this, SLOT(saveSpecificImage()));
    connect(m_Controls.btOutputtraining, SIGNAL(clicked()), this, SLOT(DoOutputTraining()));
    connect(m_Controls.btColorImage, SIGNAL(clicked()), this, SLOT(DoColorLabel()));
    connect(m_Controls.btSelectWorkingDir, SIGNAL(clicked()), this, SLOT(DoSelectWorkingDir()));
    connect(m_Controls.btResampleNonnodule, SIGNAL(clicked()), this, SLOT(DoResampleNonNodule()));
    connect(m_Controls.btResampleAll, SIGNAL(clicked()), this, SLOT(DoResampleAll()));
    connect(m_Controls.btnResampleOriginal, SIGNAL(clicked()), this, SLOT(DoResampleGroundAndTexture()));
    connect(m_Controls.btTransferToExtImage, SIGNAL(clicked()), this, SLOT(DoExternalTransform()));
    connect(m_Controls.btSelectAsExternalNodule, SIGNAL(clicked()), this, SLOT(DoSelectExternalNodule()));
    connect(m_Controls.btOutputNonExternalNodule, SIGNAL(clicked()), this, SLOT(DoOutputNonExternalNodule()));
    connect(m_Controls.btOutputBOWTraining, SIGNAL(clicked()), this, SLOT(DoOutputBOWTraining()));
  
    connect(m_Controls.btLoadWorkingDataSetConfFile, SIGNAL(clicked()), this, SLOT(LoadWorkingDataSetConfFile()));
    connect(m_Controls.btLoadCurrentWorkingData, SIGNAL(clicked()), this, SLOT(LoadCurrentWorkingListData()));
    connect(m_Controls.btLoadNextWorkingData, SIGNAL(clicked()), this, SLOT(LoadNextWorkingListData()));

    /// inti table view
    m_ReportModel = new QStandardItemModel();
    m_ReportModel->setColumnCount(4);
    m_ReportModel->setHeaderData(0, Qt::Horizontal, QString::fromLocal8Bit("ExpertId"));
    m_ReportModel->setHeaderData(1, Qt::Horizontal, QString::fromLocal8Bit("NoduleId"));
    m_ReportModel->setHeaderData(2, Qt::Horizontal, QString::fromLocal8Bit("X"));
    m_ReportModel->setHeaderData(3, Qt::Horizontal, QString::fromLocal8Bit("Y"));

    m_Controls.tableviewAllInfor->setModel(m_ReportModel);
    m_Controls.tableviewAllInfor->horizontalHeader()->setDefaultAlignment(Qt::AlignLeft);
    m_Controls.tableviewAllInfor->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_Controls.tableviewAllInfor->setEditTriggers(QAbstractItemView::NoEditTriggers);

    if (mitk::IRenderWindowPart* renderWindowPart = GetRenderWindowPart())
    {
        // let the point set widget know about the render window part (crosshair updates)
        RenderWindowPartActivated(renderWindowPart);
    }

    // setup label
    m_Controls.label->setAlignment(Qt::AlignTop);
    m_Controls.label->setWordWrap(true);

    mitk::NodePredicateDimension::Pointer dimensionPredicate = mitk::NodePredicateDimension::New(3);
    mitk::NodePredicateDataType::Pointer imagePredicate = mitk::NodePredicateDataType::New("Image");
    m_Controls.comboGround->SetDataStorage(GetDataStorage());
    m_Controls.comboGround->SetPredicate(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate));
    m_Controls.comboTexture->SetDataStorage(GetDataStorage());
    m_Controls.comboTexture->SetPredicate(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate));
    m_Controls.comboNonnodule->SetDataStorage(GetDataStorage());
    m_Controls.comboNonnodule->SetPredicate(mitk::NodePredicateAnd::New(dimensionPredicate, imagePredicate));

    m_Controls.gridLayout->setHorizontalSpacing(2);
    m_Controls.gridLayout->setVerticalSpacing(2);
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::RenderWindowPartActivated(mitk::IRenderWindowPart* renderWindowPart)
{
    if (this->m_IRenderWindowPart != renderWindowPart)
    {
        this->m_IRenderWindowPart = renderWindowPart;
    }
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::RenderWindowPartDeactivated(mitk::IRenderWindowPart* renderWindowPart)
{
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/,
    const QList<mitk::DataNode::Pointer>& nodes)
{
    // iterate all selected objects, adjust warning visibility
    foreach(mitk::DataNode::Pointer node, nodes)
    {
        if (node.IsNotNull() && dynamic_cast<mitk::Image*>(node->GetData()))
        {
            m_Controls.btLoadAnnotation->setEnabled(true);
            return;
        }
    }

    m_Controls.btLoadAnnotation->setEnabled(false);
}

//////////////////////////////////////////////////////////////////////////
mitk::Image::Pointer LIDCXMLView::getSelectedImage()
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
mitk::DataNode::Pointer LIDCXMLView::getSelectedNode()
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
void LIDCXMLView::fillDataGrid(CNoduleInfoPerSlice* pNoduleInSlice)
{
    // refresh the label
    QString qStrLabel;
    qStrLabel.sprintf("Slice Position: %d; ", pNoduleInSlice->getSlicePosition());

    // clear the model
    m_ReportModel->removeRows(0, m_ReportModel->rowCount());

    // fill grid
    int nTotalItem = 0;
    for (int i = 0; i < N_EXPERT_COUNT; i++)
    {
        int szNoduleInfo = 0;
        SNoduleInfo* pNoduleInfo = pNoduleInSlice->getNoduleInfoThisExpert(i, szNoduleInfo);
        if (NULL == pNoduleInfo)
            continue;

        for (int j = 0; j < szNoduleInfo; j++)
        {
            SNoduleInfo& sNoduleInfo = pNoduleInfo[j];

            // refresh nodule info on label
            qStrLabel += "\n";
            QString qDebugInfo;
            qDebugInfo.sprintf("ExpertIdx: %d, NoduleIdx: %d, Center(%d, %d), Size(%d, %d);", i, j,
                (sNoduleInfo.sBounds.iXMax + sNoduleInfo.sBounds.iXMin) / 2, (sNoduleInfo.sBounds.iYMax + sNoduleInfo.sBounds.iYMin) / 2,
                sNoduleInfo.sBounds.iXMax - sNoduleInfo.sBounds.iXMin, sNoduleInfo.sBounds.iYMax - sNoduleInfo.sBounds.iYMin);
            qStrLabel += qDebugInfo;

            V2IPoints& vPoints = pNoduleInfo[j].vContour;
            int nContour = vPoints.size();
            for (int k = 0; k < nContour; k++)
            {
                m_ReportModel->setItem(nTotalItem, 0, new QStandardItem(QString::number(i)));
                m_ReportModel->setItem(nTotalItem, 1, new QStandardItem(QString::number(j)));
                m_ReportModel->setItem(nTotalItem, 2, new QStandardItem(QString::number(vPoints[k].x)));
                m_ReportModel->setItem(nTotalItem, 3, new QStandardItem(QString::number(vPoints[k].y)));

                nTotalItem++;
            }
        }
    }
    m_Controls.label->setText(qStrLabel);


}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoColorLabel()
{
    QString strFormat;

    destroyContainers();
    
    //-------------------------------------------------------------------------
    // automatic assign node to combos
    mitk::DataNode::Pointer mitkNode1 = this->GetDataStorage()->GetNamedNode("normalized_texture");
    mitk::DataNode::Pointer mitkNode2 = this->GetDataStorage()->GetNamedNode("normalized_groundtruth");
    mitk::DataNode::Pointer mitkNode3 = this->GetDataStorage()->GetNamedNode("normalized_nonnodule");

    if (NULL == mitkNode1 || NULL == mitkNode2)
    {
        MITK_ERROR << "Empty texture or ground truth";
        return;
    }
    else
    {
        m_Controls.comboTexture->SetSelectedNode(mitkNode1);
        m_Controls.comboGround->SetSelectedNode(mitkNode2);
    }

    if (NULL == mitkNode3)
    { 
        MITK_WARN << "NO Nonnodule !! ";
        m_Controls.cbNononnodule->setChecked(true);
    }
    else
    {
        m_Controls.cbNononnodule->setChecked(false);
        m_Controls.comboNonnodule->SetSelectedNode(mitkNode3);
    }

    mitk::DataNode::Pointer mitkNode = m_Controls.comboGround->GetSelectedNode();
    mitk::Image::Pointer	mitkImage = dynamic_cast<mitk::Image*>(mitkNode->GetData());

    //-------------------------------------------------------------------------
    // label nodules
    if (NULL != mitkImage)
    {
        ImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage->Clone(), itkImage);

        //  cast image to binary
        typedef itk::BinaryThresholdImageFilter<ImageType, ImageType>		ThresholdFilterType;
        ThresholdFilterType::Pointer thFilter = ThresholdFilterType::New();
        thFilter->SetInput(itkImage);
        thFilter->SetLowerThreshold(1);
        thFilter->SetUpperThreshold(255);
        thFilter->SetInsideValue(1);
        thFilter->SetOutsideValue(0);
        thFilter->Update();

        typedef itk::MultiplyImageFilter<ImageType, ImageType> MultiplyCastType;
        MultiplyCastType::Pointer multiplyFilter = MultiplyCastType::New();
        multiplyFilter->SetInput(thFilter->GetOutput());
        multiplyFilter->SetConstant(255);
        multiplyFilter->Update();

        typedef itk::CastImageFilter< ImageType, BinaryImageType > ImageCastFilter;
        ImageCastFilter::Pointer castFilter = ImageCastFilter::New();
        castFilter->SetInput(multiplyFilter->GetOutput());
        castFilter->Update();

        typedef itk::BinaryImageToShapeLabelMapFilter<BinaryImageType> BinaryImageToShapeLabelMapFilterType;
        BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
        binaryImageToShapeLabelMapFilter->SetInput(castFilter->GetOutput());
        binaryImageToShapeLabelMapFilter->FullyConnectedOn();
        binaryImageToShapeLabelMapFilter->Update();

        // Loop over all of the blobs
        int szObjects = binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
        strFormat.sprintf("Total %d Nodules.", szObjects);

        MITK_INFO << strFormat.toStdString();
        m_Controls.leTab2Info->setText(strFormat);

        for (unsigned int i = 0; i < szObjects; i++)
        {
            BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
            ImageType::RegionType	region = labelObject->GetBoundingBox();
            m_vNoduleRegions.push_back(region);

            MITK_INFO << "Object " << i << region;
        }

        // Create a label image
        typedef itk::Image<unsigned char, VDIM>  LabelImageType;
        typedef itk::LabelMapToLabelImageFilter<BinaryImageToShapeLabelMapFilterType::OutputImageType, LabelImageType> LabelMapToLabelImageFilterType;
        LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
        labelMapToLabelImageFilter->SetInput(/*shapeOpeningLabelMapFilter*/binaryImageToShapeLabelMapFilter->GetOutput());
        labelMapToLabelImageFilter->Update();

        typedef itk::RGBPixel<unsigned char>		RGBPixelType;
        typedef itk::Image<RGBPixelType, VDIM>		RGBImageType;
        // Color each label/object a different color
        typedef itk::LabelToRGBImageFilter<LabelImageType, RGBImageType> RGBFilterType;
        RGBFilterType::Pointer colormapImageFilter = RGBFilterType::New();
        colormapImageFilter->SetInput(labelMapToLabelImageFilter->GetOutput());
        // 		colormapImageFilter->SetColormap(RGBFilterType::Jet);
        colormapImageFilter->Update();
        RGBImageType::Pointer rgbImage = colormapImageFilter->GetOutput();

        QString strFormat;
        strFormat.sprintf("%s_color", mitkNode->GetName().c_str());

        mitk::Image::Pointer mitkNewColorImage = mitk::Image::New();
        mitk::CastToMitkImage(rgbImage, mitkNewColorImage);
        ShowProcessedDataNode(mitkNewColorImage, strFormat.toStdString());

        // generate thumbnails
        SAFE_DELETE(m_outputNoduleTypes);
        m_outputNoduleTypes = new int[szObjects];
        memset(m_outputNoduleTypes, 0, sizeof(int)*szObjects);

        // assign prior values
        QString q_bow_output_filename;
        q_bow_output_filename.sprintf("%s/bow/training/types_config.txt", m_WorkingDir.toStdString().c_str());
        MITK_INFO << "reading config from " << q_bow_output_filename.toStdString();

        QFile qBOWInput(q_bow_output_filename);
        QTextStream tsBOWInput(&qBOWInput);
        if (qBOWInput.open(QIODevice::ReadOnly)){
            int nCount = 0;
            tsBOWInput >> nCount;
            if (nCount == szObjects){
                for (int iNodule = 0; iNodule < szObjects; iNodule++){
                    tsBOWInput >> m_outputNoduleTypes[iNodule];
                }
            }
            else{
                MITK_WARN << "prior type configs not equal with this one";
            }
            qBOWInput.close();
        }
        else{
            MITK_ERROR << "Open " << q_bow_output_filename.toStdString() << " Failed";
        }


        // enumerate all objects
        for (unsigned int i = 0; i < szObjects; i++)
        {
            BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
            ImageType::RegionType	region = labelObject->GetBoundingBox();
            ImageType::IndexType    region_index = region.GetIndex();
            ImageType::SizeType     region_size = region.GetSize();

            ImageType::IndexType center_index;
            center_index[0] = region_index[0] + region_size[0] / 2;
            center_index[1] = region_index[1] + region_size[1] / 2;
            center_index[2] = region_index[2] + region_size[2] / 2;

            RGBPixelType rgbPixel = rgbImage->GetPixel(center_index);
            QImage qThumbnail(16, 16, QImage::Format_RGB32);
            qThumbnail.fill(QColor(rgbPixel[0], rgbPixel[1], rgbPixel[2]));

            QLabel* qImageLabel = new QLabel();
            qImageLabel->setPixmap(QPixmap::fromImage(qThumbnail));

            QComboBox* qCombo = new QComboBox();
            qCombo->setProperty("nodule_index", i);
            qCombo->addItem("ISOLATED");
            qCombo->addItem("WALL");
            qCombo->addItem("VESSEL");
            qCombo->addItem("GGO");
            qCombo->addItem("TRAIL");
            qCombo->setCurrentIndex(m_outputNoduleTypes[i]);
            connect(qCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(outputNoduleTypeChanged(int)));

            QLabel* qLabel = new QLabel();
            qLabel->setText(QString("Info [%1 %2 %3],[%4 %5 %6]").arg(region_index[0]).arg(region_index[1]).arg(region_index[2]).arg(region_size[1]).arg(region_size[1]).arg(region_size[2]));

            QLabel* labelTmp1 = new QLabel();
            labelTmp1->setText("Color");
            QLabel* labelTmp2 = new QLabel();
            labelTmp2->setText("Type");

            m_Controls.gridLayout->addWidget(labelTmp1, i, 0);
            m_Controls.gridLayout->addWidget(qImageLabel, i, 1);
            m_Controls.gridLayout->addWidget(labelTmp2, i, 2);
            m_Controls.gridLayout->addWidget(qCombo, i, 3);
            m_Controls.gridLayout->addWidget(qLabel, i, 4);
        }
    }

    //-------------------------------------------------------------------------
    // label non-nodules
    mitk::DataNode::Pointer mitkNonNoduleNode = m_Controls.comboNonnodule->GetSelectedNode();
    mitk::Image::Pointer	mitkNonNoduleImage = dynamic_cast<mitk::Image*>(mitkNonNoduleNode->GetData());
    bool bNononnoduleChecked = m_Controls.cbNononnodule->isChecked();

    if (NULL != mitkNonNoduleImage && !bNononnoduleChecked)
    {
        ImageType::Pointer itkNonNoduleImage;
        mitk::CastToItkImage(mitkNonNoduleImage->Clone(), itkNonNoduleImage);

        //  cast image to binary
        typedef itk::BinaryThresholdImageFilter<ImageType, ImageType>		ThresholdFilterType;
        ThresholdFilterType::Pointer thFilter = ThresholdFilterType::New();
        thFilter->SetInput(itkNonNoduleImage);
        thFilter->SetLowerThreshold(1);
        thFilter->SetUpperThreshold(255);
        thFilter->SetInsideValue(1);
        thFilter->SetOutsideValue(0);
        thFilter->Update();

        typedef itk::MultiplyImageFilter<ImageType, ImageType> MultiplyCastType;
        MultiplyCastType::Pointer multiplyFilter = MultiplyCastType::New();
        multiplyFilter->SetInput(thFilter->GetOutput());
        multiplyFilter->SetConstant(255);
        multiplyFilter->Update();

        typedef itk::CastImageFilter< ImageType, BinaryImageType > ImageCastFilter;
        ImageCastFilter::Pointer castFilter = ImageCastFilter::New();
        castFilter->SetInput(multiplyFilter->GetOutput());
        castFilter->Update();

        typedef itk::BinaryImageToShapeLabelMapFilter<BinaryImageType> BinaryImageToShapeLabelMapFilterType;
        BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
        binaryImageToShapeLabelMapFilter->SetInput(castFilter->GetOutput());
        binaryImageToShapeLabelMapFilter->FullyConnectedOn();
        binaryImageToShapeLabelMapFilter->Update();

        // Loop over all of the blobs
        int szObjects = binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
        strFormat.sprintf("Total %d Non-nodules.", szObjects);
        MITK_INFO << strFormat.toStdString();
        QString dbgText = m_Controls.leTab2Info->text();
        dbgText += "\n";
        dbgText += strFormat;
        m_Controls.leTab2Info->setText(dbgText);

        for (unsigned int i = 0; i < szObjects; i++)
        {
            BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
            ImageType::RegionType	region = labelObject->GetBoundingBox();
            //ImageType::IndexType	region_index = region.GetIndex();
            //ImageType::SizeType     region_size = region.GetSize();
            m_NonNoduleInfoFromLabel.push_back(region);

            MITK_INFO << "Non-nodule " << i << region;
        }
    }
}

////////////////////////////////////////////////////////////////////////////
//void LIDCXMLView::saveImages()
//{
//#ifdef _OUTPUT_2D_
//    QFileDialog* openFilePath = new QFileDialog(NULL, "Select Folder", "file");
//    openFilePath->setFileMode(QFileDialog::DirectoryOnly);
//    QString dirName = openFilePath->getExistingDirectory();
//    if ("" != dirName)
//    {
//        DoResampleGroundAndTexture(dirName);
//    }
//    delete openFilePath;
//#else
//    DoResampleGroundAndTexture();
//#endif
//
//}
//
//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::outputNoduleTypeChanged(int iType)
{
    QString dbg_str;

    QComboBox* comboxBox = qobject_cast<QComboBox*>(QObject::sender());
    QVariant qNoduleIndex = comboxBox->property("nodule_index");
    int iNodule = qNoduleIndex.toInt();

    if (iType >= EOutType_COUNT){
        dbg_str.sprintf("Non-valid nodule type, Index: %d, Type: %d", iNodule, iType);
        MITK_WARN << dbg_str.toStdString();
        return;
    }

    m_outputNoduleTypes[iNodule] = iType;
    dbg_str.sprintf("Changing Index: %d, Type: %d", iNodule, iType);
    MITK_INFO << dbg_str.toStdString();
}

//////////////////////////////////////////////////////////////////////////
void iteratorTest()
{
#if 0
    ImageType::SizeType regionSize;
    regionSize[0] = 1;
    regionSize[1] = 1;
    regionSize[2] = 1;

    ImageType::IndexType regionIndex;
    regionIndex[0] = 0;
    regionIndex[1] = 0;
    regionIndex[2] = 0;

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

    ImageType::Pointer itkTexture = ImageType::New();
    mitk::CastToItkImage(mitkTexture, itkTexture);
    itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
    int nTotalOut = 0;
    // one piece of training data
    while (!iterator.IsAtEnd())
    {
        ImageType::IndexType centerIndex = iterator.GetIndex();
        for (unsigned int i = 0; i < iterator.Size(); i++)
        {
            ImageType::IndexType index = iterator.GetIndex(i);

            bool IsInBounds;
            int neighborValue = iterator.GetPixel(i, IsInBounds);
            if (IsInBounds)
            {
                nTotalOut++;
                strQDebug.sprintf("[%d, %d, %d]", index[0], index[1], index[2]);
                MITK_INFO << strQDebug;
            }
        }
        ++iterator;
    }
    MITK_INFO << "Total" << nTotalOut;
#endif

}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::buildNormalizedGaussian(double***&dwRET, int nGaussianR)
{
    int nGaussianD = 2 * nGaussianR + 1;

    double dwGaussianSigma = 0.0;
    double dwGaussianMax = 0.0;
    double dwGaussianMin = 1.0;
    for (int local_x_idx = 0; local_x_idx < nGaussianD; local_x_idx++){
        for (int local_y_idx = 0; local_y_idx < nGaussianD; local_y_idx++){
            for (int local_z_idx = 0; local_z_idx < nGaussianD; local_z_idx++){
                int local_x = (local_x_idx - nGaussianR);
                int local_y = (local_y_idx - nGaussianR);
                int local_z = (local_z_idx - nGaussianR);
                double dwThisVal = /*DW_2PI_SIGMA2*/SQRT_DW_2PI_SIGMA2 * exp(-(local_x*local_x + local_y*local_y + local_z*local_z) / (2 * DW_SIGMA*DW_SIGMA));
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
    for (int local_x_idx = 0; local_x_idx < nGaussianD; local_x_idx++){
        for (int local_y_idx = 0; local_y_idx < nGaussianD; local_y_idx++){
            for (int local_z_idx = 0; local_z_idx < nGaussianD; local_z_idx++){
                dwRET[local_x_idx][local_y_idx][local_z_idx] = (dwRET[local_x_idx][local_y_idx][local_z_idx] - dwGaussianMin) / (dwGaussianMax - dwGaussianMin);
            }
        }
    }
    // 	dwGaussianMax /= dwGaussianSigma;
    // 	dwGaussianMax /= dwGaussianMax;
    // 	return dwGaussianMax;
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::outputNoduleTrainingData(QString strSavePath, mitk::Image::Pointer mitkTexture, mitk::Image::Pointer mitkGroundTruth)
{
    // for debugging formation
    QString strQDebug;

    mitk::Image::Pointer mitkTMP = mitkGroundTruth->Clone();
    DoubleImageType::Pointer itkDoubleImage;
    mitk::CastToItkImage(mitkTMP, itkDoubleImage);
    itkDoubleImage->FillBuffer(0.0);
    mitk::Image::Pointer mitkDoubleImage;
    mitk::CastToMitkImage(itkDoubleImage, mitkDoubleImage);

    int iWidth = mitkTexture->GetDimension(0);
    int iHeight = mitkTexture->GetDimension(1);
    int iSlice = mitkTexture->GetDimension(2);

    // file header: sz_data szInput szOutput
    //     int szData = N_NODULE_D * N_NODULE_D * N_NODULE_D;
    int szInput = N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D;
    int szOutput = 1;

    if (NULL == m_outputNoduleTypes){
        MITK_WARN << "Output Types Empty...";
        return;
    }

    // collect count of each type nodule
    int szNodulesThisType[EOutType_COUNT] = { 0 };
    int szAllNodules = m_vNoduleRegions.size();
    for (int iNodule = 0; iNodule < szAllNodules; iNodule++)
    {
        ImageType::RegionType& nodule_region = m_vNoduleRegions[iNodule];
        ImageType::SizeType nodule_size = nodule_region.GetSize();
        int eOutputNoduleType = m_outputNoduleTypes[iNodule];

        int iRegionX = min((int)nodule_size[0], N_NODULE_D);
        int iRegionY = min((int)nodule_size[1], N_NODULE_D);
        int iRegionZ = min((int)nodule_size[2], N_NODULE_D);

        szNodulesThisType[eOutputNoduleType] += (iRegionX*iRegionY*iRegionZ);
    }

    // file handling
    // texture 
    QString qSTRSuffix[EOutType_COUNT] = {
        "annTrainingData_ISO.txt",
        "annTrainingData_WALL.txt",
        "annTrainingData_VESSEL.txt",
        "annTrainingData_GGO.txt",
    };
    QFile qTraining[EOutType_COUNT];
    QTextStream tsTraining[EOutType_COUNT];
    for (int i = 0; i < EOutType_COUNT; i++)
    {
        QString tmpPath = strSavePath + qSTRSuffix[i];
        qTraining[i].setFileName(tmpPath);
        MITK_INFO << "Writing to " << tmpPath;

        if (qTraining[i].open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
        {
            tsTraining[i].setDevice(&qTraining[i]);
            tsTraining[i].setRealNumberNotation(QTextStream::FixedNotation);
            // 1/2000 = 0.0005, so 4 digits after . is enough
            tsTraining[i].setRealNumberPrecision(PIXEL_PRECISION);
            tsTraining[i] << szNodulesThisType[i] << "\t" << szInput << "\t" << szOutput << "\n";
        }
        else
        {
            MITK_INFO << "Open" << tmpPath << " Failed";
            return;
        }
    }

    // transform the mitk texture
    ImageType::Pointer itkTexture;
    mitk::CastToItkImage(mitkTexture, itkTexture);

    for (int iNodule = 0; iNodule < szAllNodules; iNodule++)
    {
        ImageType::RegionType& nodule_region = m_vNoduleRegions[iNodule];
        ImageType::IndexType nodule_index = nodule_region.GetIndex();
        ImageType::SizeType nodule_size = nodule_region.GetSize();
        int eOutputNoduleType = m_outputNoduleTypes[iNodule];

        int iRegionX = min((int)nodule_size[0], N_NODULE_D);
        int iRegionY = min((int)nodule_size[1], N_NODULE_D);
        int iRegionZ = min((int)nodule_size[2], N_NODULE_D);

        ImageType::IndexType center_index;
        center_index[0] = nodule_index[0] + nodule_size[0] / 2;
        center_index[1] = nodule_index[1] + nodule_size[1] / 2;
        center_index[2] = nodule_index[2] + nodule_size[2] / 2;

        int iMaxDim = max(iRegionX, iRegionY);
        iMaxDim = max(iMaxDim, iRegionZ);

        int iGaussianR = (iMaxDim + 0.5) / 2;
        int iGaussianD = iGaussianR * 2 + 1;

        strQDebug.sprintf("working on nodule %d, type %d, center(%d,%d,%d), GaussD %d, Size(%d, %d, %d)",
            iNodule, eOutputNoduleType, center_index[0], center_index[1], center_index[2], iGaussianD, iRegionX, iRegionY, iRegionZ);
        MITK_INFO << strQDebug.toStdString();

        // build 3D gaussian kernel
        double ***dwGaussian = (double***)malloc(iGaussianD * sizeof(double**));
        for (int i = 0; i < iGaussianD; i++){
            dwGaussian[i] = (double**)malloc(iGaussianD * sizeof(double*));
            for (int j = 0; j < iGaussianD; j++){
                dwGaussian[i][j] = (double*)malloc(iGaussianD * sizeof(double));
            }
        }
        buildNormalizedGaussian(dwGaussian, iGaussianR);

        // do through mitk image iteration
        mitk::ImagePixelReadAccessor<TPixelType, 3>			ground_access(mitkGroundTruth);
        mitk::ImagePixelWriteAccessor<DoublePixelType, 3>	double_image_access(mitkDoubleImage);
        TPixelType center_value = ground_access.GetPixelByIndexSafe(/*access_index*/center_index);
        //         if (center_value > 0)
        {
            //DoubleImageType::IndexType	access_index;
            ImageType::SizeType			regionSize;
            ImageType::IndexType		regionIndex;
            ImageType::RegionType		region;
            ImageType::SizeType			radius;

            regionSize[0] = iRegionX;
            regionSize[1] = iRegionY;
            regionSize[2] = iRegionZ;
            regionIndex[0] = center_index[0] - iRegionX / 2 - 1;
            regionIndex[1] = center_index[1] - iRegionY / 2 - 1;
            regionIndex[2] = center_index[2] - iRegionZ / 2 - 1;

            if (regionIndex[0] < 0){
                regionIndex[0] = 0;
                MITK_WARN << "Iterating Region Index0 < 0";
            }
            if (regionIndex[0] > iWidth){
                regionIndex[0] = iWidth;
                MITK_WARN << "Iterating Region Index0 > Width";
            }
            if (regionIndex[1] < 0){
                regionIndex[1] = 0;
                MITK_WARN << "Iterating Region Index1 < 0";
            }
            if (regionIndex[1] > iHeight){
                regionIndex[1] = iHeight;
                MITK_WARN << "Iterating Region Index1 > Height";
            }
            if (regionIndex[2] < 0){
                regionIndex[2] = 0;
                MITK_WARN << "Iterating Region Index2 < 0";
            }
            if (regionIndex[2] > iSlice){
                regionIndex[2] = iSlice;
                MITK_WARN << "Iterating Region Index2 > Slice";
            }

            // region size
            region.SetSize(regionSize);
            // region start
            region.SetIndex(regionIndex);
            radius[0] = N_LOCAL_WINDOW_R;
            radius[1] = N_LOCAL_WINDOW_R;
            radius[2] = N_LOCAL_WINDOW_R;

            itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
            int nCounter = 0;
            while (!iterator.IsAtEnd())
            {
                ImageType::IndexType local_centerIndex = iterator.GetIndex();
                int nTotalIn = 0;
                for (unsigned int iITER = 0; iITER < iterator.Size(); ++iITER)
                {
                    ImageType::IndexType neighbour_index = iterator.GetIndex(iITER);

                    bool IsInBounds;
                    int neighborValue = iterator.GetPixel(iITER, IsInBounds);
                    if (!IsInBounds){
                        neighborValue = -F_CT_WIDTH;
                    }
                    if (neighborValue > F_CT_WIDTH){
                        neighborValue = F_CT_WIDTH;
                    }
                    if (neighborValue < -F_CT_WIDTH){
                        neighborValue = -F_CT_WIDTH;
                    }

                    double dwNormalized = ((double)(neighborValue + F_CT_WIDTH)) / ((double)(F_CT_WIDTH * 2));
                    // strQDebug.sprintf("[%d %d %d, %d]-%f", i, j, k, IsInBounds ? 1 : 0, dwNormalized);
                    // tsTraining << strQDebug << "\t";
                    tsTraining[eOutputNoduleType] << /*neighborValue*/dwNormalized << "\t";
                    nTotalIn++;
                }
                tsTraining[eOutputNoduleType] << "\n";

                TPixelType local_center_value = ground_access.GetPixelByIndexSafe(local_centerIndex);
                int local_x = local_centerIndex[0] - center_index[0] + iGaussianR;
                int local_y = local_centerIndex[1] - center_index[1] + iGaussianR;
                int local_z = local_centerIndex[2] - center_index[2] + iGaussianR;

                bool bvalidInGaussian = (local_x >= 0 && local_x < iGaussianD && local_y >= 0 && local_y < iGaussianD && local_z >= 0 && local_z < iGaussianD) ? true : false;

                if (local_center_value > 0 && bvalidInGaussian){
                    double_image_access.SetPixelByIndex(local_centerIndex, dwGaussian[local_x][local_y][local_z] * 255);

                    // to output high level result
                    QString qSTRGaussian = QString::number(dwGaussian[local_x][local_y][local_z], 'f', RESULT_PRECISION);
                    tsTraining[eOutputNoduleType] << qSTRGaussian;
                    //tsTraining << 1;
                    //tsTraining << "[" << nCounter << ", " << iGaussianR << ", " << local_x << ", " << local_y << ", " << local_z << ", " << "]" << " ---> " << dwGaussian[local_x][local_y][local_z];
                }
                else{

                    // labeled as non-nodules
                    double dwOutput = (double)MyMath::Random(0, 10000) / pow(10.0, RESULT_PRECISION) + 0.000000000001;
                    QString qSTRDWResult;
                    qSTRDWResult = QString::number(dwOutput, 'f', RESULT_PRECISION);

                    double_image_access.SetPixelByIndex(local_centerIndex, 0);
                    tsTraining[eOutputNoduleType] << qSTRDWResult;
                    // 						//tsTraining << "[" << nCounter << ", " << iGaussianR << ", " << local_x << ", " << local_y << ", " << local_z << "]" << " ---> " << 0;
                }


                tsTraining[eOutputNoduleType] << "\n";

                if (nTotalIn != szInput)
                {
                    strQDebug.sprintf("number of output training data not equal with sample size, Needed: %d, Actual: %d", szInput, nTotalIn);
                    MITK_WARN << strQDebug.toStdString();
                }

                // do not forget move to next data
                ++iterator;
                nCounter++;
            }
            //             strQDebug.sprintf("nodule %d, counter %d", iNodule, nCounter);
            //             MITK_INFO << strQDebug.toStdString();

            if (nCounter != iRegionX*iRegionY*iRegionZ)
            {
                strQDebug.sprintf("Nodule %d output not match, Needed: %d, Actual: %d", iNodule, iRegionX*iRegionY*iRegionZ, nCounter);
                MITK_WARN << strQDebug.toStdString();
            }
        }
        for (int i = 0; i < iGaussianD; i++){
            for (int j = 0; j < iGaussianD; j++){
                SAFE_DELETE(dwGaussian[i][j]);
            }
            SAFE_DELETE(dwGaussian[i]);
        }
        SAFE_DELETE(dwGaussian);
    }

    for (int i = 0; i < EOutType_COUNT; i++){
        qTraining[i].close();
    }


    strQDebug.sprintf("double_image");
    ShowProcessedDataNode(mitkDoubleImage, strQDebug.toStdString());
    // 	}
    // 	else
    // 	{
    // 		MITK_WARN << "Open " << strSavePath << " Failed...";
    // 	}
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::outputNonNoduleTrainingData(QString strSavePath, mitk::Image::Pointer mitkTexture, mitk::Image::Pointer mitkNonNodule)
{
    QString strQDebug;

    QString qSTRSuffix = "annTrainingData_NN.txt";
    QString tmpPath = strSavePath + qSTRSuffix;

    QFile qTraining(tmpPath);
    QTextStream tsTraining(&qTraining);
    MITK_INFO << "Writing to " << tmpPath;

    int iWidth = mitkTexture->GetDimension(0);
    int iHeight = mitkTexture->GetDimension(1);
    int iSlice = mitkTexture->GetDimension(2);

    ImageType::Pointer itkTexture;
    mitk::CastToItkImage(mitkTexture, itkTexture);

    ImageType::Pointer itkNonNodule;
    mitk::CastToItkImage(mitkNonNodule, itkNonNodule);


    // file header: sz_data szInput szOutput
    int szAllNonNodulesFromLabel = m_NonNoduleInfoFromLabel.size();
    int szNonNodulesFromXML = m_NonNoduleInfoFromXML.size();
    int szData = N_NON_NODULE_D * N_NON_NODULE_D * N_NON_NODULE_D;
    int szInput = N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D;
    int szOutput = 1;

    if (qTraining.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
    {
        tsTraining.setRealNumberNotation(QTextStream::FixedNotation);
        tsTraining.setRealNumberPrecision(PIXEL_PRECISION);
        if (szAllNonNodulesFromLabel == 0)
        {
            tsTraining << 0 << "\t" << szInput << "\t" << szOutput << "\n";
            qTraining.close();
            return;
        }

        tsTraining << szAllNonNodulesFromLabel * szData << "\t" << szInput << "\t" << szOutput << "\n";

        for (int iNonNodule = 0; iNonNodule < szAllNonNodulesFromLabel; iNonNodule++)
        {
            ImageType::RegionType non_nodule_region = m_NonNoduleInfoFromLabel[iNonNodule];
            ImageType::IndexType non_nodule_region_index = non_nodule_region.GetIndex();
            ImageType::SizeType non_nodule_region_size = non_nodule_region.GetSize();

            ImageType::IndexType center_index;
            center_index[0] = non_nodule_region_index[0] + non_nodule_region_size[0] / 2;
            center_index[1] = non_nodule_region_index[1] + non_nodule_region_size[1] / 2;
            center_index[2] = non_nodule_region_index[2] + non_nodule_region_size[2] / 2;

            strQDebug.sprintf("processing non-nodule %d, center(%d,%d,%d), ", iNonNodule, center_index[0], center_index[1], center_index[2]);
            MITK_INFO << strQDebug.toStdString();

            //DoubleImageType::IndexType	access_index;
            ImageType::SizeType			regionSize;
            ImageType::IndexType		regionIndex;
            ImageType::RegionType		region;
            ImageType::SizeType			radius;

            regionSize[0] = N_NON_NODULE_D;
            regionSize[1] = N_NON_NODULE_D;
            regionSize[2] = N_NON_NODULE_D;
            regionIndex[0] = center_index[0] - N_NON_NODULE_R;
            regionIndex[1] = center_index[1] - N_NON_NODULE_R;
            regionIndex[2] = center_index[2] - N_NON_NODULE_R;

            if (regionIndex[0] < 0){
                regionIndex[0] = 0;
                MITK_WARN << "Iterating Region Index0 < 0";
            }
            if (regionIndex[0] > iWidth){
                regionIndex[0] = iWidth;
                MITK_WARN << "Iterating Region Index0 > Width";
            }
            if (regionIndex[1] < 0){
                regionIndex[1] = 0;
                MITK_WARN << "Iterating Region Index1 < 0";
            }
            if (regionIndex[1] > iHeight){
                regionIndex[1] = iHeight;
                MITK_WARN << "Iterating Region Index1 > Height";
            }
            if (regionIndex[2] < 0){
                regionIndex[2] = 0;
                MITK_WARN << "Iterating Region Index2 < 0";
            }
            if (regionIndex[2] > iSlice){
                regionIndex[2] = iSlice;
                MITK_WARN << "Iterating Region Index2 > Slice";
            }

            // region size
            region.SetSize(regionSize);
            // region start
            region.SetIndex(regionIndex);
            radius[0] = N_LOCAL_WINDOW_R;
            radius[1] = N_LOCAL_WINDOW_R;
            radius[2] = N_LOCAL_WINDOW_R;

            itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
            int nCounter = 0;
            while (!iterator.IsAtEnd())
            {
                ImageType::IndexType local_centerIndex = iterator.GetIndex();
                int nTotalIn = 0;
                for (unsigned int iITER = 0; iITER < iterator.Size(); ++iITER)
                {
                    ImageType::IndexType neighbour_index = iterator.GetIndex(iITER);

                    bool IsInBounds;
                    int neighborValue = iterator.GetPixel(iITER, IsInBounds);
                    if (!IsInBounds){
                        neighborValue = -F_CT_WIDTH;
                    }
                    if (neighborValue > F_CT_WIDTH){
                        neighborValue = F_CT_WIDTH;
                    }
                    if (neighborValue < -F_CT_WIDTH){
                        neighborValue = -F_CT_WIDTH;
                    }

                    double dwNormalized = ((double)(neighborValue + F_CT_WIDTH)) / ((double)(F_CT_WIDTH * 2));
                    // strQDebug.sprintf("[%d %d %d, %d]-%f", i, j, k, IsInBounds ? 1 : 0, dwNormalized);
                    // tsTraining << strQDebug << "\t";
                    tsTraining << /*neighborValue*/dwNormalized << "\t";
                    nTotalIn++;
                }
                tsTraining << "\n";

                // labeled as non-nodules
                double dwOutput = (double)MyMath::Random(0, 10000) / pow(10.0, RESULT_PRECISION) + 0.000000000001;
                QString qSTRDWResult;
                qSTRDWResult = QString::number(dwOutput, 'f', RESULT_PRECISION);

                tsTraining << qSTRDWResult;
                tsTraining << "\n";

                if (nTotalIn != szInput)
                {
                    strQDebug.sprintf("number of output training data not equal with sample size, Needed: %d, Actual: %d", szInput, nTotalIn);
                    MITK_WARN << strQDebug;
                }

                // do not forget move to next data
                ++iterator;
                nCounter++;
            }
        }
        //              for (int i = 0; i < iGaussianD; i++){
        //                  for (int j = 0; j < iGaussianD; j++){
        //                      SAFE_DELETE(dwGaussian[i][j]);
        //                  }
        //                  SAFE_DELETE(dwGaussian[i]);
        //              }
        //              SAFE_DELETE(dwGaussian);
        //         }

        qTraining.close();
    }
    else
    {
        MITK_INFO << "Open" << tmpPath << " Failed";
        return;
    }

}

//////////////////////////////////////////////////////////////////////////
// void LIDCXMLView::outputANNTrainingData(QString strSavePath, mitk::Image* mitkTexture, mitk::Image* mitkGroundTruth)
// {
// 	// #define _DBG_OUTPUT_TRAINING_
// 
// 	// for debugging formation
// 	QString strQDebug;
// 
// 	int iWidth = mitkTexture->GetDimension(0);
// 	int iHeight = mitkTexture->GetDimension(1);
// 	int iSlice = mitkTexture->GetDimension(2);
// 
// 	// file header: sz_data szInput szOutput
// 	int szData = N_NODULE_D * N_NODULE_D * N_NODULE_D;
// 	int szInput = N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D;
// 	int szOutput = 1;
// 	// texture 
// 	QString qSTRSuffix = "annTrainingData.txt";
// 	strSavePath += qSTRSuffix;
// 	// open file writer pointers
// 	QFile qTraining(strSavePath);
// 
// 	if (qTraining.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
// 	{
// 		QTextStream tsTraining(&qTraining);
// 
// 		int szAllNodules = 0;
// 		int szNodules = m_NodulesInSlice.size();
// 		for (int i = 0; i < szNodules; i++)
// 		{
// 			CNoduleInfoPerSlice* pNodulesInSlice = m_NodulesInSlice[i];
// 			for (int j = 0; j < N_EXPERT_COUNT; j++)
// 			{// towards each expert on this slice 
// 				int szThisSlice = 0;
// 				SNoduleInfo* pInfoArray = pNodulesInSlice->getNoduleInfoThisExpert(j, szThisSlice);
// 				szAllNodules += szThisSlice;
// 			}
// 		}
// 		MITK_WARN << "All Nodule Size: " << szAllNodules;
// 
// 		// output params
// 		tsTraining << szAllNodules * szData << "\t" << szInput << "\t" << szOutput << "\n";
// 
// 		// transform the mitk texture
// 		ImageType::Pointer itkTexture;
// 		mitk::CastToItkImage(mitkTexture, itkTexture);
// 
// 		// build 2D gaussian kernel
// 		double* dwGaussian = new double[N_NODULE_D*N_NODULE_D];
// 		for (int local_x_idx = 0; local_x_idx < N_NODULE_D; local_x_idx++){
// 			for (int local_y_idx = 0; local_y_idx < N_NODULE_D; local_y_idx++){
// 				int local_x = (local_x_idx - N_NODULE_R);
// 				int local_y = (local_y_idx - N_NODULE_R);
// 				double dwThisVal = DW_2PI_SIGMA2 * exp(-(local_x*local_x + local_y*local_y) / (2 * DW_SIGMA*DW_SIGMA));
// 				dwGaussian[local_y_idx*N_NODULE_D + local_x_idx] = dwThisVal;
// 			}
// 		}
// 		// normalize 
// 		double dwGaussianSigma = 0.0;
// 		for (int iGaussian = 0; iGaussian < N_NODULE_D*N_NODULE_D; iGaussian++){
// 			dwGaussianSigma += dwGaussian[iGaussian];
// 		}
// 		for (int iGaussian = 0; iGaussian < N_NODULE_D*N_NODULE_D; iGaussian++){
// 			dwGaussian[iGaussian] /= dwGaussianSigma;
// 		}
// 
// 		// output training data
// 		for (int i = 0; i < szNodules; i++)
// 		{// towards all the nodules of this image set
// 			CNoduleInfoPerSlice*				pNodulesInSlice = m_NodulesInSlice[i];
// 			mitk::Image::ImageDataItemPointer	pThisSlice = mitkTexture->GetSliceData(pNodulesInSlice->getSlicePosition());
// 			mitk::Image::ImageDataItemPointer   pGroundThisSlice = mitkGroundTruth->GetSliceData(pNodulesInSlice->getSlicePosition());
// 
// 			// read access for ground truth
// 			mitk::ImagePixelReadAccessor<TPixelType, 2> readAccess(mitkGroundTruth, pGroundThisSlice);
// 			itk::Index<2> read_idx;
// 
// 			for (int j = 0; j < N_EXPERT_COUNT; j++)
// 			{// towards each expert on this slice 
// 				int szThisSlice = 0;
// 				SNoduleInfo* pInfoArray = pNodulesInSlice->getNoduleInfoThisExpert(j, szThisSlice);
// 				if (NULL != pInfoArray)
// 				{
// 					for (int k = 0; k < szThisSlice; k++)
// 					{
// 						SNoduleInfo&	sThisInfo = pInfoArray[k];
// 						SNoduleBounds&	sThisBounds = sThisInfo.sBounds;
// 						int iBoundsWidth = (sThisBounds.iXMax - sThisBounds.iXMin);
// 						int iBoundsHeight = (sThisBounds.iYMax - sThisBounds.iYMin);
// 						int iCenterX = ceil((sThisBounds.iXMin + sThisBounds.iXMax) / 2.0);
// 						int iCenterY = ceil((sThisBounds.iYMin + sThisBounds.iYMax) / 2.0);
// 						// bounds may be out of N_NODULE_D!!
// 						int iCenterBoundX1 = iCenterX - N_NODULE_R;
// 						int iCenterBoundX2 = iCenterX + N_NODULE_R;
// 						int iCenterBoundY1 = iCenterY - N_NODULE_R;
// 						int iCenterBoundY2 = iCenterY + N_NODULE_R;
// 
// 						ImageType::SizeType regionSize;
// 						regionSize[0] = N_NODULE_D/*iBoundsWidth*/;
// 						regionSize[1] = N_NODULE_D/*iBoundsHeight*/;
// 						regionSize[2] = N_NODULE_D;
// 
// 						ImageType::IndexType regionIndex;
// 						regionIndex[0] = iCenterX - N_NODULE_R/*sThisBounds.iXMin*/;
// 						regionIndex[1] = iCenterY - N_NODULE_R/*sThisBounds.iYMin*/;
// 						regionIndex[2] = pNodulesInSlice->getSlicePosition() - N_NODULE_R;
// 
// 						// start at index and look for neighbors at each point inside of index+size;
// 						ImageType::RegionType region;
// 						// region size
// 						region.SetSize(regionSize);
// 						// region start
// 						region.SetIndex(regionIndex);
// 
// 						ImageType::SizeType radius;
// 						radius[0] = N_LOCAL_WINDOW_R;
// 						radius[1] = N_LOCAL_WINDOW_R;
// 						radius[2] = N_LOCAL_WINDOW_R;
// 
// 						itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
// 						// one piece of training data
// 						while (!iterator.IsAtEnd())
// 						{
// 							ImageType::IndexType centerIndex = iterator.GetIndex();
// 							read_idx[0] = centerIndex[0];
// 							read_idx[1] = centerIndex[1];
// 							int centerValue = readAccess.GetPixelByIndex(read_idx);
// 							int nTotalIn = 0;
// 							for (unsigned int i = 0; i < iterator.Size(); i++)
// 							{
// 								ImageType::IndexType index = iterator.GetIndex(i);
// 
// 								bool IsInBounds;
// 								int neighborValue = iterator.GetPixel(i, IsInBounds);
// 								if (!IsInBounds)
// 								{
// 									neighborValue = -(int)F_CT_WIDTH;
// 								}
// #ifdef _DBG_OUTPUT_TRAINING_
// 								strQDebug.sprintf("[%d %d %d, %d]-%d", index[0], index[1], index[2], IsInBounds ? 1 : 0, neighborValue);
// 								tsTraining << strQDebug << "\t";
// #else
// 								tsTraining << neighborValue << "\t";
// #endif //_DBG_OUTPUT_TRAINING_
// 								nTotalIn++;
// 
// 							}
// 							tsTraining << "\n";
// #ifdef _DBG_OUTPUT_TRAINING_
// 							strQDebug.sprintf("[%d %d %d]-%d", centerIndex[0], centerIndex[1], pNodulesInSlice->getSlicePosition(), centerValue);
// 							tsTraining << strQDebug;
// #else
// 
// #ifdef ANN_OUTPUT_GAUSSIAN
// 							int local_center_x = centerIndex[0] - iCenterX;
// 							int local_center_y = centerIndex[1] - iCenterY;
// 
// 							int local_center_x_idx = local_center_x + N_NODULE_R;
// 							int local_center_y_idx = local_center_y + N_NODULE_R;
// 							// 	QString strDBG;
// 							// 	strDBG.sprintf("[%d %d %f]", local_center_x, local_center_y, dwGaussian[local_center_y_idx*N_NODULE_D + local_center_x_idx]);
// 							//	tsTraining << strDBG;
// 							tsTraining << dwGaussian[local_center_y_idx*N_NODULE_D + local_center_x_idx];
// 
// #else
// 							tsTraining << centerValue;
// #endif // ANN_OUTPUT_GAUSSIAN
// #endif //_DBG_OUTPUT_TRAINING_
// 
// 							tsTraining << "\n";
// 
// 							if (nTotalIn != szInput)
// 							{
// 								strQDebug.sprintf("number of output training data not equal with sample size, Needed: %d, Actual: %d", szInput, nTotalIn);
// 								MITK_WARN << strQDebug;
// 							}
// 
// 							// do not forget move to next data
// 							++iterator;
// 						}
// 
// 					}
// 				}
// 			}
// 		}
// 		SAFE_DELETE(dwGaussian);
// 		qTraining.close();
// 	}
// 	else
// 	{
// 		MITK_WARN << "Open " << strSavePath << " Failed...";
// 	}
// }


//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoResampleGroundAndTexture()
{

    QString strDBG;

    MITK_INFO << "DoResampleGroundAndTexture";

    mitk::Image* mitkImage = getSelectedImage();
    mitk::Image::Pointer mitkImageClone = mitkImage->Clone();

    ImageType::Pointer itkImage;
    mitk::CastToItkImage(mitkImageClone, itkImage);
    itkImage->FillBuffer(0);
    mitk::Image::Pointer mitkNewImage;
    mitk::CastToMitkImage(itkImage, mitkNewImage);

    if (NULL != mitkImage && NULL != itkImage)
    {
        mitk::PropertyList* pProList = mitkImage->GetPropertyList();
        std::string strPath;
        pProList->GetStringProperty("path", strPath);
        QString qSTRPath(strPath.c_str());

        int szNodules = m_NodulesInSlice.size();
        if (szNodules == 0)
        {
            MITK_WARN << "No ground truth loaded...";
            return;
        }

        for (int i = 0; i < szNodules; i++)
        {// towards all the nodules of this image set
            CNoduleInfoPerSlice*				pNodulesInSlice = m_NodulesInSlice[i];
            mitk::Image::ImageDataItemPointer	pThisSlice = mitkImage->GetSliceData(pNodulesInSlice->getSlicePosition());
            mitk::Image::ImageDataItemPointer	pNewSlice = mitkNewImage->GetSliceData(pNodulesInSlice->getSlicePosition());

            mitk::ImagePixelWriteAccessor<TPixelType, 2> writeAccess(mitkNewImage, pNewSlice);
            // 			mitk::ImagePixelReadAccessor<TPixelType, 2> readAccess(mitkNewImage, pNewSlice);

            itk::Index<2> write_idx;

            for (int j = 0; j < N_EXPERT_COUNT; j++)
            {// towards each expert on this slice

                int szThisSlice = 0;
                SNoduleInfo* pInfoArray = pNodulesInSlice->getNoduleInfoThisExpert(j, szThisSlice);

                if (NULL != pInfoArray)
                {
                    for (int k = 0; k < szThisSlice; k++)
                    {
                        SNoduleInfo&	sThisInfo = pInfoArray[k];
                        if (sThisInfo.eType == ENT_NonNodule)
                            continue;

                        SNoduleBounds&	sThisBounds = sThisInfo.sBounds;
                        int iBoundsWidth = (sThisBounds.iXMax - sThisBounds.iXMin);
                        int iBoundsHeight = (sThisBounds.iYMax - sThisBounds.iYMin);

                        if (iBoundsHeight <= 0 || iBoundsWidth <= 0){
                            continue;
                        }

                        int iCenterX = ceil((sThisBounds.iXMin + sThisBounds.iXMax) / 2.0);
                        int iCenterY = ceil((sThisBounds.iYMin + sThisBounds.iYMax) / 2.0);
                        // bounds may be out of N_NODULE_D!!
                        int iCenterBoundX1 = iCenterX - N_NODULE_R;
                        int iCenterBoundX2 = iCenterX + N_NODULE_R;
                        int iCenterBoundY1 = iCenterY - N_NODULE_R;
                        int iCenterBoundY2 = iCenterY + N_NODULE_R;

#ifdef _OUTPUT_2D_
                        // texture 
                        QString qSTRSuffix;
                        // folder_slicePos_expertidx_noduleidxOfThisSlice
                        qSTRSuffix.sprintf("%d_%d_%d_texture.txt", pNodulesInSlice->getSlicePosition(), j, k);
                        // open file writer pointers
                        QFile qTexture(strSaveFolder + qSTRSuffix);
                        QImage* qIMGTexture = new QImage(2 * N_NODULE_R + 1, 2 * N_NODULE_R + 1, QImage::Format_RGB32);
                        qIMGTexture->fill(QColor(0, 0, 0));

                        // ground truth
                        QString qSTRSuffix2;
                        // folder_slicePos_expertidx_noduleidxOfThisSlice
                        qSTRSuffix2.sprintf("%d_%d_%d_ground.txt", pNodulesInSlice->getSlicePosition(), j, k);
                        // open file writer pointers
                        QFile qGround(strSaveFolder + qSTRSuffix2);
                        QImage* qIMGGround = new QImage(2 * N_NODULE_R + 1, 2 * N_NODULE_R + 1, QImage::Format_RGB32);
                        qIMGGround->fill(QColor(0, 0, 0));

                        // open stream
                        if (!qTexture.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate) ||
                            !qGround.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
                        {
                            QMessageBox::warning(NULL, strSaveFolder, "can't open", QMessageBox::Yes);
                            return;
                        }
                        QTextStream tsTexture(&qTexture);
                        QTextStream tsGround(&qGround);

                        //-------------------------------------------------------------------------
                        // texture information
                        QRgb* uINTPtr = reinterpret_cast<QRgb*>(qIMGTexture->bits());
                        for (int iY = iCenterBoundY1; iY <= iCenterBoundY2; iY++)
                        {
                            for (int iX = iCenterBoundX1; iX <= iCenterBoundX2; iX++)
                            {
                                read_idx[0] = iX;
                                read_idx[1] = iY;
                                int value;
                                if (iX < 0 || iX > IMAGE_X_DIM)
                                    value = 0;
                                else if (iY < 0 || iY > IMAGE_Y_DIM)
                                    value = 0;
                                else
                                    value = readAccess.GetPixelByIndex(read_idx);

                                tsTexture << value << "\t";
                                // helper image
                                int iPixelValue;
                                if (value < -F_CT_WIDTH)
                                    iPixelValue = 0;
                                else if (value > F_CT_WIDTH)
                                    iPixelValue = 255;
                                else
                                    iPixelValue = round((value + F_CT_WIDTH) / (2 * F_CT_WIDTH) * 255);
                                *(uINTPtr++) = QColor(iPixelValue, iPixelValue, iPixelValue).rgb();
                            }
                            tsTexture << "\n";
                        }
#else

#endif

                        //-------------------------------------------------------------------------
                        // ground truth info
                        int iNewBoundX1 = std::min(iCenterBoundX1, sThisBounds.iXMin);
                        int iNewBoundX2 = std::max(iCenterBoundX2, sThisBounds.iXMax);
                        int iNewBoundY1 = std::min(iCenterBoundY1, sThisBounds.iYMin);
                        int iNewBoundY2 = std::max(iCenterBoundY2, sThisBounds.iYMax);

                        V2IPoints&	vContour = sThisInfo.vContour;
                        int			szContour = vContour.size();
                        int*		pGroundTruthMat = new int[N_NODULE_D * N_NODULE_D];
                        memset(pGroundTruthMat, 0, sizeof(int)*N_NODULE_D*N_NODULE_D);

                        // contour point look up
                        M2IPoints	mIndexedContour;
                        for_each(vContour.begin(), vContour.end(),
                            [&mIndexedContour](SPoint2I& spt)
                        {
                            mIndexedContour.insert(std::make_pair(spt.iIdx, spt));
                        });

                        //for (int iX = iCenterX - N_NODULE_R; iX <= iCenterX + N_NODULE_R; iX++)
                        for (int iX = iNewBoundX1; iX <= iNewBoundX2; iX++)
                        {
                            int iY1 = -1, iY2 = -1, iX1 = iX;

                            // step1, find top point
                            for (int iY = iNewBoundY1; iY <= iNewBoundY2; iY++)
                            {
                                int iIDX = iY * IMAGE_X_DIM + iX;
                                M2IPoints::iterator itFind = mIndexedContour.find(iIDX), itMContourEnd = mIndexedContour.end();
                                if (itFind != itMContourEnd)
                                {
                                    iY1 = iY;
                                    break;
                                }
                            }

                            // step2, find bottom point
                            for (int iY = iNewBoundY2; iY >= iNewBoundY1; iY--)
                            {
                                int iIDX = iY * IMAGE_X_DIM + iX;
                                M2IPoints::iterator itFind = mIndexedContour.find(iIDX), itMContourEnd = mIndexedContour.end();
                                if (itFind != itMContourEnd)
                                {
                                    iY2 = iY;
                                    break;
                                }
                            }

                            // step3, check and fill
                            if (iY1 >= 0 && iY2 >= 0)
                            {// inside the bound
                                // out of sample space
                                if (iX < iCenterBoundX1 || iX > iCenterBoundX2)
                                    continue;

                                // 								if (sThisBounds.iXMin == sThisBounds.iXMax && sThisBounds.iYMin == sThisBounds.iYMax)
                                // 								{
                                // 									int iY = sThisBounds.iYMin;
                                // 									iY -= iCenterY;; iX1 -= iCenterX;
                                // 									iY += N_NODULE_R; iX1 += N_NODULE_R;
                                // 
                                // 									int iIDX = iY * N_NODULE_D + iX1;
                                // 									// 									*(uINTPtr + iIDX) = QColor(255, 255, 255).rgb();
                                // 									pGroundTruthMat[iIDX] = 1;
                                // 								}
                                // 								else
                                // 								{
                                // 
                                iY1 = std::max(iY1, iCenterBoundY1);
                                iY2 = std::min(iY2, iCenterBoundY2);

                                // translate to local coordinate
                                iY1 -= iCenterY; iY2 -= iCenterY; iX1 -= iCenterX;
                                iY1 += N_NODULE_R; iY2 += N_NODULE_R; iX1 += N_NODULE_R;
                                iY1 += 1; iY2 -= 1;

                                if (iY2 >= iY1)
                                {
                                    for (int iY = iY1; iY <= iY2; iY++)
                                    {
                                        int iIDX = iY * N_NODULE_D + iX1;
                                        // 											*(uINTPtr + iIDX) = QColor(255, 255, 255).rgb();
                                        pGroundTruthMat[iIDX] = 1;
                                    }
                                }
                                // 								}
                            }
                            else
                            {// not inside, need further check
                                if (iX < iCenterBoundX1 || iX > iCenterBoundX2)
                                    continue;
                                else
                                {
                                    if (iY1 >= 0 && iY2 < 0)
                                    {
                                        iY1++;
                                        iY2 = iCenterBoundY2;
                                    }
                                    else if (iY1 < 0 && iY2 >= 0)
                                    {
                                        iY2--;
                                        iY1 = iCenterBoundY1;
                                    }
                                    else
                                    {// both are negative, not found
                                        continue;
                                        //iY1 = iCenterBoundY1;
                                        //iY2 = iCenterBoundY2;
                                    }


                                    // translate to local coordinate
                                    iY1 -= iCenterY; iY2 -= iCenterY; iX1 -= iCenterX;
                                    iY1 += N_NODULE_R; iY2 += N_NODULE_R; iX1 += N_NODULE_R;
                                    if (iY2 >= iY1)
                                    {
                                        // iY1 and iY2 out of current index, fill the whole line
                                        for (int iY = iY1; iY <= iY2; iY++)
                                        {
                                            int iIDX = iY * N_NODULE_D + iX1;
                                            //*(uINTPtr + iIDX) = QColor(255, 255, 255).rgb();
                                            pGroundTruthMat[iIDX] = 1;
                                        }
                                    }
                                }
                            }
                        }

#ifdef _OUTPUT_2D_
                        // write...
                        uINTPtr = reinterpret_cast<QRgb*>(qIMGGround->bits());
                        for (int iYMat = 0; iYMat < N_NODULE_D; iYMat++)
                        {
                            for (int iXMat = 0; iXMat < N_NODULE_D; iXMat++)
                            {
                                int iMatIdx = iYMat * N_NODULE_D + iXMat;
                                if (pGroundTruthMat[iMatIdx] > 0)
                                {
                                    *(uINTPtr + iMatIdx) = QColor(255, 255, 255).rgb();
                                    tsGround << 1.00000000 << "\t";
                                }
                                else
                                {
                                    *(uINTPtr + iMatIdx) = QColor(0, 0, 0).rgb();
                                    tsGround << 0.00000000 << "\t";
                                }
                            }
                            tsGround << "\n";
                        }
                        SAFE_DELETE(pGroundTruthMat);

                        // save and clean
                        qIMGTexture->save(strSaveFolder + qSTRSuffix + ".png");
                        qIMGGround->save(strSaveFolder + qSTRSuffix2 + ".png");
                        SAFE_DELETE(qIMGTexture);
                        SAFE_DELETE(qIMGGround);

                        // close file handler
                        qTexture.close();
                        qGround.close();

#else
                        for (int iYMat = 0; iYMat < N_NODULE_D; iYMat++)
                        {
                            for (int iXMat = 0; iXMat < N_NODULE_D; iXMat++)
                            {
                                int iMatIdx = iYMat * N_NODULE_D + iXMat;
                                int iRealX = iXMat - N_NODULE_R + iCenterX;
                                int iRealY = iYMat - N_NODULE_R + iCenterY;
                                write_idx[0] = iRealX;
                                write_idx[1] = iRealY;
                                int val = writeAccess.GetPixelByIndex(write_idx);

                                if (pGroundTruthMat[iMatIdx] > 0)
                                    writeAccess.SetPixelByIndex(write_idx, val + 1);

                                //                                 strDBG.sprintf("Writing %d %d %d", iRealX, iRealY, pNodulesInSlice->getSlicePosition());
                                //                                 MITK_INFO << strDBG.toStdString();
                            }
                        }
#endif

                    }
                }
            }
        }
#ifdef _OUTPUT_2D_
        SAFE_DELETE(pWorkDir);
#endif
    }

    // 	SAFE_DELETE(dwGaussian);

#ifdef _OUTPUT_2D_
#else
    mitk::Image::Pointer normGroundImage = normalizeImage(mitkNewImage);
    ImageType::Pointer itkNormGround;
    mitk::CastToItkImage(normGroundImage, itkNormGround);

    //  cast image to binary
    typedef itk::BinaryThresholdImageFilter<ImageType, ImageType>		ThresholdFilterType;
    ThresholdFilterType::Pointer thFilter = ThresholdFilterType::New();
    thFilter->SetInput(itkNormGround);
    thFilter->SetLowerThreshold(2);
    thFilter->SetUpperThreshold(255);
    thFilter->SetInsideValue(1);
    thFilter->SetOutsideValue(0);
    thFilter->Update();

    mitk::Image::Pointer mitkNormThreImage;
    mitk::CastToMitkImage(thFilter->GetOutput(), mitkNormThreImage);
    // allocate a new node and show
    mitk::DataNode::Pointer newNode = mitk::DataNode::New();
    newNode->SetData(mitkNormThreImage);
    newNode->SetProperty("name", mitk::StringProperty::New("normalized_groundtruth"));
    // add result to data tree
    this->GetDataStorage()->Add(newNode);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();

    mitk::Image::Pointer mitkNormTexImage = normalizeImage(mitkImage);
    // allocate a new node and show
    mitk::DataNode::Pointer newNode2 = mitk::DataNode::New();
    newNode2->SetData(mitkNormTexImage);
    newNode2->SetProperty("name", mitk::StringProperty::New("normalized_texture"));
    // add result to data tree
    this->GetDataStorage()->Add(newNode2);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();

#endif
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoExternalTransform()
{
    mitk::Image::Pointer mitkImage = getSelectedImage();
    ImageType::Pointer itkImage = NULL;
    if (NULL == mitkImage)
    {
        MITK_WARN << "Select Data Node First";
        return;
    }
    mitk::CastToItkImage(mitkImage, itkImage);

    m_vExternalRegions.clear();
    m_vExternalNoduleIndex.clear();

    typedef itk::MultiplyImageFilter<ImageType, ImageType> MultiplyCastType;
    MultiplyCastType::Pointer multiplyFilter = MultiplyCastType::New();
    multiplyFilter->SetInput(itkImage);
    multiplyFilter->SetConstant(255);
    multiplyFilter->Update();

    typedef itk::ConnectedComponentImageFilter <ImageType, ImageType> ConnectedComponentImageFilterType_3D;
    ConnectedComponentImageFilterType_3D::Pointer connectedComponentImageFilter_3D = ConnectedComponentImageFilterType_3D::New();
    connectedComponentImageFilter_3D->SetInput(multiplyFilter->GetOutput());
    connectedComponentImageFilter_3D->Update();

    // shape analysis
    typedef itk::LabelImageToShapeLabelMapFilter<ImageType> ShapeLabelMapFilter_3D;
    ShapeLabelMapFilter_3D::Pointer shape_filter_3d = ShapeLabelMapFilter_3D::New();
    shape_filter_3d->SetInput(connectedComponentImageFilter_3D->GetOutput());
    shape_filter_3d->Update();

    typedef itk::ShapeLabelObject<TPixelType, VDIM>	ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType >	ShapeLabelMapType;
    ShapeLabelMapType::Pointer shape_labelMap = shape_filter_3d->GetOutput();

    int szLabelMap = shape_labelMap->GetNumberOfLabelObjects();
    for (unsigned int n = 0; n < szLabelMap; ++n){
        ShapeLabelObjectType*	labelObject = shape_labelMap->GetNthLabelObject(n);
        ImageType::RegionType	bound = labelObject->GetBoundingBox();
        m_vExternalRegions.push_back(bound);
    }

    QString qSTRDBG;
    qSTRDBG.sprintf("Done! With %d Labels", szLabelMap);
    m_Controls.lbNoduleIndex->setText(qSTRDBG);
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoSelectExternalNodule()
{
    mitk::Image::Pointer	mitkFeatureImage = getSelectedImage();
    mitk::BaseGeometry::Pointer geometry = mitkFeatureImage->GetGeometry();

    QString strDBG;
    mitk::Point3D crossPositionInWorldCoordinates = m_IRenderWindowPart->GetSelectedPosition();
    strDBG.sprintf("Mouse Corr (%f %f %f)", crossPositionInWorldCoordinates[0], crossPositionInWorldCoordinates[1], crossPositionInWorldCoordinates[2]);
    MITK_INFO << strDBG.toStdString();

    int szStoredLabels = m_vExternalRegions.size();
    for (int i = 1; i < szStoredLabels; i++)
    {
        ImageType::RegionType& thisRegion = m_vExternalRegions[i];
        ImageType::IndexType index = thisRegion.GetIndex();
        ImageType::SizeType size = thisRegion.GetSize();

        mitk::Point3D cornerPoint1InWorldCoordinates;
        geometry->IndexToWorld(index, cornerPoint1InWorldCoordinates);

        if (crossPositionInWorldCoordinates[0] >= cornerPoint1InWorldCoordinates[0] && crossPositionInWorldCoordinates[0] <= cornerPoint1InWorldCoordinates[0] + size[0] &&
            crossPositionInWorldCoordinates[1] >= cornerPoint1InWorldCoordinates[1] && crossPositionInWorldCoordinates[1] <= cornerPoint1InWorldCoordinates[1] + size[1] &&
            crossPositionInWorldCoordinates[2] >= cornerPoint1InWorldCoordinates[2] && crossPositionInWorldCoordinates[2] <= cornerPoint1InWorldCoordinates[2] + size[2])
        {
            int szExternalNoduleIndex = m_vExternalNoduleIndex.size();
            bool bDuplicate = false;
            for (int j = 0; j < szExternalNoduleIndex; j++)
            {
                if (m_vExternalNoduleIndex[j] == i)
                {
                    bDuplicate = true;
                    break;
                }
            }
            if (!bDuplicate)
                m_vExternalNoduleIndex.push_back(i);

            QString qTMP;
            szExternalNoduleIndex = m_vExternalNoduleIndex.size();
            for (int j = 0; j < szExternalNoduleIndex; j++)
            {
                qTMP.append(QString::number(m_vExternalNoduleIndex[j], 10));
                qTMP.append(" ");
            }
            m_Controls.lbNoduleIndex->setText(qTMP);

            strDBG.sprintf("%d", i);
            QMessageBox msgBox;
            msgBox.setText(strDBG);
            msgBox.exec();
        }
    }
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::showNextImage()
{
    //     mitk::Image* pSelected = getSelectedImage();
    //     if (NULL != pSelected)
    {
        if (m_NodulesInSlice.empty())
            return;
        int szAllDrawableLabel = m_NodulesInSlice.size();
        if (iCurrengShownImage >= szAllDrawableLabel - 1)
            iCurrengShownImage = szAllDrawableLabel - 1;
        else
            iCurrengShownImage++;

        fillDataGrid(m_NodulesInSlice[iCurrengShownImage]);
        CMyDrawableNodule* pLabel = buildDrawableLabel(m_bgImage/*pSelected*/, m_NodulesInSlice[iCurrengShownImage]);
        m_Controls.scrollArea->setWidget(pLabel);
    }
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::showPriorImage()
{
    //     mitk::Image* pSelected = getSelectedImage();
    if (NULL != m_bgImage/*pSelected*/)
    {
        if (m_NodulesInSlice.empty())
            return;
        int szAllDrawableLabel = m_NodulesInSlice.size();
        if (iCurrengShownImage <= 0)
            iCurrengShownImage = 0;
        else
            iCurrengShownImage--;

        fillDataGrid(m_NodulesInSlice[iCurrengShownImage]);
        CMyDrawableNodule* pLabel = buildDrawableLabel(m_bgImage/*pSelected*/, m_NodulesInSlice[iCurrengShownImage]);
        m_Controls.scrollArea->setWidget(pLabel);
    }

}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::destroyContainers()
{
    MITK_INFO << "Destroying Internal Container";

    int szAllDrawableNodules = m_NodulesInSlice.size();
    for (int i = 0; i < szAllDrawableNodules; i++)
        SAFE_DELETE(m_NodulesInSlice[i]);
    m_NodulesInSlice.clear();

    m_NonNoduleInfoFromXML.clear();
    m_NonNoduleInfoFromLabel.clear();
    m_vNoduleRegions.clear();

    iCurrengShownImage = -1;

    SAFE_DELETE(m_outputNoduleTypes);

    destroyUIComponents();
}


//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::destroyUIComponents()
{
    // ui restore
    m_ReportModel->removeRows(0, m_ReportModel->rowCount());
    QLabel* pBlankLabel = new QLabel();
    m_Controls.scrollArea->setWidget(pBlankLabel);

    // remove widgets from grid
    int szItems = m_Controls.gridLayout->count();
    std::vector<QLayoutItem*> vItems;
    for (int i = 0; i < szItems; i++){
        vItems.push_back(m_Controls.gridLayout->itemAt(i));
    }
    for (int i = 0; i < szItems; i++)
    {
        m_Controls.gridLayout->removeItem(vItems[i]);
        vItems[i]->widget()->deleteLater();
    }
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoLoadAnnotations()
{
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) return;

    mitk::DataNode* node = nodes.front();

    if (!node)
    {
        // Nothing selected. Inform the user and return
        QMessageBox::information(NULL, "Error", "Please load and select an image before starting image processing.");
        return;
    }

    // a node itself is not very useful, we need its data item (the image)
    mitk::BaseData*		data = node->GetData();
    mitk::PropertyList* pPropertyList = node->GetData()->GetPropertyList();
    std::string			qstrPath = pPropertyList->GetProperty("path")->GetValueAsString();
    if (qstrPath == "")
        return;
    qstrPath = qstrPath.substr(0, qstrPath.find_last_of("/"));

    // test
    mitk::BaseGeometry::ConstPointer  geometry = data->GetGeometry();
    if (NULL != geometry)
    {
        mitk::BoundingBox::BoundsArrayType bounds = geometry->GetBounds();
        mitk::Point3D cornerPoint1InIndexCoordinates;
        cornerPoint1InIndexCoordinates[0] = bounds[0];
        cornerPoint1InIndexCoordinates[1] = bounds[2];
        cornerPoint1InIndexCoordinates[2] = bounds[4];

        mitk::Point3D cornerPoint2InIndexCoordinates;
        cornerPoint2InIndexCoordinates[0] = bounds[1];
        cornerPoint2InIndexCoordinates[1] = bounds[3];
        cornerPoint2InIndexCoordinates[2] = bounds[5];

        mitk::Point3D cornerPoint1InWorldCoordinates;
        mitk::Point3D cornerPoint2InWorldCoordinates;
        geometry->IndexToWorld(cornerPoint1InIndexCoordinates, cornerPoint1InWorldCoordinates);
        geometry->IndexToWorld(cornerPoint2InIndexCoordinates, cornerPoint2InWorldCoordinates);
        QString strDBG;
        strDBG.sprintf("corner 1:(%f, %f, %f), corner 2:(%f, %f, %f)",
            cornerPoint1InWorldCoordinates[0], cornerPoint1InWorldCoordinates[1], cornerPoint1InWorldCoordinates[2],
            cornerPoint2InWorldCoordinates[0], cornerPoint2InWorldCoordinates[1], cornerPoint2InWorldCoordinates[2]);
        MITK_INFO << strDBG;
    }


    // filter out the xml configuration
    QDir *dataPath = new QDir(qstrPath.c_str());
    QStringList filter;
    filter << "*.xml";
    dataPath->setNameFilters(filter);
    QList<QFileInfo> *fileInfo = new QList<QFileInfo>(dataPath->entryInfoList(filter));
    int nXMLFileCount = fileInfo->count();
    if (nXMLFileCount == 0)
    {
        MITK_WARN << "No XML Configuration In " << qstrPath;
        return;
    }
    else if (nXMLFileCount > 1)
    {
        MITK_WARN << "Multiple XML Configuration In " << qstrPath << " Using The First!";
    }
    QString strXMLConf = fileInfo->at(0).filePath();
    MITK_INFO << "Loading XML Configuration " << strXMLConf;

    if (data)
    {
        // test if this data item is an image or not (could also be a surface or something totally different)
        mitk::Image* image = dynamic_cast<mitk::Image*>(data);
        if (image)
        {
            // clone base data 
            m_bgImage = image->Clone();

            // zero, destroy
            destroyContainers();

            VVSNoduleInfo vTmpNodules;
            // first, load
            LoadAnnotationXML(strXMLConf, vTmpNodules);

            // second, analysis
            collectInfoFromXML(image, vTmpNodules);

            // last, label
            LabelImageUsingXML(image, vTmpNodules);

            // non nodule
            collectNonNoduleInfo(vTmpNodules, m_NonNoduleInfoFromXML);
        }
    }

    m_Controls.label->setText("Processing XML Done!!");
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoNormalization()
{
    mitk::Image*		mitkImage = getSelectedImage();
    mitk::DataNode*		mitkNode = getSelectedNode();

    QString qNodeName;
    qNodeName.sprintf("%s_resample", mitkNode->GetName().c_str());
    mitk::Image::Pointer pResult = normalizeImage(mitkImage);
    ShowProcessedDataNode(pResult, qNodeName.toStdString());
}

//////////////////////////////////////////////////////////////////////////
mitk::Image::Pointer LIDCXMLView::normalizeImage(mitk::Image::Pointer mitkImage)
{
    if (NULL != mitkImage)
    {
        typedef itk::ResampleImageFilter< ImageType, ImageType >			ResampleImageFilterType;
        typedef itk::LinearInterpolateImageFunction< ImageType, double >	LinearInterpolatorType;

        ImageType::Pointer itkImage;
        mitk::CastToItkImage(mitkImage->Clone(), itkImage);

        ResampleImageFilterType::Pointer	resampler = ResampleImageFilterType::New();
        LinearInterpolatorType::Pointer		interpolator = LinearInterpolatorType::New();
        resampler->SetInterpolator(interpolator);
        resampler->SetInput(itkImage);
        resampler->SetOutputOrigin(itkImage->GetOrigin());

        ImageType::SizeType			input_size = itkImage->GetLargestPossibleRegion().GetSize();
        ImageType::SpacingType		input_spacing = itkImage->GetSpacing();
        ImageType::SpacingValueType input_min_spacing = F_NORMALIZE_SPACE;

        ImageType::SizeType		output_size;
        ImageType::SpacingType	output_spacing;

        output_size[0] = input_size[0] * (input_spacing[0] / input_min_spacing);
        output_size[1] = input_size[1] * (input_spacing[1] / input_min_spacing);
        output_size[2] = input_size[2] * (input_spacing[2] / input_min_spacing);
        output_spacing[0] = input_min_spacing;
        output_spacing[1] = input_min_spacing;
        output_spacing[2] = input_min_spacing;

        QString debugSTR;
        debugSTR.sprintf("Resample...Original [Spacing(%f, %f, %f), Size(%f, %f, %f)], Output [Spacing(%f, %f, %f), Size(%f, %f, %f)]",
            input_spacing[0], input_spacing[1], input_spacing[2], input_size[0], input_size[1], input_size[2],
            output_spacing[0], output_spacing[1], output_spacing[2], output_size[0], output_size[1], output_size[2]);
        MITK_INFO << debugSTR;

        resampler->SetSize(output_size);
        resampler->SetOutputSpacing(output_spacing);
        resampler->SetOutputDirection(itkImage->GetDirection());
        resampler->Update();

        // show results
        mitk::Image::Pointer pResult;
        mitk::CastToMitkImage(resampler->GetOutput(), pResult);

        return pResult;
    }
    return NULL;
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::ShowProcessedDataNode(ImageType::Pointer itkImage, std::string strNodeName, bool bSetDefaultProperties /*= false*/)
{
    mitk::Image::Pointer resultImage = mitk::Image::New();
    mitk::CastToMitkImage(itkImage, resultImage);
    ShowProcessedDataNode(resultImage, strNodeName, bSetDefaultProperties);
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::ShowProcessedDataNode(mitk::Image::Pointer mitkImage, std::string strNodeName, bool bSetDefaultProperties /*= false*/)
{
    // allocate a new node and show
    mitk::DataNode::Pointer newNode = mitk::DataNode::New();
    newNode->SetData(mitkImage);
    newNode->SetProperty("name", mitk::StringProperty::New(strNodeName));

    if (bSetDefaultProperties)
    {
        // set some properties
        newNode->SetProperty("binary", mitk::BoolProperty::New(true));
        newNode->SetProperty("color", mitk::ColorProperty::New(1.0, 0.0, 0.0));
        newNode->SetProperty("volumerendering", mitk::BoolProperty::New(true));
        newNode->SetProperty("layer", mitk::IntProperty::New(1));
    }
    newNode->SetProperty("opacity", mitk::FloatProperty::New(0.5));

    // add result to data tree
    this->GetDataStorage()->Add(newNode);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();

}

//////////////////////////////////////////////////////////////////////////
bool LIDCXMLView::setLastOpenPath(QString q_path)
{
    berry::IPreferencesService* prefService = mitk::PluginActivator::GetInstance()->GetPreferencesService();

    berry::IPreferences::Pointer  prefs = berry::IPreferences::Pointer(0);
    if (prefService)
        prefs = prefService->GetSystemPreferences()->Node("/General");

    if (prefs.IsNotNull())
    {
        prefs->Put("LastFileOpenPath", q_path);
        prefs->Flush();
        return true;
    }

    return false;
}

//////////////////////////////////////////////////////////////////////////
QString LIDCXMLView::getLastOpenPath()
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

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::collectInfoFromXML(mitk::Image* mitkImage, VVSNoduleInfo& vvNoduleInfo)
{
    MITK_INFO << "Collecting information from xml";

    mitk::Image::ImageDataItemPointer	volume = mitkImage->GetVolumeData();
    mitk::ImageDescriptor*				imgDescriptor = mitkImage->GetImageDescriptor();
    mitk::BaseGeometry*					imggeometry = mitkImage->GetGeometry();

    // process
    int szExperts = /*m_NoduleInExpert*/vvNoduleInfo.size();
    int nTotalItem = 0;
    for (int i = 0; i < szExperts; i++)
    {
        VSNoduleInfo& vTmpNodule = /*m_NoduleInExpert*/vvNoduleInfo[i];
        int szVTmpNodule = vTmpNodule.size();
        for (int j = 0; j < szVTmpNodule; j++)
        {
            SNoduleInfo& sInfor = vTmpNodule[j];
            V2IPoints& vContour = sInfor.vContour;
            int nContour = vContour.size();
            int xBMin = 10000, xBMax = 0, yBMin = 10000, yBMax = 0;
            /// find binding box
            for (int k = 0; k < nContour; k++)
            {
                SPoint2I& spThisPt = vContour[k];
                if (spThisPt.x < xBMin)
                    xBMin = spThisPt.x;
                if (spThisPt.x > xBMax)
                    xBMax = spThisPt.x;
                if (spThisPt.y < yBMin)
                    yBMin = spThisPt.y;
                if (spThisPt.y > yBMax)
                    yBMax = spThisPt.y;
            }
            // update expert id
            sInfor.nExpertId = i;

            // update bounds
            sInfor.sBounds.iXMax = xBMax;
            sInfor.sBounds.iXMin = xBMin;
            sInfor.sBounds.iYMax = yBMax;
            sInfor.sBounds.iYMin = yBMin;

            // look up slice index
            mitk::Point3D pt3DWorld;
            itk::Index<3> pt3DUnit;
            pt3DWorld[0] = 1.0;
            pt3DWorld[1] = 1.0;
            pt3DWorld[2] = sInfor.dwZPosition;
            imggeometry->WorldToIndex(pt3DWorld, pt3DUnit);
            //MITK_INFO << "World Position: " << sInfor.dwZPosition;
            sInfor.nSlicePos = pt3DUnit[2];
        }
    }
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::collectNonNoduleInfo(VVSNoduleInfo& vvNoduleInfo, VSNoduleInfo& vOutNonNodules)
{
    QString qDBG;

    // each expert
    int nAllExperts = /*m_NoduleInExpert*/vvNoduleInfo.size();

    typedef std::map<double, SNoduleInfo> SliceNoduleMap;
    SliceNoduleMap slice_nodule_map;
    int nAllNonNodules = 0;
    for (int i = 0; i < nAllExperts; i++)
    {
        VSNoduleInfo& vThisExpert = /*m_NoduleInExpert*/vvNoduleInfo[i];
        int nThisExpert = vThisExpert.size();
        for (int j = 0; j < nThisExpert; j++)
        {
            SNoduleInfo& sThisExpThisInfo = vThisExpert[j];
            switch (sThisExpThisInfo.eType)
            {
            case ENT_Nodule:
                // not process
                break;

            case ENT_NonNodule:
                slice_nodule_map.insert(std::make_pair(sThisExpThisInfo.nSlicePos, sThisExpThisInfo));
                nAllNonNodules++;
                break;
            }
        }
    }//end each expert
    qDBG.sprintf("Non-nodules before/after (%d/%d) filtering.", nAllNonNodules, slice_nodule_map.size());
    MITK_INFO << qDBG;

    SliceNoduleMap::iterator it = slice_nodule_map.begin(), itEnd = slice_nodule_map.end();
    for (; it != itEnd; it++){
        vOutNonNodules.push_back(it->second);
    }
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::LabelImageUsingXML(mitk::Image* mitkImage, VVSNoduleInfo& vvNoduleInfo)
{
    MITK_INFO << "Labeling image using loaded xml";

    mitk::Image::ImageDataItemPointer	volume = mitkImage->GetVolumeData();
    mitk::ImageDescriptor*				imgDescriptor = mitkImage->GetImageDescriptor();
    mitk::BaseGeometry*					imggeometry = mitkImage->GetGeometry();

    int	nDims = imgDescriptor->GetNumberOfDimensions();
    if (nDims == 3)
    {
        const unsigned int* piDim = imgDescriptor->GetDimensions();
        mitk::PixelType pIxelType = mitkImage->GetPixelType();
        unsigned int iWidth = piDim[0];
        unsigned int iHeight = piDim[1];
        MITK_INFO << "(width, height, count)" << piDim[0] << ", " << piDim[1] << ", " << piDim[2];
        int nSlices = piDim[2];

        // nodule infor helper
        typedef std::map<int, CNoduleInfoPerSlice*> NoduleInfoPerSliceMap;
        NoduleInfoPerSliceMap mNodulesPerSlice;

        // each expert
        int nAllExperts = /*m_NoduleInExpert*/vvNoduleInfo.size();
        for (int i = 0; i < nAllExperts; i++)
        {
            VSNoduleInfo& vThisExpert = /*m_NoduleInExpert*/vvNoduleInfo[i];
            int nThisExpert = vThisExpert.size();
            for (int j = 0; j < nThisExpert; j++)
            {
                SNoduleInfo& sThisExpThisInfo = vThisExpert[j];
                switch (sThisExpThisInfo.eType)
                {
                case ENT_Nodule:
                {
                    V2IPoints& vContour = sThisExpThisInfo.vContour;
                    int nContour = vContour.size();
                    // find or add a new drawable nodule
                    NoduleInfoPerSliceMap::iterator itAllDrawableEnd = mNodulesPerSlice.end();
                    NoduleInfoPerSliceMap::iterator itFind = mNodulesPerSlice.find(sThisExpThisInfo.nSlicePos);
                    if (itFind != itAllDrawableEnd)
                    {// an existing drawable nodule
                        CNoduleInfoPerSlice* cThis = itFind->second;
                        cThis->addNoduleInfo(i, sThisExpThisInfo);
                    }
                    else
                    {
                        CNoduleInfoPerSlice* pNodule = new CNoduleInfoPerSlice();
                        pNodule->addNoduleInfo(i, sThisExpThisInfo);
                        mNodulesPerSlice.insert(std::make_pair(sThisExpThisInfo.nSlicePos, pNodule));
                    }
                }
                break;

                case ENT_NonNodule:
                    // no process
                    break;
                }
            }
        }//end each expert

        // transfer data to member variable
        for_each(mNodulesPerSlice.begin(), mNodulesPerSlice.end(),
            [this](pair<int, CNoduleInfoPerSlice*> it)
        {
            m_NodulesInSlice.push_back(it.second);
        });

        mNodulesPerSlice.clear();
    }
    else
    {
        MITK_WARN << "Dimension not match. 3D volume needed";
    }
}

//////////////////////////////////////////////////////////////////////////
CMyDrawableNodule* LIDCXMLView::buildDrawableLabel(mitk::Image* mitkImage, CNoduleInfoPerSlice* sThisExpThisInfo)
{
    // acquire base image data
    //	mitk::Image::ImageDataItemPointer pThisSlice = mitkImage->GetSliceData(sThisExpThisInfo->getSlicePosition());
    // 	mitk::ImageVtkWriteAccessor *pVTKWrite = pThisSlice->GetVtkImageAccessor(mitkImage);
    // 	pNodule->setBaseImageFromVtkData(pVTKWrite->GetVtkImageData());
    std::string pType = mitkImage->GetPixelType().GetComponentTypeAsString();
    CMyDrawableNodule* pNodule = new CMyDrawableNodule();
    if (pType == "short")
    {
        ImageType::Pointer		itknewimage = ImageType::New();
        mitk::CastToItkImage(mitkImage, itknewimage);
        mitk::Image::Pointer	mitknewimage = mitk::Image::New();
        mitk::CastToMitkImage(itknewimage, mitknewimage);
        pNodule->setBaseImageData(mitknewimage, sThisExpThisInfo->getSlicePosition());
    }
    else
    {
        pNodule->setBaseImageData(mitkImage, sThisExpThisInfo->getSlicePosition());
    }

    // build drawable label
    for (int i = 0; i < N_EXPERT_COUNT; i++)
    {
        int szThisNodule = 0;
        SNoduleInfo* pNoduleInfo = sThisExpThisInfo->getNoduleInfoThisExpert(i, szThisNodule);
        if (NULL == pNoduleInfo)
            continue;

        for (int j = 0; j < szThisNodule; j++)
        {
            pNodule->addNoduleInfo(i, pNoduleInfo[j]);
        }
    }
    return pNodule;
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoOutputNonExternalNodule()
{
    QString strQDebug;

    QString qSTRSuffix = "annTrainingData_ExtNN.txt";
    QString tmpPath = m_WorkingDir + qSTRSuffix;

    QFile qTraining(tmpPath);
    QTextStream tsTraining(&qTraining);
    MITK_INFO << "Writing to " << tmpPath;

    mitk::DataStorage::Pointer pDataStorage = GetDataStorage();
    mitk::DataNode::Pointer pSelectedNode = pDataStorage->GetNamedNode("image_extract_bg_-2048");
    mitk::BaseData*	data = pSelectedNode->GetData();
    mitk::Image::Pointer mitkTexture = NULL;
    if (data){
        mitkTexture = dynamic_cast<mitk::Image*>(data);
    }
    else{
        MITK_WARN << "Not Found image_extract_bg_-2048";
        return;
    }

    int iWidth = mitkTexture->GetDimension(0);
    int iHeight = mitkTexture->GetDimension(1);
    int iSlice = mitkTexture->GetDimension(2);

    ImageType::Pointer itkTexture;
    mitk::CastToItkImage(mitkTexture, itkTexture);

    // file header: sz_data szInput szOutput
    int szExternalRegions = m_vExternalRegions.size();
    int szNoduleExternal = m_vExternalNoduleIndex.size();
    int szNonNoduleExternal = szExternalRegions - szNoduleExternal;
    int szData = N_NON_NODULE_D * N_NON_NODULE_D * N_NON_NODULE_D;
    int szInput = N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D;
    int szOutput = 1;

    strQDebug.sprintf("Ext total %d, Nodule %d, Non-Nodule %d", szExternalRegions, szNoduleExternal, szNonNoduleExternal);
    MITK_INFO << strQDebug.toStdString();

    int*nSignExternalRegions = new int[szExternalRegions];
    memset(nSignExternalRegions, 0, szExternalRegions*sizeof(int));
    for (int i = 0; i < szNoduleExternal; i++){
        nSignExternalRegions[m_vExternalNoduleIndex[i]] = -1;
    }
    //     for (int i = 0; i < szExternalRegions; i++)
    //         MITK_INFO << nSignExternalRegions[i];

    if (qTraining.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
    {
        tsTraining.setRealNumberNotation(QTextStream::FixedNotation);
        tsTraining.setRealNumberPrecision(PIXEL_PRECISION);
        if (szExternalRegions == 0)
        {
            tsTraining << 0 << "\t" << szInput << "\t" << szOutput << "\n";
            qTraining.close();
            return;
        }

        tsTraining << szNonNoduleExternal * szData << "\t" << szInput << "\t" << szOutput << "\n";

        for (int iNonNodule = 0; iNonNodule < szExternalRegions; iNonNodule++)
        {
            if (nSignExternalRegions[iNonNodule] < 0)
                continue;

            ImageType::RegionType non_nodule_region = m_vExternalRegions[iNonNodule];
            ImageType::IndexType non_nodule_region_index = non_nodule_region.GetIndex();
            ImageType::SizeType non_nodule_region_size = non_nodule_region.GetSize();

            ImageType::IndexType center_index;
            center_index[0] = non_nodule_region_index[0] + non_nodule_region_size[0] / 2;
            center_index[1] = non_nodule_region_index[1] + non_nodule_region_size[1] / 2;
            center_index[2] = non_nodule_region_index[2] + non_nodule_region_size[2] / 2;

            strQDebug.sprintf("processing non-nodule %d, center(%d,%d,%d), ", iNonNodule, center_index[0], center_index[1], center_index[2]);
            MITK_INFO << strQDebug.toStdString();

            //DoubleImageType::IndexType	access_index;
            ImageType::SizeType			regionSize;
            ImageType::IndexType		regionIndex;
            ImageType::RegionType		region;
            ImageType::SizeType			radius;

            regionSize[0] = N_NON_NODULE_D;
            regionSize[1] = N_NON_NODULE_D;
            regionSize[2] = N_NON_NODULE_D;
            regionIndex[0] = center_index[0] - N_NON_NODULE_R;
            regionIndex[1] = center_index[1] - N_NON_NODULE_R;
            regionIndex[2] = center_index[2] - N_NON_NODULE_R;

            if (regionIndex[0] < 0){
                regionIndex[0] = 0;
                MITK_WARN << "Iterating Region Index0 < 0";
            }
            if (regionIndex[0] > iWidth){
                regionIndex[0] = iWidth;
                MITK_WARN << "Iterating Region Index0 > Width";
            }
            if (regionIndex[1] < 0){
                regionIndex[1] = 0;
                MITK_WARN << "Iterating Region Index1 < 0";
            }
            if (regionIndex[1] > iHeight){
                regionIndex[1] = iHeight;
                MITK_WARN << "Iterating Region Index1 > Height";
            }
            if (regionIndex[2] < 0){
                regionIndex[2] = 0;
                MITK_WARN << "Iterating Region Index2 < 0";
            }
            if (regionIndex[2] > iSlice){
                regionIndex[2] = iSlice;
                MITK_WARN << "Iterating Region Index2 > Slice";
            }

            // region size
            region.SetSize(regionSize);
            // region start
            region.SetIndex(regionIndex);
            radius[0] = N_LOCAL_WINDOW_R;
            radius[1] = N_LOCAL_WINDOW_R;
            radius[2] = N_LOCAL_WINDOW_R;

            itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
            int nCounter = 0;
            while (!iterator.IsAtEnd())
            {
                ImageType::IndexType local_centerIndex = iterator.GetIndex();
                int nTotalIn = 0;
                for (unsigned int iITER = 0; iITER < iterator.Size(); ++iITER)
                {
                    ImageType::IndexType neighbour_index = iterator.GetIndex(iITER);

                    bool IsInBounds;
                    int neighborValue = iterator.GetPixel(iITER, IsInBounds);
                    if (!IsInBounds){
                        neighborValue = -F_CT_WIDTH;
                    }
                    if (neighborValue > F_CT_WIDTH){
                        neighborValue = F_CT_WIDTH;
                    }
                    if (neighborValue < -F_CT_WIDTH){
                        neighborValue = -F_CT_WIDTH;
                    }

                    double dwNormalized = ((double)(neighborValue + F_CT_WIDTH)) / ((double)(F_CT_WIDTH * 2));
                    // strQDebug.sprintf("[%d %d %d, %d]-%f", i, j, k, IsInBounds ? 1 : 0, dwNormalized);
                    // tsTraining << strQDebug << "\t";
                    tsTraining << /*neighborValue*/dwNormalized << "\t";
                    nTotalIn++;
                }
                tsTraining << "\n";

                // labeled as non-nodules
                double dwOutput = (double)MyMath::Random(0, 10000) / pow(10.0, RESULT_PRECISION) + 0.000000000001;
                QString qSTRDWResult;
                qSTRDWResult = QString::number(dwOutput, 'f', RESULT_PRECISION);

                tsTraining << qSTRDWResult;
                tsTraining << "\n";

                if (nTotalIn != szInput)
                {
                    strQDebug.sprintf("number of output training data not equal with sample size, Needed: %d, Actual: %d", szInput, nTotalIn);
                    MITK_WARN << strQDebug;
                }

                // do not forget move to next data
                ++iterator;
                nCounter++;
            }
            if (nCounter != szData)
            {
                strQDebug.sprintf("Nodule: %d, Total Need: %d, Actual: %d", iNonNodule, szData, nCounter);
                MITK_WARN << strQDebug;
            }
        }
        qTraining.close();
    }
    else
    {
        MITK_INFO << "Open" << tmpPath << " Failed";
        return;
    }

    SAFE_DELETE(nSignExternalRegions);
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::LoadAnnotationXML(QString filename, VVSNoduleInfo& vvNoduleInfo)
{
    MITK_INFO << "Loading XML Conf: " << filename;

    if (filename.isEmpty())
    {
        return;
    }
    else
    {
        QFile file(filename);
        QXmlStreamReader streamXML;
        if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::critical(NULL, tr("Error"),
                tr("Cannot read file %1").arg(filename));
            return;
        }
        streamXML.setDevice(&file);

        while (!streamXML.atEnd())
        {
            streamXML.readNext();
            if (streamXML.isStartElement() && streamXML.name() == "readingSession")
            {
                VSNoduleInfo vInfoThisRad;
                SNoduleInfo sNodInfo;
                SPoint2I	spNodePoint;

                while (!streamXML.atEnd())
                {
                    streamXML.readNext();
                    // nodule, extract roi
                    if (streamXML.isStartElement() && streamXML.name() == "unblindedReadNodule")
                    {
                        sNodInfo.clear();
                        while (!streamXML.atEnd())
                        {
                            streamXML.readNext();
                            if (streamXML.isStartElement() && streamXML.name() == "roi")
                            {
                                sNodInfo.eType = ENT_Nodule;

                                QString strzPos, strUID, strInclusion, strXCoord, strYCoord;
                                streamXML.readNext();

                                while (!streamXML.atEnd() && !streamXML.isStartElement() && streamXML.name() != "imageZposition")
                                    streamXML.readNext();
                                strzPos = streamXML.readElementText();
                                sNodInfo.dwZPosition = strzPos.toDouble();

                                while (!streamXML.atEnd() && !streamXML.isStartElement() && streamXML.name() != "imageSOP_UID")
                                    streamXML.readNext();
                                strUID = streamXML.readElementText();
                                sNodInfo.strUID = strUID.toStdString();

                                while (!streamXML.atEnd() && !streamXML.isStartElement() && streamXML.name() != "inclusion")
                                    streamXML.readNext();
                                strInclusion = streamXML.readElementText();
                                if (strInclusion == "TRUE")
                                    sNodInfo.bInclusion = true;
                                else
                                    sNodInfo.bInclusion = false;

                                while (!streamXML.atEnd())
                                {
                                    streamXML.readNext();
                                    if (streamXML.isStartElement() && streamXML.name() == "edgeMap")
                                    {
                                        while (!streamXML.atEnd() && streamXML.name() != "xCoord")
                                            streamXML.readNext();
                                        strXCoord = streamXML.readElementText();
                                        while (!streamXML.atEnd() && streamXML.name() != "yCoord")
                                            streamXML.readNext();
                                        strYCoord = streamXML.readElementText();
                                        spNodePoint.x = strXCoord.toInt();
                                        spNodePoint.y = strYCoord.toInt();
                                        spNodePoint.iIdx = spNodePoint.y*IMAGE_X_DIM + spNodePoint.x;
                                        sNodInfo.vContour.push_back(spNodePoint);
                                    }
                                    else if (streamXML.isEndElement() && streamXML.name() == "roi")
                                    {
                                        break;
                                    }
                                }
                                vInfoThisRad.push_back(sNodInfo);
                            }
                            else if (streamXML.isEndElement() && streamXML.name() == "unblindedReadNodule")
                            {
                                break;
                            }
                        }
                    }
                    // non-nodules
                    else if (streamXML.isStartElement() && streamXML.name() == "nonNodule")
                    {

                        sNodInfo.clear();
                        QString strzPos, strUID, strInclusion, strXCoord, strYCoord;
                        streamXML.readNext();

                        while (!streamXML.atEnd() && streamXML.name() != "imageZposition")
                            streamXML.readNext();
                        strzPos = streamXML.readElementText();
                        sNodInfo.dwZPosition = strzPos.toDouble();

                        while (!streamXML.atEnd() && streamXML.name() != "imageSOP_UID")
                            streamXML.readNext();
                        strUID = streamXML.readElementText();
                        sNodInfo.strUID = strUID.toStdString();

                        while (!streamXML.atEnd() && streamXML.name() != "xCoord")
                            streamXML.readNext();
                        strXCoord = streamXML.readElementText();
                        while (!streamXML.atEnd() && streamXML.name() != "yCoord")
                            streamXML.readNext();
                        strYCoord = streamXML.readElementText();
                        spNodePoint.x = strXCoord.toInt();
                        spNodePoint.y = strYCoord.toInt();
                        spNodePoint.iIdx = spNodePoint.y*IMAGE_X_DIM + spNodePoint.x;
                        sNodInfo.vContour.push_back(spNodePoint);

                        sNodInfo.eType = ENT_NonNodule;
                        vInfoThisRad.push_back(sNodInfo);
                    }
                    else if (streamXML.isEndElement() && streamXML.name() == "readingSession")
                    {
                        /*m_NoduleInExpert*/vvNoduleInfo.push_back(vInfoThisRad);
                        break;
                    }
                }
            }
        }
        file.close();
    }
}

////////////////////////////////////////////////////////////////////////////
//void LIDCXMLView::saveSpecificImage()
//{
//	QString strImageIndex = m_Controls.txtSpecificImage->text();
//	if (strImageIndex.isEmpty())
//		QMessageBox::warning(NULL, "Warning", "Specify the Image to be Saved !", QMessageBox::Yes);
//
//	QFileDialog* openFilePath = new QFileDialog(NULL, "Select Folder", "file");
//	openFilePath->setFileMode(QFileDialog::DirectoryOnly);
//	QString fileName = openFilePath->getExistingDirectory();
//
//	if (!fileName.isNull())
//	{
//		mitk::Image* pSelectedImage = getSelectedImage();
//		if (NULL != pSelectedImage)
//		{
//			mitk::BaseGeometry* pGeometry = pSelectedImage->GetGeometry();
//			mitk::Image::ImageDataItemPointer pThisSlice = pSelectedImage->GetSliceData(strImageIndex.toInt());
//			int iXSIZE = pSelectedImage->GetDimension(0);
//			int iYSIZE = pSelectedImage->GetDimension(1);
//
//			// 
//			fileName = fileName + "/" + strImageIndex;
//
//			if (NULL != pThisSlice)
//			{
//				QFile qTextFile(fileName + ".txt");
//				// open stream
//				if (!qTextFile.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
//				{
//					QMessageBox::warning(NULL, fileName, "can't open", QMessageBox::Yes);
//					return;
//				}
//
//				QTextStream qTxtStream(&qTextFile);
//				QImage* qIMGTexture = new QImage(iXSIZE, iYSIZE, QImage::Format_RGB32);
//				QRgb* uINTPtr = reinterpret_cast<QRgb*>(qIMGTexture->bits());
//
//				mitk::ImagePixelReadAccessor<int, 2> readAccess(pSelectedImage, pThisSlice);
//				itk::Index<2> idx;
//				for (int iY = 0; iY < iYSIZE; iY++)
//				{
//					for (int iX = 0; iX < iXSIZE; iX++)
//					{
//						idx[0] = iX;
//						idx[1] = iY;
//						int value = readAccess.GetPixelByIndex(idx);
//
//						// helper image
//						int iPixelValue;
//						if (value < -F_CT_WIDTH)
//							iPixelValue = 0;
//						else if (value > F_CT_WIDTH)
//							iPixelValue = 255;
//						else
//							iPixelValue = round((value + F_CT_WIDTH) / (2 * F_CT_WIDTH) * 255);
//
//						qTxtStream << iPixelValue << "\t";
//						*(uINTPtr++) = QColor(iPixelValue, iPixelValue, iPixelValue).rgb();
//					}
//					qTxtStream << "\n";
//				}
//				qTextFile.flush();
//				qTextFile.close();
//				qIMGTexture->save(fileName + ".png");
//			}
//		}
//	}
//}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoOutputTraining()
{
    mitk::DataNode::Pointer textureNode = m_Controls.comboTexture->GetSelectedNode();
    mitk::DataNode::Pointer maskNode = m_Controls.comboGround->GetSelectedNode();
    mitk::DataNode::Pointer nonnoduleNode = m_Controls.comboNonnodule->GetSelectedNode();

    if (NULL == textureNode || NULL == maskNode || NULL == nonnoduleNode){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    mitk::Image::Pointer	mitkTexture = dynamic_cast<mitk::Image*> (textureNode->GetData());
    mitk::Image::Pointer	mitkMask = dynamic_cast<mitk::Image*>(maskNode->GetData());
    mitk::Image::Pointer	mitkNonNodule = dynamic_cast<mitk::Image*>(nonnoduleNode->GetData());

    // check if images are valid
    if ((!mitkTexture) || (!mitkMask) || (!mitkNonNodule) || (mitkTexture->IsInitialized() == false) || (mitkMask->IsInitialized() == false) || (mitkNonNodule->IsInitialized() == false))
    {
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    //     QString working_dir = getLastOpenPath();
    if (m_WorkingDir == "" || m_WorkingDir == "\\" || m_WorkingDir == "/")
    {
        MITK_WARN << "Select Working Dir First";
        return;
    }

    outputNoduleTrainingData(m_WorkingDir, mitkTexture, mitkMask);
    outputNonNoduleTrainingData(m_WorkingDir, mitkTexture, mitkNonNodule);
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoSelectWorkingDir()
{
    //     QString tmp = QFileDialog::getExistingDirectory(
    //         NULL, "Select Working Directory", "E:\\LIDC\\LIDC-IDRI", QFileDialog::ShowDirsOnly);
    // 
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
void LIDCXMLView::DoResampleNonNodule()
{
    QString strDBG;

    mitk::Image* mitkImage = getSelectedImage();
    mitk::Image::Pointer mitkImageClone = mitkImage->Clone();

    ImageType::Pointer itkImageClone;
    mitk::CastToItkImage(mitkImageClone, itkImageClone);
    itkImageClone->FillBuffer(0);
    mitk::Image::Pointer mitkTobeResmapleImage;
    mitk::CastToMitkImage(itkImageClone, mitkTobeResmapleImage);

    mitk::Image::Pointer mitkResampledImage = normalizeImage(mitkTobeResmapleImage);
    ImageType::Pointer itkImage;
    mitk::CastToItkImage(mitkResampledImage, itkImage);

    if (NULL != mitkImage && NULL != itkImage)
    {
        ImageType::SizeType			input_size = itkImageClone->GetLargestPossibleRegion().GetSize();
        ImageType::SpacingType		input_spacing = itkImageClone->GetSpacing();
        ImageType::SpacingValueType input_min_spacing = F_NORMALIZE_SPACE;
        double dwScaleRate[3] = { 0.0 };

        dwScaleRate[0] = (input_spacing[0] / input_min_spacing);
        dwScaleRate[1] = (input_spacing[1] / input_min_spacing);
        dwScaleRate[2] = (input_spacing[2] / input_min_spacing);

        int szNonNodules = m_NonNoduleInfoFromXML.size();
        if (szNonNodules == 0)
        {
            MITK_WARN << "No non-nodules loaded......";
            return;
        }

        for (int i = 0; i < szNonNodules; i++)
        {
            SNoduleInfo& non_noduleinfo = m_NonNoduleInfoFromXML[i];
            V2IPoints& vcontour = non_noduleinfo.vContour;
            if (vcontour.size() == 0){
                strDBG.sprintf("Non-nodule %d is with empty contour", i);
                MITK_WARN << strDBG.toStdString();
                continue;
            }

            SPoint2I& spt = vcontour[0];

            ImageType::IndexType regionIndex;
            regionIndex[0] = round(spt.x * dwScaleRate[0]);
            regionIndex[1] = round(spt.y * dwScaleRate[1]);
            regionIndex[2] = round(non_noduleinfo.nSlicePos * dwScaleRate[2]);

            ImageType::SizeType regionSize;
            regionSize[0] = 1;
            regionSize[1] = 1;
            regionSize[2] = 1;

            strDBG.sprintf("working on nonnodule %d,(%d %d %d)->(%d %d %d)", i,
                spt.x, spt.y, non_noduleinfo.nSlicePos,
                regionIndex[0], regionIndex[1], regionIndex[2]);

            MITK_INFO << strDBG.toStdString();

            // start at index and look for neighbors at each point inside of index+size;
            ImageType::RegionType region;
            // region size
            region.SetSize(regionSize);
            // region start
            region.SetIndex(regionIndex);

            ImageType::SizeType radius;
            radius[0] = N_NON_NODULE_R;
            radius[1] = N_NON_NODULE_R;
            radius[2] = N_NON_NODULE_R;

            itk::NeighborhoodIterator<ImageType> iterator(radius, itkImage, region);
            while (!iterator.IsAtEnd())
            {
                ImageType::IndexType centerIndex = iterator.GetIndex();
                int nCount = 0;
                for (unsigned int iITER = 0; iITER < iterator.Size(); ++iITER)
                {
                    ImageType::IndexType neighbour_index = iterator.GetIndex(iITER);
                    bool IsInBounds;
                    int neighborValue = iterator.GetPixel(iITER, IsInBounds);
                    if (IsInBounds){
                        //itkImage->SetPixel(neighbour_index, 1);
                        iterator.SetPixel(iITER, 1);
                    }
                    nCount++;
                }
                //writeAccess.SetPixelByIndexSafe(centerIndex, 255);
                if (nCount != N_NON_NODULE_D*N_NON_NODULE_D*N_NON_NODULE_D)
                {
                    strDBG.sprintf("Wrong iterate count, want %d, actual %d", N_NODULE_D*N_NODULE_D*N_NODULE_D, nCount);
                    MITK_WARN << strDBG.toStdString();
                }

                ++iterator;
            }
        }

        // 	SAFE_DELETE(dwGaussian);

        mitk::Image::Pointer mitkNewImage;
        mitk::CastToMitkImage(itkImage, mitkNewImage);
        //         // allocate a new node and show
        mitk::DataNode::Pointer newNode = mitk::DataNode::New();
        newNode->SetData(mitkNewImage/*normalizedNonnodule*/);
        newNode->SetProperty("name", mitk::StringProperty::New("normalized_nonnodule"));
        // add result to data tree
        this->GetDataStorage()->Add(newNode);
        mitk::RenderingManager::GetInstance()->RequestUpdateAll();

        // 	QString qNewMitkImageName = strSaveFolder + "noduleGroundTruth.nrrd";
        // 	mitk::IOUtil::Save(mitkNewImage, qNewMitkImageName.toStdString());

    }
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoResampleAll()
{
    DoResampleGroundAndTexture();
    DoResampleNonNodule();
}


//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::DoOutputBOWTraining()
{
    mitk::DataNode::Pointer textureNode = m_Controls.comboTexture->GetSelectedNode();
    mitk::DataNode::Pointer maskNode = m_Controls.comboGround->GetSelectedNode();
    mitk::DataNode::Pointer nonnoduleNode = m_Controls.comboNonnodule->GetSelectedNode();

    if (NULL == textureNode || NULL == maskNode || NULL == nonnoduleNode){
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    mitk::Image::Pointer	mitkTexture = dynamic_cast<mitk::Image*> (textureNode->GetData());
    mitk::Image::Pointer	mitkMask = dynamic_cast<mitk::Image*>(maskNode->GetData());
    mitk::Image::Pointer	mitkNonNodule = dynamic_cast<mitk::Image*>(nonnoduleNode->GetData());

    // check if images are valid
    if ((!mitkTexture) || (!mitkMask) || (!mitkNonNodule) || (mitkTexture->IsInitialized() == false) || (mitkMask->IsInitialized() == false) || (mitkNonNodule->IsInitialized() == false))
    {
        MITK_WARN << "At least one of the input images are broken or not initialized. Returning";
        return;
    }

    //     QString working_dir = getLastOpenPath();
    if (m_WorkingDir == "" || m_WorkingDir == "\\" || m_WorkingDir == "/")
    {
        MITK_WARN << "Select Working Dir First";
        return;
    }

    QString q_tmp_work_dir = m_WorkingDir;
    QDir q_dir;
    q_tmp_work_dir += "bow";
    if (!q_dir.exists(q_tmp_work_dir)){
        q_dir.mkdir(q_tmp_work_dir);
    }

    // arrange dirs
    q_tmp_work_dir += "\\training";
    if (!q_dir.exists(q_tmp_work_dir)){
        q_dir.mkdir(q_tmp_work_dir);
    }
    else{
        clearDirectory(q_tmp_work_dir);
    }

    // dir has been cleared in last step
    QString q_this_path;
    for (int i = 0; i < EOutType_COUNT; i++)
    {
        q_this_path.sprintf("%s\\%s", q_tmp_work_dir.toStdString().c_str(), STR_OUTPUT_TYPE[i].toStdString().c_str());
        q_dir.mkdir(q_this_path);
    }

    // make non-nodule dirs
    q_this_path.sprintf("%s\\%s", q_tmp_work_dir.toStdString().c_str(), "NONNODULE");
    q_dir.mkdir(q_this_path);

    // last element used to store non-nodule infor.
    int szNodulesThisType[EOutType_COUNT + 1] = { 0 };
    outputNoduleBOWTrainingData(q_tmp_work_dir, mitkTexture, mitkMask, szNodulesThisType);
    outputNonNoduleBOWTrainingData(q_tmp_work_dir, mitkTexture, mitkNonNodule, szNodulesThisType[EOutType_COUNT]);

    // write config for better extraction
    QFile qFileConfig(q_tmp_work_dir + "\\config.txt");
    QTextStream tsConfig(&qFileConfig);
    if (qFileConfig.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate)){
        for (int i = 0; i < EOutType_COUNT + 1; i++)
            tsConfig << szNodulesThisType[i] << "\n";
        qFileConfig.close();
    }
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::outputNoduleBOWTrainingData(QString strSavePath, mitk::Image::Pointer mitkTexture, mitk::Image::Pointer mitkGroundTruth, int* szNodulesThisType)
{
    // for debugging formation
    QString strQDebug;

    mitk::Image::Pointer mitkTMP = mitkGroundTruth->Clone();
    DoubleImageType::Pointer itkDoubleImage;
    mitk::CastToItkImage(mitkTMP, itkDoubleImage);
    itkDoubleImage->FillBuffer(0.0);
    mitk::Image::Pointer mitkDoubleImage;
    mitk::CastToMitkImage(itkDoubleImage, mitkDoubleImage);

    int iWidth = mitkTexture->GetDimension(0);
    int iHeight = mitkTexture->GetDimension(1);
    int iSlice = mitkTexture->GetDimension(2);

    if (NULL == m_outputNoduleTypes){
        MITK_WARN << "Output Types Empty...";
        return;
    }

    // collect count of each type nodule
    //     int szNodulesThisType[N_NODULE_TYPES] = { 0 };
    int szAllNodules = m_vNoduleRegions.size();

    // transform the mitk texture
    ImageType::Pointer itkTexture;
    mitk::CastToItkImage(mitkTexture, itkTexture);

    for (int i = 0; i < EOutType_COUNT; i++)
        szNodulesThisType[i] = 0;

    for (int iNodule = 0; iNodule < szAllNodules; iNodule++)
    {
        ImageType::RegionType& nodule_region = m_vNoduleRegions[iNodule];
        ImageType::IndexType nodule_index = nodule_region.GetIndex();
        ImageType::SizeType nodule_size = nodule_region.GetSize();
        int eOutputNoduleType = m_outputNoduleTypes[iNodule];
        szNodulesThisType[eOutputNoduleType]++;

        QString q_bow_output_filename;
        q_bow_output_filename.sprintf("%s\\%s\\%d.txt", strSavePath.toStdString().c_str(),
            STR_OUTPUT_TYPE[eOutputNoduleType].toStdString().c_str(),
            szNodulesThisType[eOutputNoduleType]);

        MITK_INFO << "Writing to " << q_bow_output_filename.toStdString();

        QFile qBOWOutput(q_bow_output_filename);
        QTextStream tsBOWOutput(&qBOWOutput);
        if (qBOWOutput.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
        {
            tsBOWOutput.setRealNumberNotation(QTextStream::FixedNotation);
            // 1/2000 = 0.0005, so 4 digits after . is enough
            tsBOWOutput.setRealNumberPrecision(PIXEL_PRECISION);
        }
        else
        {
            MITK_ERROR << "Open " << q_bow_output_filename.toStdString() << " Failed";
            return;
        }

        ImageType::IndexType center_index;
        center_index[0] = nodule_index[0] + nodule_size[0] / 2;
        center_index[1] = nodule_index[1] + nodule_size[1] / 2;
        center_index[2] = nodule_index[2] + nodule_size[2] / 2;

        // do through mitk image iteration
        mitk::ImagePixelReadAccessor<TPixelType, 3>			ground_access(mitkGroundTruth);
        mitk::ImagePixelWriteAccessor<DoublePixelType, 3>	double_image_access(mitkDoubleImage);

        //DoubleImageType::IndexType	access_index;
        ImageType::SizeType			regionSize;
        ImageType::IndexType		regionIndex;
        ImageType::RegionType		region;
        ImageType::SizeType			radius;

        regionSize[0] = /*N_NODULE_D*/N_NODULE_BOW_D;
        regionSize[1] = N_NODULE_BOW_D;
        regionSize[2] = N_NODULE_BOW_D;
        regionIndex[0] = center_index[0] - N_NODULE_BOW_R - 1;
        regionIndex[1] = center_index[1] - N_NODULE_BOW_R - 1;
        regionIndex[2] = center_index[2] - N_NODULE_BOW_R - 1;

        if (regionIndex[0] < 0){
            regionIndex[0] = 0;
            MITK_WARN << "Iterating Region Index0 < 0";
        }
        if (regionIndex[0] > iWidth){
            regionIndex[0] = iWidth;
            MITK_WARN << "Iterating Region Index0 > Width";
        }
        if (regionIndex[1] < 0){
            regionIndex[1] = 0;
            MITK_WARN << "Iterating Region Index1 < 0";
        }
        if (regionIndex[1] > iHeight){
            regionIndex[1] = iHeight;
            MITK_WARN << "Iterating Region Index1 > Height";
        }
        if (regionIndex[2] < 0){
            regionIndex[2] = 0;
            MITK_WARN << "Iterating Region Index2 < 0";
        }
        if (regionIndex[2] > iSlice){
            regionIndex[2] = iSlice;
            MITK_WARN << "Iterating Region Index2 > Slice";
        }

        // region size
        region.SetSize(regionSize);
        // region start
        region.SetIndex(regionIndex);
        radius[0] = 0;
        radius[1] = 0;
        radius[2] = 0;

        itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
        int nCounter = 0;
        while (!iterator.IsAtEnd())
        {
            ImageType::IndexType local_centerIndex = iterator.GetIndex();
            if (local_centerIndex[0] < 0 || local_centerIndex[1] < 0 || local_centerIndex[2] < 0 ||
                local_centerIndex[0] >= iWidth || local_centerIndex[1] >= iHeight || local_centerIndex[2] >= iSlice)
            {
                tsBOWOutput << 0.0 << "\n";
                ++iterator;
                nCounter++;
                continue;
            }

            int nTotalIn = 0;
            for (unsigned int iITER = 0; iITER < iterator.Size(); ++iITER)
            {
                ImageType::IndexType neighbour_index = iterator.GetIndex(iITER);

                bool IsInBounds;
                int neighborValue = iterator.GetPixel(iITER, IsInBounds);
                if (!IsInBounds){
                    neighborValue = -F_CT_WIDTH;
                }
                if (neighborValue > F_CT_WIDTH){
                    neighborValue = F_CT_WIDTH;
                }
                if (neighborValue < -F_CT_WIDTH){
                    neighborValue = -F_CT_WIDTH;
                }

                double dwNormalized = ((double)(neighborValue + F_CT_WIDTH)) / ((double)(F_CT_WIDTH * 2));
                tsBOWOutput << dwNormalized << "\n";
                nTotalIn++;
            }

            ++iterator;
            nCounter++;
        }
        qBOWOutput.close();
    }

    //-------------------------------------------------------------------------
    // output config
    QString q_bow_output_filename;
    q_bow_output_filename.sprintf("%s\\types_config.txt", strSavePath.toStdString().c_str());
    MITK_INFO << "Writing config to " << q_bow_output_filename.toStdString();

    QFile qBOWOutput(q_bow_output_filename);
    QTextStream tsBOWOutput(&qBOWOutput);
    if (qBOWOutput.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate)){
        tsBOWOutput.setRealNumberNotation(QTextStream::FixedNotation);
        // 1/2000 = 0.0005, so 4 digits after . is enough
        tsBOWOutput.setRealNumberPrecision(PIXEL_PRECISION);
    }
    else{
        MITK_ERROR << "Open " << q_bow_output_filename.toStdString() << " Failed";
        return;
    }
    
    tsBOWOutput << szAllNodules << "\n";
    for (int iNodule = 0; iNodule < szAllNodules; iNodule++)
    {
        int eOutputNoduleType = m_outputNoduleTypes[iNodule];
        tsBOWOutput << eOutputNoduleType << "\n";
    }
    qBOWOutput.close();
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::outputNonNoduleBOWTrainingData(QString strSavePath, mitk::Image::Pointer mitkTexture, mitk::Image::Pointer mitkNonNodule, int& szNonnodules)
{
    QString strQDebug;


    int iWidth = mitkTexture->GetDimension(0);
    int iHeight = mitkTexture->GetDimension(1);
    int iSlice = mitkTexture->GetDimension(2);

    ImageType::Pointer itkTexture;
    mitk::CastToItkImage(mitkTexture, itkTexture);

    ImageType::Pointer itkNonNodule;
    mitk::CastToItkImage(mitkNonNodule, itkNonNodule);


    // file header: sz_data szInput szOutput
    int szAllNonNodulesFromLabel = m_NonNoduleInfoFromLabel.size();
    int szNonNodulesFromXML = m_NonNoduleInfoFromXML.size();
//     int szData = N_NON_NODULE_D * N_NON_NODULE_D * N_NON_NODULE_D;
//     int szInput = N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D * N_LOCAL_WINDOW_D;
//    int szOutput = 1;
  
    szNonnodules = 0;
    for (int iNonNodule = 0; iNonNodule < szAllNonNodulesFromLabel; iNonNodule++)
    {
        szNonnodules++;

        QString tmpPath;
        tmpPath.sprintf("%s\\NONNODULE\\%d.txt", strSavePath.toStdString().c_str(), szNonnodules);

        QFile qTraining(tmpPath);
        QTextStream tsTraining(&qTraining);
        MITK_INFO << "Writing to " << tmpPath;

        if (qTraining.open(QIODevice::ReadWrite | QIODevice::Text | QIODevice::Truncate))
        {
            tsTraining.setRealNumberNotation(QTextStream::FixedNotation);
            tsTraining.setRealNumberPrecision(PIXEL_PRECISION);

            ImageType::RegionType non_nodule_region = m_NonNoduleInfoFromLabel[iNonNodule];
            ImageType::IndexType non_nodule_region_index = non_nodule_region.GetIndex();
            ImageType::SizeType non_nodule_region_size = non_nodule_region.GetSize();

            ImageType::IndexType center_index;
            center_index[0] = non_nodule_region_index[0] + non_nodule_region_size[0] / 2;
            center_index[1] = non_nodule_region_index[1] + non_nodule_region_size[1] / 2;
            center_index[2] = non_nodule_region_index[2] + non_nodule_region_size[2] / 2;

            // strQDebug.sprintf("processing non-nodule %d, center(%d,%d,%d), ", iNonNodule, center_index[0], center_index[1], center_index[2]);
            // MITK_INFO << strQDebug.toStdString();

            ImageType::SizeType			regionSize;
            ImageType::IndexType		regionIndex;
            ImageType::RegionType		region;
            ImageType::SizeType			radius;

            regionSize[0] = N_NODULE_BOW_D/*N_NODULE_D*/;
            regionSize[1] = N_NODULE_BOW_D;
            regionSize[2] = N_NODULE_BOW_D;
            regionIndex[0] = center_index[0] - N_NODULE_BOW_R;
            regionIndex[1] = center_index[1] - N_NODULE_BOW_R;
            regionIndex[2] = center_index[2] - N_NODULE_BOW_R;

            if (regionIndex[0] < 0){
                regionIndex[0] = 0;
                MITK_WARN << "Iterating Region Index0 < 0";
            }
            if (regionIndex[0] >= iWidth){
                regionIndex[0] = iWidth;
                MITK_WARN << "Iterating Region Index0 > Width";
            }
            if (regionIndex[1] < 0){
                regionIndex[1] = 0;
                MITK_WARN << "Iterating Region Index1 < 0";
            }
            if (regionIndex[1] >= iHeight){
                regionIndex[1] = iHeight;
                MITK_WARN << "Iterating Region Index1 > Height";
            }
            if (regionIndex[2] < 0){
                regionIndex[2] = 0;
                MITK_WARN << "Iterating Region Index2 < 0";
            }
            if (regionIndex[2] >= iSlice){
                regionIndex[2] = iSlice;
                MITK_WARN << "Iterating Region Index2 > Slice";
            }

            // region size
            region.SetSize(regionSize);
            // region start
            region.SetIndex(regionIndex);
            radius[0] = 0;
            radius[1] = 0;
            radius[2] = 0;

            itk::NeighborhoodIterator<ImageType> iterator(radius, itkTexture, region);
            int nCounter = 0;
            while (!iterator.IsAtEnd())
            {
                ImageType::IndexType local_centerIndex = iterator.GetIndex();

                // over boundary points
                if (local_centerIndex[0] >= iWidth || local_centerIndex[0] < 0 || 
                    local_centerIndex[1] >= iHeight|| local_centerIndex[1] < 0 || 
                    local_centerIndex[2] >= iSlice || local_centerIndex[2] < 0 )
                {
                    tsTraining << 0.0 << "\n";
                    ++iterator;
                    nCounter++;
                    continue;
                }

                unsigned int szIterrator = iterator.Size();
                int nTotalIn = 0;
                for (unsigned int iITER = 0; iITER < szIterrator; ++iITER)
                {
                    ImageType::IndexType neighbour_index = iterator.GetIndex(iITER);

                    bool IsInBounds;
                    int neighborValue = iterator.GetPixel(iITER, IsInBounds);
                    if (!IsInBounds){
                        neighborValue = -F_CT_WIDTH;
                    }
                    if (neighborValue > F_CT_WIDTH){
                        neighborValue = F_CT_WIDTH;
                    }
                    if (neighborValue < -F_CT_WIDTH){
                        neighborValue = -F_CT_WIDTH;
                    }

                    double dwNormalized = ((double)(neighborValue + F_CT_WIDTH)) / ((double)(F_CT_WIDTH * 2));
                    tsTraining << dwNormalized << "\n";
                    nTotalIn++;
                }

                ++iterator;
                nCounter++;
            }
            qTraining.close();
        }
        else
        {
            MITK_INFO << "Open" << tmpPath << " Failed";
            return;
        }
    }
}
//////////////////////////////////////////////////////////////////////////
bool LIDCXMLView::clearDirectory(QString path)
{
    QDir dir(path);
    QFileInfoList fileList;
    QFileInfo curFile;
    if (!dir.exists())  { return false; }//�ļ����棬�򷵻�false
    fileList = dir.entryInfoList(QDir::Dirs | QDir::Files
        | QDir::Readable | QDir::Writable
        | QDir::Hidden | QDir::NoDotAndDotDot
        , QDir::Name);
    while (fileList.size() > 0)//��������
    {
        int infoNum = fileList.size();
        for (int i = infoNum - 1; i >= 0; i--)
        {
            curFile = fileList[i];
            if (curFile.isFile())//������ļ���ɾ���ļ�
            {
                QFile fileTemp(curFile.filePath());
                fileTemp.remove();
                fileList.removeAt(i);
            }
            if (curFile.isDir())//������ļ���
            {
                QDir dirTemp(curFile.filePath());
                QFileInfoList fileList1 = dirTemp.entryInfoList(QDir::Dirs | QDir::Files
                    | QDir::Readable | QDir::Writable
                    | QDir::Hidden | QDir::NoDotAndDotDot
                    , QDir::Name);
                if (fileList1.size() == 0)//�²�û���ļ����ļ���
                {
                    dirTemp.rmdir(".");
                    fileList.removeAt(i);
                }
                else//�²����ļ��л��ļ�
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


////////////////////////////////////////////////////////////////////////
void LIDCXMLView::LoadWorkingDataSetConfFile()
{  // Ask the user for a list of files to open
    QString fileName = QFileDialog::getOpenFileName(NULL, "Open",
        "E:\\LIDC",
        tr("Confs (*.txt)"));

    if (!fileName.isEmpty())
    {
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return;

        m_Controls.cbWorkingDataSetList->clear();

        while (!file.atEnd()) {
            QByteArray line = file.readLine();
            QString data = QString::fromStdString(line.toStdString());
            // NOTE: PYTHON adds \n to line end, must remove it !!!!
            data.remove('\n');
            m_Controls.cbWorkingDataSetList->addItem(data);
        }
    }

}

////////////////////////////////////////////////////////////////////////
void LIDCXMLView::LoadCurrentWorkingListData()
{
    QString q_path = m_Controls.cbWorkingDataSetList->currentText();

    setLastOpenPath(q_path);

    QString q_path_1 = q_path + "\\ground.mitk";
    QString q_path_2 = q_path + "\\main.mitk";

    // prefer lighter data ...
    if (QFile::exists(q_path_1) == true){
        loadWorkingListData(q_path_1);
    }
    else if (QFile::exists(q_path_2) == true){
        loadWorkingListData(q_path_2);
    }
    else{
        MITK_ERROR << "Neihter ground nor main .mitk found in " << q_path;
        return;
    }

    //post operations
    QString tmp = getLastOpenPath();
    if (tmp != ""){
        m_WorkingDir = tmp;
    }
    m_WorkingDir += "/";
    m_Controls.lbWorkingDir->setText(m_WorkingDir);

    DoColorLabel();

}

////////////////////////////////////////////////////////////////////////
void LIDCXMLView::LoadNextWorkingListData()
{
    int total = m_Controls.cbWorkingDataSetList->count();
    int current_idx = m_Controls.cbWorkingDataSetList->currentIndex();
    if (current_idx >= total-1){
        MITK_WARN << "The Last Piece of Data !!";
        return;
    }

    m_Controls.cbWorkingDataSetList->setCurrentIndex(current_idx + 1);
    QString q_path = m_Controls.cbWorkingDataSetList->currentText();
    setLastOpenPath(q_path);

    QString q_path_1 = q_path + "\\ground.mitk";
    QString q_path_2 = q_path + "\\main.mitk";

    // prefer lighter data ...
    if (QFile::exists(q_path_1) == true){
        loadWorkingListData(q_path_1);
    }
    else if (QFile::exists(q_path_2) == true){
        loadWorkingListData(q_path_2);
    }
    else{
        MITK_ERROR << "Neihter ground nor main .mitk found in " << q_path;
        return;
    }

    //post operations
    QString tmp = getLastOpenPath();
    if (tmp != ""){
        m_WorkingDir = tmp;
    }
    m_WorkingDir += "/";
    m_Controls.lbWorkingDir->setText(m_WorkingDir);

    DoColorLabel();

}


//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::closeWorkspace()
{
    try
    {
        ctkPluginContext* context = mitk::PluginActivator::GetContext();
        mitk::IDataStorageService* dss = 0;
        ctkServiceReference dsRef = context->getServiceReference<mitk::IDataStorageService>();
        if (dsRef)
        {
            dss = context->getService<mitk::IDataStorageService>(dsRef);
        }

        if (!dss)
        {
            MITK_WARN << "IDataStorageService service not available. Unable to close project.";
            context->ungetService(dsRef);
            return;
        }

        mitk::IDataStorageReference::Pointer dataStorageRef = dss->GetActiveDataStorage();
        if (dataStorageRef.IsNull())
        {
            // No active data storage set (i.e. not editor with a DataStorageEditorInput is active).
            dataStorageRef = dss->GetDefaultDataStorage();
        }

        mitk::DataStorage::Pointer dataStorage = dataStorageRef->GetDataStorage();
        if (dataStorage.IsNull())
        {
            MITK_WARN << "No data storage available. Cannot close project.";
            return;
        }

        //check if we got the default datastorage and if there is anything else then helper object in the storage
        if (dataStorageRef->IsDefault() &&
            dataStorage->GetSubset(mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object", mitk::BoolProperty::New(true))))->empty())
        {
            return;
        }

        /* Remove everything */
        mitk::DataStorage::SetOfObjects::ConstPointer nodesToRemove = dataStorage->GetAll();
        dataStorage->Remove(nodesToRemove);

        // Remove the datastorage from the data storage service
        dss->RemoveDataStorageReference(dataStorageRef);
    }
    catch (std::exception& e)
    {
        MITK_ERROR << "Exception caught during closing project: " << e.what();
        QMessageBox::warning(NULL, "Error", QString("An error occurred during Close Project: %1").arg(e.what()));
    }


}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::loadWorkingListData(QString file_path)
{
    closeWorkspace();

    loadFiles(file_path);
}

//////////////////////////////////////////////////////////////////////////
void LIDCXMLView::loadFiles(const QString fileNames)
// //! [UtilLoadFiles]
{
    if (fileNames.isEmpty())
        return;

    mitk::IDataStorageReference::Pointer dataStorageRef;

    {
        ctkPluginContext* context = mitk::PluginActivator::GetContext();
        mitk::IDataStorageService* dss = 0;
        ctkServiceReference dsRef = context->getServiceReference<mitk::IDataStorageService>();
        if (dsRef)
        {
            dss = context->getService<mitk::IDataStorageService>(dsRef);
        }

        if (!dss)
        {
            QString msg = "IDataStorageService service not available. Unable to open files.";
            MITK_WARN << msg.toStdString();
            QMessageBox::warning(QApplication::activeWindow(), "Unable to open files", msg);
            return;
        }

        // Get the active data storage (or the default one, if none is active)
        dataStorageRef = dss->GetDataStorage();
        context->ungetService(dsRef);
    }

    mitk::DataStorage::Pointer dataStorage = dataStorageRef->GetDataStorage();

    // Do the actual work of loading the data into the data storage

    // Turn off ASSERT
#if defined(_MSC_VER) && !defined(NDEBUG) && defined(_DEBUG) && defined(_CRT_ERROR)
    int lastCrtReportType = _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_DEBUG);
#endif

    mitk::DataStorage::SetOfObjects::Pointer data;
    try
    {
        data = QmitkIOUtil::Load(fileNames, *dataStorage);
    }
    catch (const mitk::Exception& e)
    {
        MITK_INFO << e;
        return;
    }
    const bool dsmodified = !data->empty();

    // Set ASSERT status back to previous status.
#if defined(_MSC_VER) && !defined(NDEBUG) && defined(_DEBUG) && defined(_CRT_ERROR)
    if (lastCrtReportType)
        _CrtSetReportMode(_CRT_ASSERT, lastCrtReportType);
#endif

    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(dataStorage);
}


