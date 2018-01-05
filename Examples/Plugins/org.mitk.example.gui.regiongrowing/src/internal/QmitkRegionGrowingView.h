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


#ifndef QmitkRegionGrowingView_h
#define QmitkRegionGrowingView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include "ui_QmitkRegionGrowingViewControls.h"

//! [includes]
#include "mitkPointSet.h"
#include "mitkIRenderWindowPartListener.h"
#include "TypeDefs.h"
#include "mitkImage.h"

#include <vector>

class QmitkPointListWidget;
class vtkRenderer;

// 
//-------------------------------------------------------------------------
class QmitkRegionGrowingView : public QmitkAbstractView, public mitk::IRenderWindowPartListener
{
    Q_OBJECT

public:
    typedef std::vector<ImageType::IndexType>	VIndexes;
    typedef std::vector<double>					VIndexValue;

    typedef std::vector<int>                    VecInt;
    typedef std::vector<VecInt>                 VVecInt;

public:
    struct SVisualizeNodule
    {
        ImageType::RegionType	nodule_region;
        VIndexes				nodule_indexes;
        VIndexValue				nodule_index_values;
        int						nodule_size;

        SVisualizeNodule()
        {
            nodule_size = -1;
        }

        SVisualizeNodule(const SVisualizeNodule& otherNodule)
        {
            nodule_region = otherNodule.nodule_region;
            nodule_indexes.insert(nodule_indexes.begin(), otherNodule.nodule_indexes.begin(), otherNodule.nodule_indexes.end());
            nodule_index_values.insert(nodule_index_values.begin(), otherNodule.nodule_index_values.begin(), otherNodule.nodule_index_values.end());
            nodule_size = otherNodule.nodule_size;
        }
        ~SVisualizeNodule()
        {
            nodule_indexes.clear();
            nodule_index_values.clear();
            nodule_size = -1;
        }
    };
    typedef std::vector<SVisualizeNodule>	VVisualizeNodules;
    enum EANNResultsType
    {
        EANNRT_ISO =  0,
        EANNRT_WALL,
        EANNRT_VESSEL,
        EANNRT_GGO,
        EANNRT_COUNT,
    };
    
public:
    static const std::string VIEW_ID;
    QmitkRegionGrowingView();

    protected slots:
    /// \brief Called when the user clicks the GUI button
    void DoRegionGrowing();
    void DoHessian();
    void DoThreshold();
    // 	void DoNormalization();
    // 	void DoMethodTest();
    // 	void DoMethodTest2();
    void DoSelectRegion();
    void DoFillHole();
    void DoExtractRealData();
    void DoErosion();
    void DoColorLabel();
    void DoColorLabel_2();
    void DoOutputForecast();
    // 	void DoCrop();
    void DoShowAnnResults();
    // 	void DoLoadIndexes();
    // 	void DoLoadIndexValues();
    // 	void DoRollingBall();
    void DoSelectWorkingDir();
    void DoThinning();
    void DoOSTUThreshold();
    void DoKMeans();
    void DoSelectLabel();
    void DoMedian();
    void DoDilation();
    void DoFillBorder();
    void DoTransferToMask();
    void DoChangeBGColor();
    void DoTakeScreenShots();
    void DoOutputScoringResults();
    void DoRefineRegionGrowing();
    void DoReduceFaces();
    void DoOutputBOWSamples();

    protected slots:
    void DoMethodTest1();
    void DoMethodTest2();
    void DoMethodTest3();
    void DoMethodTest4();

protected:
    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;

    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source,
        const QList<mitk::DataNode::Pointer>& nodes) override;

    //! [render-window-part-listener]
    void RenderWindowPartActivated(mitk::IRenderWindowPart* renderWindowPart) override;
    void RenderWindowPartDeactivated(mitk::IRenderWindowPart* renderWindowPart) override;
    //! [render-window-part-listener]

    Ui::QmitkRegionGrowingViewControls m_Controls;

private:
    // get current selected as mitk::image
    mitk::Image*	getSelectedImage();
    mitk::DataNode* getSelectedNode();
    // 	QString			getCurrentDataPath();
    QString         getLastOpenPath();

private:
    void ShowProcessedDataNode(mitk::Image::Pointer itkImage, std::string strNodeName, bool bSetDefaultProperties = false, mitk::DataNode::Pointer nodeParent = NULL);
    void ShowProcessedDataNode(ImageType::Pointer itkImage, std::string strNodeName, bool bSetDefaultProperties = false, mitk::DataNode::Pointer nodeParent = NULL);

private:
    bool LoadIndexValues(QString strFile, VVisualizeNodules& regions);
    bool LoadIndexes(QString strFile, ImageType::SizeType& image_size, VVisualizeNodules& regions);
    void SaveIndexes(QString strFile, ImageType::SizeType& image_size, VVisualizeNodules& regions);
    void SaveRegionParams(QString strFile, VVisualizeNodules& regions);
    bool LoadRegionParams(QString strFile, VVisualizeNodules& regions);
    void TakeScreenshot(vtkRenderer* renderer, unsigned int magnificationFactor, QString fileName);
    void SelectBackgroundColor();
    void OutputScoringResultImages(VVecInt& vVNodules);
    // clear all contents in this directory, root dir not included
    bool ClearDirectory(QString path);

private:
    bool normalizeRegionValues(VIndexValue& vIndexValues, VVisualizeNodules& vOutVisNodules);
    void collectRegions(BinaryImageType::Pointer itkBinaryImage, VVisualizeNodules& vOutput);
//     void buildNormalizedGaussian(double***&dwRET, int nGaussianR, bool bSigmaAdjustable=false);
    void buildNormalizedGaussian(double***&dwRET, int nGaussianRX, int nGaussianRY, int nGaussianRZ, bool bSigmaAdjustable);

    void showOneAnnResult(QString basepath, 
        EANNResultsType result_type, mitk::DataNode::Pointer baseNode, VecInt& vTotalNoduleIndex,
        bool bOutputGaussian = false, bool bOutputGGroundtruth = false, int nShownResults = 4);

    void ExportROIImage(QString outputFileName, 
        mitk::Image::Pointer& texture_image, mitk::Image::Pointer& result_image, 
        ImageType::RegionType& region_roi, ImageType::IndexType& offset);

private:
    void OutputBOWSampleData();

private:
    virtual void NodeAdded(const mitk::DataNode* node);

private:
    /// \brief This is the actual seed point data object
    mitk::PointSet::Pointer m_PointSet;
    QmitkPointListWidget* m_PointListWidget;
    ImageType::RegionType mCurrentRegion;

    QString m_WorkingDir;
    typedef std::vector<ImageType::RegionType> VLabeledRegions;
    VLabeledRegions m_vLabeledRegions;
    mitk::IRenderWindowPart* m_IRenderWindowPart;

    // screen shot
    QColor m_BackgroundColor;

    VVecInt m_vvTotalNoduleIndex;
};

#endif // QmitkRegionGrowingView_h
