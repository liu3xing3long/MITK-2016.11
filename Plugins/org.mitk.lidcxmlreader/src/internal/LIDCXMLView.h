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


#ifndef LIDCXMLView_h
#define LIDCXMLView_h

#include "berryISelectionListener.h"
#include "QmitkAbstractView.h"
#include "QStandardItemModel"
#include "ui_LIDCXMLViewControls.h"
#include "itkImage.h"
#include "QImage"
#include <map>
#include "mitkImage.h"

#include "Typedefs.h"
#include "NoduleInfoPerSlice.h"
#include "mitkIRenderWindowPartListener.h"

class vtkImageData;
class CMyDrawableNodule;


/**
  \brief LIDCXMLView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
  */
class LIDCXMLView : public QmitkAbstractView, public mitk::IRenderWindowPartListener
{
	// this is needed for all Qt objects that should have a Qt meta-object
	// (everything that derives from QObject and wants to have signal/slots)
	Q_OBJECT

public:
	LIDCXMLView();
	~LIDCXMLView();

public:
	static const std::string VIEW_ID;

protected:
    void RenderWindowPartActivated(mitk::IRenderWindowPart* renderWindowPart) override;
    void RenderWindowPartDeactivated(mitk::IRenderWindowPart* renderWindowPart) override;

public slots :
	/// \brief Do image processing
	void DoLoadAnnotations();
	void DoNormalization();
	void DoColorLabel();
	void DoOutputTraining();
    void DoOutputBOWTraining();
    void DoSelectWorkingDir();
    void DoResampleNonNodule();
    void DoResampleAll();
    void DoResampleGroundAndTexture();

    void DoExternalTransform();
    void DoSelectExternalNodule();
    void DoOutputNonExternalNodule();

    /// \brief Load selected annotation xml file
	void LoadAnnotationXML(QString filename, VVSNoduleInfo& vvNoduleInfo);

	/// show image function
	void showPriorImage();
	void showNextImage();

    // 
    void LoadWorkingDataSetConfFile();
    void LoadCurrentWorkingListData();
    void LoadNextWorkingListData();

    ///
    void outputNoduleTypeChanged(int iType);

// 	/// 
// 	void saveSpecificImage();

protected:

	virtual void CreateQtPartControl(QWidget *parent) override;

	virtual void SetFocus() override;

	/// \brief called by QmitkFunctionality when DataManager's selection has changed
	virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source,
		const QList<mitk::DataNode::Pointer>& nodes) override;

	Ui::LIDCXMLViewControls m_Controls;

private:
	CMyDrawableNodule*		buildDrawableLabel(mitk::Image* mitkImage, CNoduleInfoPerSlice* sThisExpThisInfo);
    void					LabelImageUsingXML(mitk::Image* mitkImage, VVSNoduleInfo& vvNoduleInfo);
    void					collectInfoFromXML(mitk::Image* mitkImage, VVSNoduleInfo& vvNoduleInfo);
    void                    collectNonNoduleInfo(VVSNoduleInfo& vvNoduleInfo, VSNoduleInfo& vOutNonNodules);
	void					destroyContainers();
    void                    destroyUIComponents();
	void					fillDataGrid(CNoduleInfoPerSlice* pNoduleInSlice);
    void                    closeWorkspace();
    void                    loadWorkingListData(QString file_path);
    void                    loadFiles(const QString fileNames);


private:
    void					outputNoduleTrainingData(QString strSavePath, mitk::Image::Pointer mitkTexture, mitk::Image::Pointer mituGroundTruth);
    void					outputNonNoduleTrainingData(QString strSavePath, mitk::Image::Pointer mitkTexture, mitk::Image::Pointer mitkNonNodule);
    void                    outputNoduleBOWTrainingData(QString strSavePath, mitk::Image::Pointer mitkTexture, mitk::Image::Pointer mitkGroundTruth, int* szNodulesThisType);
    void					outputNonNoduleBOWTrainingData(QString strSavePath, mitk::Image::Pointer mitkTexture, mitk::Image::Pointer mitkNonNodule, int& szNonnodules);
    
    mitk::Image::Pointer    normalizeImage(mitk::Image::Pointer mitkImage);
    void                    buildNormalizedGaussian(double***&dwRET, int nGaussianR);

private:
	mitk::Image::Pointer	getSelectedImage();
	mitk::DataNode::Pointer getSelectedNode();
	void ShowProcessedDataNode(mitk::Image::Pointer itkImage, std::string strNodeName, bool bSetDefaultProperties = false);
	void ShowProcessedDataNode(ImageType::Pointer itkImage, std::string strNodeName, bool bSetDefaultProperties = false);
    QString getLastOpenPath();
    bool setLastOpenPath(QString q_path);
    bool clearDirectory(QString path);

private:
// 	VVSNoduleInfo m_NoduleInExpert;
	QStandardItemModel  *m_ReportModel;
    mitk::IRenderWindowPart* m_IRenderWindowPart;

private:
    // used from xml
    VSNoduleInfo m_NonNoduleInfoFromXML;
	typedef vector<CNoduleInfoPerSlice*> VNoduleInfoPerSlice;
	VNoduleInfoPerSlice m_NodulesInSlice;
	int iCurrengShownImage;
	typedef vector<ImageType::RegionType> VNoduleCenter;
	VNoduleCenter m_vNoduleRegions;
    // used from direct labeling
    VNoduleCenter m_NonNoduleInfoFromLabel;
    QString m_WorkingDir; 
    int*   m_outputNoduleTypes;

    //////////////////////////////////////////////////////////////////////////
    VNoduleCenter m_vExternalRegions;
    typedef vector<int> VInt;
    VInt m_vExternalNoduleIndex;

    mitk::Image::Pointer m_bgImage;

};
#endif // LIDCXMLView_h
