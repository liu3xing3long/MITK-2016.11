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

#ifndef mitkBinaryThresholdTool_h_Included
#define mitkBinaryThresholdTool_h_Included

#include "mitkAutoSegmentationTool.h"
#include "mitkCommon.h"
#include "mitkDataNode.h"
#include <MitkSegmentationExports.h>

#include <itkImage.h>

namespace us
{
  class ModuleResource;
}

namespace mitk
{
  /**
  \brief Calculates the segmented volumes for binary images.

  \ingroup ToolManagerEtAl
  \sa mitk::Tool
  \sa QmitkInteractiveSegmentation

  Last contributor: $Author$
  */
  class MITKSEGMENTATION_EXPORT BinaryThresholdTool : public AutoSegmentationTool
  {
  public:
    Message3<double, double, bool> IntervalBordersChanged;
    Message1<double> ThresholdingValueChanged;

    mitkClassMacro(BinaryThresholdTool, AutoSegmentationTool);
    itkFactorylessNewMacro(Self) itkCloneMacro(Self)

      virtual const char **GetXPM() const override;
    us::ModuleResource GetIconResource() const override;
    virtual const char *GetName() const override;

    virtual void Activated() override;
    virtual void Deactivated() override;

    virtual void SetThresholdValue(double value);
    virtual void AcceptCurrentThresholdValue();
    virtual void CancelThresholding();

  protected:
    BinaryThresholdTool(); // purposely hidden
    virtual ~BinaryThresholdTool();

    void SetupPreviewNode();

    void CreateNewSegmentationFromThreshold(DataNode *node);

    void OnRoiDataChanged();
    void UpdatePreview();

    template <typename TPixel, unsigned int VImageDimension>
    void ITKThresholding(itk::Image<TPixel, VImageDimension> *originalImage,
                         mitk::Image *segmentation,
                         double thresholdValue,
                         unsigned int timeStep);
    template <typename TPixel, unsigned int VImageDimension>
    void ITKThresholdingOldBinary(itk::Image<TPixel, VImageDimension> *originalImage,
                                  mitk::Image *segmentation,
                                  double thresholdValue,
                                  unsigned int timeStep);

    DataNode::Pointer m_ThresholdFeedbackNode;
    DataNode::Pointer m_OriginalImageNode;
    DataNode::Pointer m_NodeForThresholding;

    double m_SensibleMinimumThresholdValue;
    double m_SensibleMaximumThresholdValue;
    double m_CurrentThresholdValue;
    bool m_IsFloatImage;

    bool m_IsOldBinary = false;
  };

} // namespace

#endif
