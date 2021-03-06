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

#ifndef __mitkLabelSetImageToSurfaceThreadedFilter_H_
#define __mitkLabelSetImageToSurfaceThreadedFilter_H_

#include "mitkSegmentationSink.h"
#include "mitkSurface.h"
#include <MitkMultilabelExports.h>

namespace mitk
{
  class MITKMULTILABEL_EXPORT LabelSetImageToSurfaceThreadedFilter : public SegmentationSink
  {
  public:
    mitkClassMacro(LabelSetImageToSurfaceThreadedFilter, SegmentationSink)
      mitkAlgorithmNewMacro(LabelSetImageToSurfaceThreadedFilter);

  protected:
    LabelSetImageToSurfaceThreadedFilter(); // use smart pointers
    virtual ~LabelSetImageToSurfaceThreadedFilter();

    virtual void Initialize(const NonBlockingAlgorithm *other = nullptr) override;
    virtual bool ReadyToRun() override;

    virtual bool ThreadedUpdateFunction() override; // will be called from a thread after calling StartAlgorithm

    virtual void ThreadedUpdateSuccessful() override; // will be called from a thread after calling StartAlgorithm

  private:
    int m_RequestedLabel;
    Surface::Pointer m_Result;
  };

} // namespace

#endif // __mitkLabelSetImageToSurfaceThreadedFilter_H_
