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

#ifndef IMAGECHANNELSELECTOR_H_HEADER_INCLUDED_C1E4F4E7
#define IMAGECHANNELSELECTOR_H_HEADER_INCLUDED_C1E4F4E7

#include "mitkSubImageSelector.h"
#include <MitkCoreExports.h>

namespace mitk
{
  //##Documentation
  //## @brief Provides access to a channel of the input image
  //##
  //## If the input is generated by a ProcessObject, only the required data is
  //## requested.
  //## @ingroup Process
  class MITKCORE_EXPORT ImageChannelSelector : public SubImageSelector
  {
  public:
    mitkClassMacro(ImageChannelSelector, SubImageSelector);

    itkFactorylessNewMacro(Self) itkCloneMacro(Self)

      itkGetConstMacro(ChannelNr, int);
    itkSetMacro(ChannelNr, int);

  protected:
    ImageChannelSelector();

    virtual ~ImageChannelSelector();

    virtual void GenerateOutputInformation() override;

    virtual void GenerateInputRequestedRegion() override;

    virtual void GenerateData() override;

    int m_ChannelNr;
  };

} // namespace mitk

#endif /* IMAGECHANNELSELECTOR_H_HEADER_INCLUDED_C1E4F4E7 */
