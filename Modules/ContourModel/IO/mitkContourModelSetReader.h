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
#ifndef _MITK_ContourModelSetReader__H_
#define _MITK_ContourModelSetReader__H_

// MITK
#include <mitkAbstractFileReader.h>
#include <mitkBaseData.h>
#include <mitkContourModel.h>
#include <mitkContourModelSet.h>
#include <mitkMimeType.h>

#include <stack>
#include <string>
#include <tinyxml.h>
#include <vtkXMLParser.h>

namespace mitk
{
  /**
   * @brief
   * @ingroup MitkContourModelModule
  */
  class ContourModelSetReader : public mitk::AbstractFileReader
  {
  public:
    ContourModelSetReader(const ContourModelSetReader &other);

    ContourModelSetReader();

    virtual ~ContourModelSetReader();

    using AbstractFileReader::Read;
    virtual std::vector<itk::SmartPointer<BaseData>> Read() override;

  protected:
  private:
    ContourModelSetReader *Clone() const override;

    us::ServiceRegistration<mitk::IFileReader> m_ServiceReg;
  };
}

#endif