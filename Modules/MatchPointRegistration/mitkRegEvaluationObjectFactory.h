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

#ifndef __MITKRegEvaluationObjectFactory_h
#define __MITKRegEvaluationObjectFactory_h

#include <mitkCoreObjectFactory.h>
#include "MitkMatchPointRegistrationExports.h"

namespace mitk {

  /** Factory that registers everything (the mapper) needed for handling
   RegEvaluationObject instances in MITK.*/
  class RegEvaluationObjectFactory : public mitk::CoreObjectFactoryBase
  {
  public:
    mitkClassMacro(RegEvaluationObjectFactory,CoreObjectFactoryBase);
    itkNewMacro(RegEvaluationObjectFactory);

    ~RegEvaluationObjectFactory();

    virtual void SetDefaultProperties(mitk::DataNode* node);
    virtual const char* GetFileExtensions();
    virtual mitk::CoreObjectFactoryBase::MultimapType GetFileExtensionsMap();
    virtual const char* GetSaveFileExtensions();
    virtual mitk::CoreObjectFactoryBase::MultimapType GetSaveFileExtensionsMap();
    virtual mitk::Mapper::Pointer CreateMapper(mitk::DataNode* node, MapperSlotId slotId);
    void RegisterIOFactories();
  protected:
    std::string m_FileExtensions;
    RegEvaluationObjectFactory();
  };

}

#endif
