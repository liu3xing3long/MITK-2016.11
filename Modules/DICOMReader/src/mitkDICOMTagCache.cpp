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

#include "mitkDICOMTagCache.h"

mitk::DICOMTagCache::DICOMTagCache()
:itk::Object()
{
}

mitk::DICOMTagCache::DICOMTagCache( const DICOMTagCache&)
:itk::Object()
{
}

mitk::DICOMTagCache::~DICOMTagCache()
{
}

void mitk::DICOMTagCache::SetInputFiles(const StringList& filenames)
{
  m_InputFilenames = filenames;
  this->Modified();
}
