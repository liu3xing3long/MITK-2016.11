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

#include "mapDeploymentDLLHelper.h"
#include "mapContinuousElements.h"
#include "mapITKRigid3DClosedFormRegistrationAlgorithmTemplate.h"
#include "mapConfigure.h"

#include "MITK_Rigid_closedform_points_default_ProfileResource.h"

typedef map::core::continuous::Elements<3>::InternalPointSetType PointSetType;
typedef map::algorithm::boxed::ITKRigid3DClosedFormRegistrationAlgorithmTemplate<PointSetType, ::map::algorithm::MITK_Rigid_closedform_points_defaultUIDPolicy>::Type
AlgorithmType;

mapDeployAlgorithmMacro(AlgorithmType);
