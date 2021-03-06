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

#include "org_mitk_example_gui_customviewer_Activator.h"

#include "CustomViewer.h"
#include "DicomPerspective.h"
#include "ViewerPerspective.h"

ctkPluginContext *org_mitk_example_gui_customviewer_Activator::PluginContext = nullptr;

void org_mitk_example_gui_customviewer_Activator::start(ctkPluginContext *context)
{
  BERRY_REGISTER_EXTENSION_CLASS(CustomViewer, context)
  BERRY_REGISTER_EXTENSION_CLASS(ViewerPerspective, context)
  BERRY_REGISTER_EXTENSION_CLASS(DicomPerspective, context)
  PluginContext = context;
}

void org_mitk_example_gui_customviewer_Activator::stop(ctkPluginContext *context)
{
  Q_UNUSED(context)

  PluginContext = nullptr;
}

ctkPluginContext *org_mitk_example_gui_customviewer_Activator::GetPluginContext()
{
  return PluginContext;
}
