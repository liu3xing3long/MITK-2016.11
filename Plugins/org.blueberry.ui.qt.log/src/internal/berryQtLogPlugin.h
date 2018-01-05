/*===================================================================

BlueBerry Platform

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef BERRYLOGPLUGIN_H_
#define BERRYLOGPLUGIN_H_

#include <ctkPluginActivator.h>
#include "berryQtPlatformLogModel.h"


namespace berry {

class QtLogPlugin : public QObject, public ctkPluginActivator
{
  Q_OBJECT
  Q_PLUGIN_METADATA(IID "org_blueberry_ui_qt_log")
  Q_INTERFACES(ctkPluginActivator)

public:

  QtLogPlugin();

  void start(ctkPluginContext* context) override;
  void stop(ctkPluginContext* context) override;

  static QtLogPlugin* GetInstance();

  QtPlatformLogModel* GetLogModel() const;

  ctkPluginContext* GetContext() const;

private:
  static QtLogPlugin* instance;
  QtPlatformLogModel* m_LogModel;
  ctkPluginContext* m_Context;

};

}

#endif /*BERRYLOGPLUGIN_H_*/
