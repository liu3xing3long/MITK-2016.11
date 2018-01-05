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


#include "org_mitk_lidcxmlreader_Activator.h"

#include <QtPlugin>

#include "LIDCXMLView.h"

namespace mitk {

    org_mitk_lidcxmlreader_Activator* org_mitk_lidcxmlreader_Activator::m_Instance = nullptr;
    ctkPluginContext* org_mitk_lidcxmlreader_Activator::m_Context = nullptr;

    void org_mitk_lidcxmlreader_Activator::start(ctkPluginContext* context)
    {
        this->m_Instance = this;
        this->m_Context = context;
        
        BERRY_REGISTER_EXTENSION_CLASS(LIDCXMLView, context)

        this->m_PrefServiceTracker.reset(new ctkServiceTracker<berry::IPreferencesService*>(context));
        this->m_PrefServiceTracker->open();
    }

    void org_mitk_lidcxmlreader_Activator::stop(ctkPluginContext* context)
    {
        Q_UNUSED(context)
        
        this->m_PrefServiceTracker.reset();
        this->m_Context = nullptr;
        this->m_Instance = nullptr;
    }
    ctkPluginContext* org_mitk_lidcxmlreader_Activator::GetContext()
    {
        return m_Context;
    }

    org_mitk_lidcxmlreader_Activator *org_mitk_lidcxmlreader_Activator::GetInstance()
    {
        return m_Instance;
    }

    berry::IPreferencesService* org_mitk_lidcxmlreader_Activator::GetPreferencesService()
    {
        return m_PrefServiceTracker->getService();
    }
}

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
Q_EXPORT_PLUGIN2(org_mitk_lidcxmlreader, mitk::org_mitk_lidcxmlreader_Activator)
#endif
