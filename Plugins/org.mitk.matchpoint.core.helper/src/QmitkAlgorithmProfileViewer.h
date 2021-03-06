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


#ifndef QMITK_ALGORITHM_PROFILE_VIEWER_H
#define QMITK_ALGORITHM_PROFILE_VIEWER_H

#include <mapDeploymentDLLInfo.h>

#include "org_mitk_matchpoint_core_helper_Export.h"

#include "ui_QmitkAlgorithmProfileViewer.h"
#include <QWidget>

/**
 * \class QmitkAlgorithmProfileViewer
 * \brief Widget that views the information and profile of an algorithm stored in an DLLInfo object.
 */
class MITK_MATCHPOINT_CORE_HELPER_EXPORT QmitkAlgorithmProfileViewer : public QWidget, private Ui::QmitkAlgorithmProfileViewer
{
  Q_OBJECT

public:
  QmitkAlgorithmProfileViewer(QWidget* parent=0);

  /**
   * \brief Updates the widget according to the new info.
   * \param pointer to the info instance.
   * \remark The DLLInfo is not stored internally or as reference
   * to update the widget you must use the updateInfo() method.
   */
  void updateInfo(const map::deployment::DLLInfo* newInfo);

public Q_SLOTS:
  /**
    * \brief Slot that can be used to trigger updateInfo();
    */
  void OnInfoChanged(const map::deployment::DLLInfo* newInfo);
};

#endif // QmitkAlgorithmProfileViewer_H
