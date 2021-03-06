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

#ifndef QmitkFRPROptimizerViewWidgetHIncluded
#define QmitkFRPROptimizerViewWidgetHIncluded

#include "MitkRigidRegistrationUIExports.h"
#include "QmitkRigidRegistrationOptimizerGUIBase.h"
#include "ui_QmitkFRPROptimizerControls.h"
#include <itkArray.h>
#include <itkImage.h>
#include <itkObject.h>

/*!
* \brief Widget for rigid registration
*
* Displays options for rigid registration.
*/
class MITKRIGIDREGISTRATIONUI_EXPORT QmitkFRPROptimizerView : public QmitkRigidRegistrationOptimizerGUIBase
{
public:
  QmitkFRPROptimizerView(QWidget *parent = nullptr, Qt::WindowFlags f = nullptr);
  ~QmitkFRPROptimizerView();

  virtual mitk::OptimizerParameters::OptimizerType GetOptimizerType() override;

  virtual itk::Object::Pointer GetOptimizer() override;

  virtual itk::Array<double> GetOptimizerParameters() override;

  virtual void SetOptimizerParameters(itk::Array<double> metricValues) override;

  virtual void SetNumberOfTransformParameters(int transformParameters) override;

  virtual QString GetName() override;

  virtual void SetupUI(QWidget *parent) override;

protected:
  Ui::QmitkFRPROptimizerControls m_Controls;

  int m_NumberTransformParameters;
};

#endif
