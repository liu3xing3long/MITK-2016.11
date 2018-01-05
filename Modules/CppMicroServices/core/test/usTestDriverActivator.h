/*=============================================================================

  Library: CppMicroServices

  Copyright (c) German Cancer Research Center,
    Division of Medical and Biological Informatics

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

=============================================================================*/

#ifndef USTESTDRIVERACTIVATOR_H
#define USTESTDRIVERACTIVATOR_H

#include <usModuleActivator.h>

US_BEGIN_NAMESPACE

class TestDriverActivator : public ModuleActivator
{
public:

  TestDriverActivator();

  static bool LoadCalled();

  void Load(ModuleContext*) override;

  void Unload(ModuleContext* ) override;

private:

  static TestDriverActivator* m_Instance;
  bool m_LoadCalled;
};

US_END_NAMESPACE

#endif // USTESTDRIVERACTIVATOR_H
