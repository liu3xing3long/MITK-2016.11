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

#ifndef __BERRY_AND_EXPRESSION_H__
#define __BERRY_AND_EXPRESSION_H__

#include "berryCompositeExpression.h"

namespace berry {

class AndExpression : public CompositeExpression {

public:

  bool operator==(const Object* object) const override;

  EvaluationResult::ConstPointer Evaluate(IEvaluationContext* context) const override;
};

}  // namespace berry

#endif // __BERRY_AND_EXPRESSION_H__
