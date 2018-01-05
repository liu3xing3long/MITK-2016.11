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
#/* **************************************************************************
#  *                                                                          *
#  *     (C) Copyright Paul Mensonides 2002.
#  *     Distributed under the Boost Software License, Version 1.0. (See
#  *     accompanying file LICENSE_1_0.txt or copy at
#  *     http://www.boost.org/LICENSE_1_0.txt)
#  *                                                                          *
#  ************************************************************************** */
#
#/* See http://www.boost.org for most recent version. */
#
#ifndef MITK_PREPROCESSOR_SEQ_FOR_EACH_PRODUCT_HPP
#define MITK_PREPROCESSOR_SEQ_FOR_EACH_PRODUCT_HPP
#
#include "mitkPPArithmeticDec.h"
#include "mitkPPConfig.h"
#include "mitkPPControlIf.h"
#include "mitkPPRepetitionFor.h"
#include "mitkPPSeq.h"
#include "mitkPPSeqSize.h"
#include "mitkPPTupleElem.h"
#include "mitkPPTupleRem.h"
#
#/* MITK_PP_SEQ_FOR_EACH_PRODUCT */
#
#define MITK_PP_SEQ_FOR_EACH_PRODUCT(macro, sets) MITK_PP_SEQ_FOR_EACH_PRODUCT_E(MITK_PP_FOR, macro, sets)
#
#/* MITK_PP_SEQ_FOR_EACH_PRODUCT_R */
#
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_R(r, macro, sets) MITK_PP_SEQ_FOR_EACH_PRODUCT_E(MITK_PP_FOR_##r, macro, sets)
#
#if ~MITK_PP_CONFIG_FLAGS() & MITK_PP_CONFIG_EDG()
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_E(impl, macro, sets)                                                              \
  impl((MITK_PP_SEQ_HEAD(sets)(nil), MITK_PP_SEQ_TAIL(sets)(nil), (nil), macro),                                       \
       MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                                 \
       MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                                 \
       MITK_PP_SEQ_FOR_EACH_PRODUCT_M_0)
#else
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_E(impl, macro, sets) MITK_PP_SEQ_FOR_EACH_PRODUCT_E_I(impl, macro, sets)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_E_I(impl, macro, sets)                                                            \
  impl((MITK_PP_SEQ_HEAD(sets)(nil), MITK_PP_SEQ_TAIL(sets)(nil), (nil), macro),                                       \
       MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                                 \
       MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                                 \
       MITK_PP_SEQ_FOR_EACH_PRODUCT_M_0)
#endif
#
#if MITK_PP_CONFIG_FLAGS() & MITK_PP_CONFIG_STRICT()
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_P(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_P_I data
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_P_I(cset, rset, res, macro) MITK_PP_DEC(MITK_PP_SEQ_SIZE(cset))
#else
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_P(r, data) MITK_PP_DEC(MITK_PP_SEQ_SIZE(MITK_PP_TUPLE_ELEM(4, 0, data)))
#endif
#
#if ~MITK_PP_CONFIG_FLAGS() & MITK_PP_CONFIG_MWCC()
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_O(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_O_I data
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_O_I(cset, rset, res, macro) (MITK_PP_SEQ_TAIL(cset), rset, res, macro)
#else
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_O(r, data)                                                                        \
  (MITK_PP_SEQ_TAIL(MITK_PP_TUPLE_ELEM(4, 0, data)),                                                                   \
   MITK_PP_TUPLE_ELEM(4, 1, data),                                                                                     \
   MITK_PP_TUPLE_ELEM(4, 2, data),                                                                                     \
   MITK_PP_TUPLE_ELEM(4, 3, data))
#endif
#
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, i)                                                                        \
  MITK_PP_IF(MITK_PP_DEC(MITK_PP_SEQ_SIZE(MITK_PP_TUPLE_ELEM(4, 1, data))),                                            \
             MITK_PP_SEQ_FOR_EACH_PRODUCT_N_##i,                                                                       \
             MITK_PP_SEQ_FOR_EACH_PRODUCT_I)
#
#if ~MITK_PP_CONFIG_FLAGS() & MITK_PP_CONFIG_EDG()
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_I(r, data)                                                                        \
  MITK_PP_SEQ_FOR_EACH_PRODUCT_I_I(r,                                                                                  \
                                   MITK_PP_TUPLE_ELEM(4, 0, data),                                                     \
                                   MITK_PP_TUPLE_ELEM(4, 1, data),                                                     \
                                   MITK_PP_TUPLE_ELEM(4, 2, data),                                                     \
                                   MITK_PP_TUPLE_ELEM(4, 3, data))
#else
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_I(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_I_IM(r, MITK_PP_TUPLE_REM_4 data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_I_IM(r, im) MITK_PP_SEQ_FOR_EACH_PRODUCT_I_I(r, im)
#endif
#
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_I_I(r, cset, rset, res, macro)                                                    \
  macro(r, MITK_PP_SEQ_TAIL(res(MITK_PP_SEQ_HEAD(cset))))
#
#if ~MITK_PP_CONFIG_FLAGS() & MITK_PP_CONFIG_MWCC()
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data) MITK_PP_SEQ_FOR_EACH_PRODUCT_H_I data
#else
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data)                                                                           \
  MITK_PP_SEQ_FOR_EACH_PRODUCT_H_I(MITK_PP_TUPLE_ELEM(4, 0, data),                                                     \
                                   MITK_PP_TUPLE_ELEM(4, 1, data),                                                     \
                                   MITK_PP_TUPLE_ELEM(4, 2, data),                                                     \
                                   MITK_PP_TUPLE_ELEM(4, 3, data))
#endif
#
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_H_I(cset, rset, res, macro)                                                       \
  (MITK_PP_SEQ_HEAD(rset)(nil), MITK_PP_SEQ_TAIL(rset), res(MITK_PP_SEQ_HEAD(cset)), macro)
#
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_0(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 0)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_1(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 1)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_2(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 2)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_3(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 3)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_4(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 4)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_5(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 5)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_6(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 6)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_7(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 7)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_8(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 8)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_9(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 9)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_10(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 10)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_11(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 11)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_12(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 12)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_13(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 13)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_14(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 14)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_15(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 15)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_16(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 16)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_17(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 17)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_18(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 18)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_19(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 19)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_20(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 20)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_21(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 21)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_22(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 22)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_23(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 23)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_24(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 24)(r, data)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_M_25(r, data) MITK_PP_SEQ_FOR_EACH_PRODUCT_C(data, 25)(r, data)
#
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_0(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_1)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_1(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_2)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_2(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_3)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_3(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_4)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_4(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_5)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_5(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_6)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_6(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_7)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_7(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_8)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_8(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_9)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_9(r, data)                                                                      \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_10)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_10(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_11)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_11(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_12)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_12(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_13)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_13(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_14)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_14(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_15)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_15(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_16)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_16(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_17)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_17(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_18)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_18(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_19)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_19(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_20)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_20(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_21)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_21(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_22)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_22(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_23)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_23(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_24)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_24(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_25)
#define MITK_PP_SEQ_FOR_EACH_PRODUCT_N_25(r, data)                                                                     \
  MITK_PP_FOR_##r(MITK_PP_SEQ_FOR_EACH_PRODUCT_H(data),                                                                \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_P,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_O,                                                                      \
                  MITK_PP_SEQ_FOR_EACH_PRODUCT_M_26)
#
#endif
