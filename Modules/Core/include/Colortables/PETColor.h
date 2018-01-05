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

#ifndef _lut_PETColor_h_
#define _lut_PETColor_h_

static const int PETColor[256][3] = {
  {0, 0, 0},       {0, 2, 1},       {0, 4, 3},       {0, 6, 5},       {0, 8, 7},       {0, 10, 9},
  {0, 12, 11},     {0, 14, 13},     {0, 16, 15},     {0, 18, 17},     {0, 20, 19},     {0, 22, 21},
  {0, 24, 23},     {0, 26, 25},     {0, 28, 27},     {0, 30, 29},     {0, 32, 31},     {0, 34, 33},
  {0, 36, 35},     {0, 38, 37},     {0, 40, 39},     {0, 42, 41},     {0, 44, 43},     {0, 46, 45},
  {0, 48, 47},     {0, 50, 49},     {0, 52, 51},     {0, 54, 53},     {0, 56, 55},     {0, 58, 57},
  {0, 60, 59},     {0, 62, 61},     {0, 65, 63},     {0, 67, 65},     {0, 69, 67},     {0, 71, 69},
  {0, 73, 71},     {0, 75, 73},     {0, 77, 75},     {0, 79, 77},     {0, 81, 79},     {0, 83, 81},
  {0, 85, 83},     {0, 87, 85},     {0, 89, 87},     {0, 91, 89},     {0, 93, 91},     {0, 95, 93},
  {0, 97, 95},     {0, 99, 97},     {0, 101, 99},    {0, 103, 101},   {0, 105, 103},   {0, 107, 105},
  {0, 109, 107},   {0, 111, 109},   {0, 113, 111},   {0, 115, 113},   {0, 117, 115},   {0, 119, 117},
  {0, 121, 119},   {0, 123, 121},   {0, 125, 123},   {0, 128, 125},   {1, 126, 127},   {3, 124, 129},
  {5, 122, 131},   {7, 120, 133}, // ok
  {136, 118, 135}, {138, 116, 137}, {140, 114, 139}, {142, 112, 141}, {144, 110, 143}, {146, 108, 145},
  {148, 106, 147}, {150, 104, 149}, {152, 102, 151}, {154, 100, 153}, {156, 98, 155},  {158, 96, 157},
  {160, 94, 159},  {162, 92, 161},  {164, 90, 163},  {166, 88, 165},  {168, 86, 167},  {170, 84, 169},
  {172, 82, 171},  {47, 80, 173},   {49, 78, 175},   {51, 76, 177},   {53, 74, 179},   {55, 72, 181},
  {57, 70, 183}, // ok
  {59, 68, 185},   {61, 66, 187},   {63, 64, 189},   {65, 63, 191},   {67, 61, 193},   {69, 59, 195},
  {71, 57, 197},   {73, 55, 199},   {75, 53, 201},   {77, 51, 203},   {79, 49, 205},   {81, 47, 207},
  {83, 45, 209},   {85, 43, 211},   {86, 41, 213},   {88, 39, 215},   {90, 37, 217},   {92, 35, 219},
  {94, 33, 221},   {96, 31, 223},   {98, 29, 225}, // ok
  {100, 27, 227},  {102, 25, 229},  {104, 23, 231},  {106, 21, 233},  {108, 19, 235},  {110, 17, 237},
  {112, 15, 239},  {114, 13, 241},  {116, 11, 243},  {118, 9, 245},   {120, 7, 247},   {122, 5, 249},
  {124, 3, 251},   {126, 1, 253},   {128, 0, 255}, // ok
  {130, 2, 252},   {132, 4, 248},   {134, 6, 244},   {136, 8, 240},   {138, 10, 236},  {140, 12, 232},
  {142, 14, 228},  {144, 16, 224},  {146, 18, 220},  {148, 20, 216},  {150, 22, 212},  {152, 24, 208},
  {154, 26, 204},  {156, 28, 200},  {158, 30, 196},  {160, 32, 192},  {162, 34, 188},  {164, 36, 184},
  {166, 38, 180},  {168, 40, 176},  {170, 42, 172},  {171, 44, 168},  {173, 46, 164},  {175, 48, 160},
  {177, 50, 156},  {179, 52, 152},  {181, 54, 148},  {183, 56, 144},  {185, 58, 140},  {187, 60, 136},
  {189, 62, 132},  {191, 64, 128}, // ok
  {193, 66, 124},  {195, 68, 120},  {197, 70, 116},  {199, 72, 112},  {201, 74, 108},  {203, 76, 104},
  {205, 78, 100},  {207, 80, 96},   {209, 82, 92},   {211, 84, 88},   {213, 86, 84},   {215, 88, 80},
  {217, 90, 76},   {219, 92, 72},   {221, 94, 68},   {223, 96, 64},   {225, 98, 60},   {227, 100, 56},
  {229, 102, 52},  {231, 104, 48},  {233, 106, 44},  {235, 108, 40},  {237, 110, 36},  {239, 112, 32},
  {241, 114, 28},  {243, 116, 24},  {245, 118, 20},  {247, 120, 16},  {249, 122, 12},  {251, 124, 8},
  {253, 126, 4},   {255, 128, 0}, // ok
  {255, 130, 4},   {255, 132, 8},   {255, 134, 12},  {255, 136, 16},  {255, 138, 20},  {255, 140, 24},
  {255, 142, 28},  {255, 144, 32},  {255, 146, 36},  {255, 148, 40},  {255, 150, 44},  {255, 152, 48},
  {255, 154, 50},  {255, 156, 54},  {255, 158, 58},  {255, 160, 62},  {255, 162, 64},  {255, 164, 68},
  {255, 166, 72},  {255, 168, 76},  {255, 170, 80},  {255, 172, 85},  {255, 174, 89},  {255, 176, 93},
  {255, 178, 97},  {255, 180, 101}, {255, 182, 105}, {255, 184, 109}, {255, 186, 113}, {255, 188, 117},
  {255, 190, 121}, {255, 192, 125}, {255, 194, 129}, {255, 196, 133}, {255, 198, 137}, {255, 200, 141},
  {255, 202, 145}, {255, 204, 149}, {255, 206, 153}, {255, 208, 157}, {255, 210, 161}, {255, 212, 165},
  {255, 214, 170}, {255, 216, 174}, {255, 218, 178}, {255, 220, 182}, {255, 222, 186}, {255, 224, 190},
  {255, 226, 194}, {255, 228, 198}, {255, 230, 202}, {255, 232, 206}, {255, 234, 210}, {255, 236, 214},
  {255, 238, 218}, {255, 240, 222}, {255, 242, 226}, {255, 244, 230}, {255, 246, 234}, {255, 248, 238},
  {255, 250, 242}, {255, 252, 250}, {255, 255, 255}};

#endif
