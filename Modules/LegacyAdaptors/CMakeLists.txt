mitkFunctionCheckCompilerFlags("-Wno-deprecated-declarations" CMAKE_CXX_FLAGS)

MITK_CREATE_MODULE(
  DEPENDS MitkCore MitkIpPic
  PACKAGE_DEPENDS PRIVATE ITK|ITKIOImageBase
  )
