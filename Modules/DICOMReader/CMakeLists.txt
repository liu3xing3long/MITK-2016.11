MITK_CREATE_MODULE(
  DEPENDS MitkCore
  PACKAGE_DEPENDS
    PUBLIC tinyxml
    PRIVATE ITK|ITKIOImageBase+ITKIOGDCM DCMTK
)

add_subdirectory(test)
