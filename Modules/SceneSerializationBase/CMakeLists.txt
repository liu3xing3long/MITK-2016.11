MITK_CREATE_MODULE(
  DEPENDS MitkCore
  PACKAGE_DEPENDS PUBLIC tinyxml
  WARNINGS_AS_ERRORS
)
add_subdirectory(test)

