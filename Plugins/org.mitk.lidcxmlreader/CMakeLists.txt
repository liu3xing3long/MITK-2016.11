project(org_mitk_lidcxmlreader)

link_libraries(libITKLabelMap-4.9.so)

mitk_create_plugin(
  EXPORT_DIRECTIVE LIDCXMLREADER_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt
)
