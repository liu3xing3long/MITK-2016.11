project(org_mitk_gui_qt_classificationsegmentation)

mitk_create_plugin(
  EXPORT_DIRECTIVE CLASSIFICATIONSEGMENTATION_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkCLVigraRandomForest MitkCLUtilities
)
