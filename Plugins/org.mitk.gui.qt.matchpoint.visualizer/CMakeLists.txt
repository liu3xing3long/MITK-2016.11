project(org_mitk_gui_qt_matchpoint_visualizer)

mitk_create_plugin(
  EXPORT_DIRECTIVE MITK_MITK_GUI_QT_REGVIS_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkMatchPointRegistration
)
