project(org_mitk_gui_qt_matchpoint_algorithm_control)

mitk_create_plugin(
  EXPORT_DIRECTIVE MITK_MITK_GUI_QT_ALG_CTRL_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkMatchPointRegistration MitkMatchPointRegistrationUI
  PACKAGE_DEPENDS CTK|CTKCore
)
