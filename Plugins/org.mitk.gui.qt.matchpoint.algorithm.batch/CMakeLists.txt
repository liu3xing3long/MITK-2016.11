project(org_mitk_gui_qt_matchpoint_algorithm_batch)

mitk_create_plugin(
  EXPORT_DIRECTIVE MITK_GUI_QT_MATCHPOINT_ALG_BATCH_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkMatchPointRegistration MitkMatchPointRegistrationUI
  PACKAGE_DEPENDS CTK|CTKCore
)