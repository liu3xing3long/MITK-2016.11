project(org_mitk_gui_qt_matchpoint_algorithm_browser)

mitk_create_plugin(
  EXPORT_DIRECTIVE MITKGUI_QT_MATCHPOINT_ALG_BROWSER_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkMatchPointRegistration MitkMatchPointRegistrationUI
  PACKAGE_DEPENDS CTK|CTKCore
)
