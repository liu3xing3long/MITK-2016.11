MITK_CREATE_MODULE( IGTUIConfigurationWidgetRegistry
  INCLUDE_DIRS
    PRIVATE src/ConfigurationWidgetRegistry
  DEPENDS PRIVATE MitkIGTUI
  WARNINGS_AS_ERRORS
  AUTOLOAD_WITH MitkIGTUI
)