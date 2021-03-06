set(SRC_CPP_FILES
)

set(INTERNAL_CPP_FILES
  ChangeTextToLowerCase.cpp
  ChangeTextToUpperCase.cpp
  org_mitk_example_gui_extensionpointcontribution_Activator.cpp
)

set(MOC_H_FILES
  src/internal/ChangeTextToLowerCase.h
  src/internal/ChangeTextToUpperCase.h
  src/internal/org_mitk_example_gui_extensionpointcontribution_Activator.h
)

set(CACHED_RESOURCE_FILES
  plugin.xml
)

set(CPP_FILES )

foreach(file ${SRC_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/${file})
endforeach(file ${SRC_CPP_FILES})

foreach(file ${INTERNAL_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/internal/${file})
endforeach(file ${INTERNAL_CPP_FILES})
