project(org_blueberry_core_runtime_tests)

mitk_create_plugin(
  EXPORT_DIRECTIVE BERRY_RUNTIME_TESTS
  TEST_PLUGIN
)

target_link_libraries(${PROJECT_NAME} optimized CppUnit debug CppUnitd)

MACRO_TEST_PLUGIN()
