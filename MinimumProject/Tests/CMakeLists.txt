include_directories("${PROJECT_SOURCE_DIR}/Library")

add_executable(Test Test.cpp)
target_link_libraries(Test Class)