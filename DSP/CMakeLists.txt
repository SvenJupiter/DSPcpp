# define Library "DSPcpp"
add_library(DSPcpp SHARED "")
target_include_directories(DSPcpp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(DSPcpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(DSPcpp PRIVATE -D BUILD_SHARDED_DSP_LIB=1)
target_compile_definitions(DSPcpp INTERFACE -D USE_SHARDED_DSP_LIB=1)
target_compile_options(DSPcpp PRIVATE -Werror -Wall -Wextra)
# target_compile_options(DSPcpp PRIVATE -Wall)


# add library subfolders
cmake_policy(SET CMP0076 NEW)
add_subdirectory(src)
