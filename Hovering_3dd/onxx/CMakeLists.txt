cmake_minimum_required(VERSION 3.10)

# Set the project name
project(ONNXRuntimeExample)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(ONNXRuntimeExample test.cpp)

# Manually specify the include directory for the ONNX Runtime headers
target_include_directories(ONNXRuntimeExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/onnxruntime-linux-x64-1.16.3/include)

# Set the path to the ONNX Runtime library
set(ONNX_RUNTIME_LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/onnxruntime-linux-x64-1.16.3/lib/libonnxruntime.so)

# Link against the ONNX Runtime library
target_link_libraries(ONNXRuntimeExample PRIVATE ${ONNX_RUNTIME_LIB_PATH})

