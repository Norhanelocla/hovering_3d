#include <iostream>
#include <vector>
#include "onnxruntime-linux-x64-1.16.3/include/onnxruntime_cxx_api.h"

int main() {
    // Initialize ONNX Runtime
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "test");
    Ort::SessionOptions session_options;

    // Path to your ONNX model
    const char* model_path = "student.onnx";

    try {
        // Load ONNX model
        Ort::Session session(env, model_path, session_options);

        // Use the specific input and output node names from Netron
        const char* input_node_name = "input_8"; // The input node name
        const char* output_node_name = "tf.math.multiply_3"; // The output node name

        // Define input tensor shape: {1, 55}
        std::vector<int64_t> input_tensor_shape = {1, 55};

    // Create input tensor with example data
         std::vector<float> input_tensor_values = {
        -2.0025043, 3.5976875, 2.2145584, -0.79325545, 1.7449694, 0.61773235,
        0.0, 0.0, 0.0, -1.5, -0.26175457, 0.96513444, 0.0, 0.42433208, 0.9055067,
        0.0, 0.5450073, -0.8384313, 0.09647848, -0.11284962, 0.988917, 0.09706222,
        0.07902301, 0.0, 0.0, -0.42105064, -3.8156188, -0.18604296, -1.9376003,
        3.9941428, -4.241599, -4.605695, -0.36507428, 4.404779, -3.60964, 4.5535803,
        -4.8436255, 3.9975975, -3.8799942, 3.6966925, -2.3749971, 4.652264, -4.653452,
        -4.3407617, -2.5666244, 4.6044726, -4.638184, 3.7637677, -4.7817354, 2.7948966,
        -4.3707013, 4.137228, -4.978357, 2.336006, -3.3592284
    };

        Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
        Ort::Value input_tensor = Ort::Value::CreateTensor<float>(memory_info, input_tensor_values.data(), input_tensor_values.size(), input_tensor_shape.data(), input_tensor_shape.size());

        // Run the model
        auto output_tensors = session.Run(Ort::RunOptions{nullptr}, &input_node_name, &input_tensor, 1, &output_node_name, 1);

        // Process and print the output
        float* output_data = output_tensors.front().GetTensorMutableData<float>();
        size_t output_data_length = output_tensors.front().GetTensorTypeAndShapeInfo().GetElementCount();

        std::cout << "Output: ";
        for (size_t i = 0; i < output_data_length; ++i) {
            std::cout << output_data[i] << " ";
        }
        std::cout << std::endl;

    } catch (const Ort::Exception& exception) {
        std::cerr << "ONNX Runtime error: " << exception.what() << std::endl;
        return -1;
    }

    return 0;
}
