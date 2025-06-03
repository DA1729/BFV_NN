#include "./bfv_lib/bfv_lib.h"
#include "./libs/mod_tools.h"
#include <iostream>
#include <random>
#include <vector>
#include <cmath>

class bfv_neural_network{
    private:
        bfv* evaluator;
        int input_size;
        int output_size;

        // encrypted weight and bias polynomials
        std::vector<std::vector<std::vector<poly_utils_1>>> enc_weights;
        std::vector<std::vector<poly_utils_1>> enc_biases;

        // helper functions
        std::vector<poly_utils_1> enc_zero(){
            poly_utils_1 zero_poly = evaluator -> int_encode(0);
            return evaluator->encryption(zero_poly);
        }

        std::vector<poly_utils_1> add_encrypted(const std::vector<poly_utils_1>& ct1, const std::vector<poly_utils_1>& ct2){
            return evaluator-> homomorphic_addition(ct1, ct2);
        }

        std::vector<poly_utils_1> mul_encrypted(const std::vector<poly_utils_1>& ct1, const std::vector<poly_utils_1>& ct2){
            std::vector<poly_utils_1> result = evaluator->homomorphic_multiplication(ct1, ct2);

            return evaluator->relinearization_1(result);
        }

    public:
        bfv_neural_network(bfv* evaluator, int input_size, int output_size)
            : evaluator(evaluator), input_size(input_size), output_size(output_size) {

                // initializing weights and biases
                enc_weights.resize(output_size);
                enc_biases.resize(output_size);

                for (size_t i = 0; i < output_size; i++){
                    enc_weights[i].resize(input_size);
                }
            }

        void init_network(){
            std::random_device rd;
            std::mt19937 g(rd());
            std::uniform_int_distribution<int64_t> weight_d(-10, 10);
            std::uniform_int_distribution<int64_t> bias_d(-5, 5);

            std::cout << "initializing network parameters" << std::endl;
            std::cout << "input size: " << input_size << std::endl;
            std::cout << "outptu size: " << output_size << std::endl;

            // generating and ecrypting weights

            for (size_t i = 0; i < output_size; i++){
                for (size_t j = 0; j < input_size; j++){
                    int64_t weight = weight_d(g);
                    poly_utils_1 weight_poly = evaluator->int_encode(weight);
                    enc_weights[i][j] = evaluator->encryption(weight_poly);
                }

                int64_t bias = bias_d(g);
                poly_utils_1 bias_poly = evaluator->int_encode(bias);
                enc_biases[i] = evaluator->encryption(bias_poly);
            }

            std::cout << "network parameters encrypted" << std::endl;
            std::cout << std::endl;
        }

        // below function is for testing purposes
        void set_weights(const std::vector<std::vector<int64_t>>& weights, const std::vector<int64_t>& biases){
            for (size_t i = 0; i < output_size; i++){
                for (size_t j = 0; j < input_size; j++){
                       poly_utils_1 weight_poly = evaluator->int_encode(weights[i][j]);
                       enc_weights[i][j] = evaluator->encryption(weight_poly);
                }

                poly_utils_1 bias_poly = evaluator -> int_encode(biases[i]);
                enc_biases[i] = evaluator->encryption(bias_poly);
            }

            // printing out weights and biases matrices
            std::cout<<"weights matrix" << std::endl;
            for (size_t i = 0; i < output_size; i++){
                std::cout << "[";
                for (size_t j = 0; j < input_size; j++){
                    std::cout<<weights[i][j];
                    if (j < input_size - 1) std::cout << ", ";
                }
                std::cout << "]" << std::endl;
            }

            std::cout << "biases" << std::endl;
            std::cout << "[";
            for (size_t i = 0; i < output_size; i++){
                std::cout<<biases[i];
                if (i < output_size - 1) std::cout<<", ";
            }

            std::cout << "]" << std::endl;
            std::cout << std::endl;
        }


        // forward pass through the layer
        std::vector<std::vector<poly_utils_1>> forward_pass(const std::vector<std::vector<poly_utils_1>>& enc_input){
            std::cout<< "---performing forward pass on encrypted input---" << std::endl;

            std::vector<std::vector<poly_utils_1>> output(output_size);

            for (size_t i = 0; i < output_size; i++){
                std::cout << "computing output of neuron " << i << std::endl;

                output[i] = enc_biases[i];

                for (size_t j = 0; j < input_size; j++){
                    std::vector<poly_utils_1> weighted_in = mul_encrypted(enc_weights[i][j], enc_input[j]);

                    output[i] = add_encrypted(output[i], weighted_in);
                }
            }

            std::cout<<"forward pass completed" << std::endl;
            std::cout<<std::endl;

            return output;
        }

        std::vector<int64_t> decrypt_output(const std::vector<std::vector<poly_utils_1>>& enc_output){
            std::vector<int64_t> dec_output(output_size);

            std::cout<<"---decrypting neural network output---" << std::endl;

            for (size_t i = 0; i < output_size; i++){
                poly_utils_1 dec_poly = evaluator->decryption(enc_output[i]);
                dec_output[i] = evaluator->int_decode(dec_poly);
                std::cout<<"output" << i << ": " << dec_output[i] << std::endl;
            }

            std::cout << std::endl;
            return dec_output;
        }

        std::vector<int64_t> expected_output(const std::vector<int64_t>& inputs,
                                                const std::vector<std::vector<int64_t>>& weights,
                                                const std::vector<int64_t>& biases){
            std::vector<int64_t> expected_out(output_size);

            std::cout<<"---computing the expected output---" << std::endl;

            for (size_t i = 0; i < output_size; i++){
                expected_out[i] = biases[i];
                for (size_t j = 0; j < input_size; j++){
                    expected_out[i] += weights[i][j] * inputs[j];
                }

                std::cout << "expected output " << i << ": " << expected_out[i] << std::endl;
            }

            std::cout<<std::endl;

            return expected_out;
        }

};


int main(void){
    std::cout << "=== BFV encrypted single layer neural network demo ===" << std::endl;
    std::cout << std::endl;

    int64_t t = 16;
    int64_t n = 1024;
    int64_t q = 132120577;
    int64_t psi = 73993;
    
    int64_t psiv = mod_inv(psi, q);
    int64_t w = mod_pow(psi, 2, q);
    int64_t wv = mod_inv(w, q);
    
    double mu = 0.0;
    double sigma = 0.5 * 3.2;
    
    // generating polynomial arithmetic tables
    std::vector<int64_t> w_table(n, 1);
    std::vector<int64_t> wv_table(n, 1);
    std::vector<int64_t> psi_table(n, 1);
    std::vector<int64_t> psiv_table(n, 1);
    
    for (int64_t i = 1; i < n; i++) {
        w_table[i] = (w_table[i-1] * w) % q;
        wv_table[i] = (wv_table[i-1] * wv) % q;
        psi_table[i] = (psi_table[i-1] * psi) % q;
        psiv_table[i] = (psiv_table[i-1] * psiv) % q;
    }
    
    ntt_params qnp;
    qnp.w = w_table;
    qnp.w_inv = wv_table;
    qnp.psi = psi_table;
    qnp.psi_inv = psiv_table;
    
    // bfv evaluator
    bfv Evaluator(n, q, t, mu, sigma, qnp);

    // generating keys
    std::cout << "---generating keys---" << std::endl;
    Evaluator.secret_key_gen();
    Evaluator.public_key_gen();
    Evaluator.eval_key_gen_1();
    std::cout << "---BFV keys generated successfully---" << std::endl;

    // neural network parameters
    int input_size = 3;
    int output_size = 2;

    // creating neural network
    bfv_neural_network nn(&Evaluator, input_size, output_size);

    // custom weights and biases for testing purposes

    std::vector<std::vector<int64_t>> weights = {
        {-1, 2, -3},
        {-2, 1, 1}
    };

    std::vector<int64_t> biases = {1, 2};

    nn.set_weights(weights, biases);

    // input
    std::vector<int64_t> input = {1, 2, 3};

    std::cout << "---test input:---" << std::endl;
    std::cout << "[";
    for (size_t i = 0; i < input_size; i++){
        std::cout<<input[i];
        if (i < input_size-1) std::cout << ", ";
    }

    std::cout << "]" << std::endl;
    std::cout << std::endl;

    // encrypting input
    std::cout<<"---encrypting the input data---" << std::endl;
    std::vector<std::vector<poly_utils_1>> enc_input(input_size);
    for (size_t i = 0; i < input_size; i++){
        poly_utils_1 input_poly = Evaluator.int_encode(input[i]);
        enc_input[i] = Evaluator.encryption(input_poly);
    }

    std::cout << "---input encrypted successfully---" << std::endl;
    std::cout << std::endl;

    // forward pass
    std::vector<std::vector<poly_utils_1>> enc_output = nn.forward_pass(enc_input);

    // decryption and verification
    std::vector<int64_t> dec_output = nn.decrypt_output(enc_output);
    std::vector<int64_t> expected_output = nn.expected_output(input, weights, biases);

    // verification
    std::cout << "---verification---" << std::endl;
    bool all_correct = true;

    for (size_t i = 0; i < output_size; i++){
        bool correct = (dec_output[i] == expected_output[i]);
        std::cout << "* output " << i << ": encrypted=" << dec_output[i] 
                  << ", expected=" << expected_output[i] 
                  << " -> " << (correct ? "CORRECT" : "INCORRECT") << std::endl;
        if (!correct) all_correct = false;  
    }
    std::cout << std::endl;
    if (all_correct) {
        std::cout << "SUCCESS: encrypted neural network computation matches expected results!" << std::endl;
    } else {
        std::cout << "ERROR: encrypted neural network computation has errors." << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << "=== demo complete ===" << std::endl;

    return 0;

}
