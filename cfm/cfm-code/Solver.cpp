//
// Created by feiw on 15/04/18.
//

#include "Solver.h"


Sgd::Sgd(float learning_rate){
    this->learning_rate = learning_rate;
}

void Sgd::adjustWeights(std::vector<float> &grads,
                        std::set<unsigned int> &used_idxs,
                        boost::shared_ptr<Param> param) {

    for (auto &used_idx: used_idxs)
        (*(param->getWeightsPtr()))[used_idx] += learning_rate * grads[used_idx];
}

Momentum::Momentum(unsigned int length, float learning_rate, float momentum) {
    prev_v = std::vector<float>(length, 0.0);
    this->learning_rate = learning_rate;
    this->momentum = momentum;
}

void Momentum::adjustWeights(std::vector<float> &grads,
                             std::set<unsigned int> &used_idxs,
                             boost::shared_ptr<Param> param) {

    for (auto &used_idx: used_idxs) {
        float v = momentum * prev_v[used_idx] + learning_rate * grads[used_idx];
        (*(param->getWeightsPtr()))[used_idx] += v;
        prev_v[used_idx] = v;
    }
}

Adam::Adam(unsigned int length,
           float learning_rate,
           float beta_1,
           float beta_2,
           float eps,
           bool use_amsgrad) {

    first_moment_vector = std::vector<float>(length, 0.0);
    second_moment_vector = std::vector<float>(length, 0.0);
    this->learning_rate = learning_rate;
    this->beta_1 = beta_1;
    this->beta_2 = beta_2;
    this->eps = eps;
    this->iteration_count = 0;
    this->use_amsgrad = use_amsgrad;
    if (this->use_amsgrad)
        this->second_moment_max_vector = std::vector<float>(length, 0.0);
}


void Adam::adjustWeights(std::vector<float> &grads,
                         std::set<unsigned int> &used_idxs,
                         boost::shared_ptr<Param> param) {

    // Adam use one base iterator
    iteration_count += 1;
    for (auto &used_idx: used_idxs) {
        // Update biased first moment estimate
        // m_t = beta_1 * m_{t-1} + ( 1 - beta_1 ) * g_t
        float m_t = first_moment_vector[used_idx] * beta_1 + (1.0 - beta_1) * grads[used_idx];
        first_moment_vector[used_idx] = m_t;

        // Update biased second raw moment estimate
        // v_t = beta_2 * v_{t-1} + ( 1 - beta_2 ) * g_t^2
        float v_t = beta_2 * second_moment_vector[used_idx]
                     + (1.0 - beta_2) * std::pow(grads[used_idx], 2);
        second_moment_vector[used_idx] = v_t;

        if(this->use_amsgrad){
            // Referrence ON THE CONVERGENCE PROOF OF AMSGRAD https://arxiv.org/pdf/1904.03590.pdf
            // select max between v_t and all v_t-1 before it
            v_t = std::max(v_t, second_moment_max_vector[used_idx]);
            // keep track of v_hat which is the max of v_t some far
            second_moment_max_vector[used_idx] = v_t;
        }

        // Compute bias-corrected first moment estimate
        // hat_m_t = m_t / ( 1 - beta_1 ^ t)
        float m_hat = m_t / (1.0 - std::pow(beta_1, iteration_count));
        // Compute bias-corrected second raw moment estimate
        // hat_v_t = v_t / ( 1 - beta_2 ^ t)
        float v_hat = v_t / (1.0 - std::pow(beta_2, iteration_count));

        // Update parameters
        // theta_t = theta_{t-1} - alpha * m_hat / ( sqrt(v_hat) + eps)
        (*(param->getWeightsPtr()))[used_idx] += learning_rate * m_hat / (std::sqrt(v_hat) + eps);
    }
}

void AdamW::adjustWeights(std::vector<float> &grads,
                          std::set<unsigned int> &used_idxs,
                          boost::shared_ptr<Param> param) {

    std::vector<unsigned int> bias_index;
    param->getBiasIndexes(bias_index);
    // Adam use one base iterator
    iteration_count += 1;
    for (auto &used_idx: used_idxs) {
        // Update biased first moment estimate
        // m_t = beta_1 * m_{t-1} + ( 1 - beta_1 ) * g_t
        float m_t = first_moment_vector[used_idx] * beta_1 + (1.0 - beta_1) * grads[used_idx];
        first_moment_vector[used_idx] = m_t;

        // Update biased second raw moment estimate
        // v_t = beta_2 * v_{t-1} + ( 1 - beta_2 ) * g_t^2
        float v_t = beta_2 * second_moment_vector[used_idx]
                     + (1.0 - beta_2) * std::pow(grads[used_idx], 2);
        second_moment_vector[used_idx] = v_t;

        if(this->use_amsgrad){
            // Referrence ON THE CONVERGENCE PROOF OF AMSGRAD https://arxiv.org/pdf/1904.03590.pdf
            // select max between v_t and all v_t-1 before it
            v_t = std::max(v_t, second_moment_max_vector[used_idx]);
            // keep track of v_hat which is the max of v_t some far
            second_moment_max_vector[used_idx] = v_t;
        }

        // Compute bias-corrected first moment estimate
        // hat_m_t = m_t / ( 1 - beta_1 ^ t)
        float m_hat = m_t / (1.0 - std::pow(beta_1, iteration_count));
        // Compute bias-corrected second raw moment estimate
        // hat_v_t = v_t / ( 1 - beta_2 ^ t)
        float v_hat = v_t / (1.0 - std::pow(beta_2, iteration_count));

        // Update parameters
        // theta_t = theta_{t-1} - alpha * m_hat / ( sqrt(v_hat) + eps)
        if(std::find(bias_index.begin(), bias_index.end(),used_idx) != bias_index.end())
            (*(param->getWeightsPtr()))[used_idx] =  (1.0 - w) *  (*(param->getWeightsPtr()))[used_idx] + learning_rate * m_hat / (std::sqrt(v_hat) + eps);
        else
            (*(param->getWeightsPtr()))[used_idx] =  (*(param->getWeightsPtr()))[used_idx] + learning_rate * m_hat / (std::sqrt(v_hat) + eps);
    }
}

AdaDelta::AdaDelta(unsigned int length,
                   float learning_rate,
                   float decay_rate,
                   float eps) {

    mean_squared_delta_x = std::vector<float>(length, 0.0);
    mean_squared_gradients = std::vector<float>(length, 0.0);
    this->learning_rate = learning_rate;
    this->decay_rate = decay_rate;
    this->eps = eps;
    this->iteration_count = 0;
}


void AdaDelta::adjustWeights(std::vector<float> &grads,
                             std::set<unsigned int> &used_idxs,
                             boost::shared_ptr<Param> param) {
    
    // TODO: MAKE SURE THIS WORKS
    iteration_count += 1;
    for (auto &used_idx: used_idxs) {
        // Accumulate Gradient
        // E[g^2]_t = decay_rate * E[g^2]_{t-1} + ( 1 - decay_rate ) * grads_t^2
        mean_squared_gradients[used_idx] =
                decay_rate * mean_squared_gradients[used_idx] + (1 - decay_rate) * grads[used_idx] * grads[used_idx];
        // Compute Update
        // RMS[Delta_x]_{t-1} = sqrt(E[Delta_x^2]_{t-1} + eps)
        // NOTE mean_squared_delta_x has not update yet
        float rms_delta_x = sqrt(mean_squared_delta_x[used_idx] + eps);
        // RMS[g]_t = sqrt(E[g^2]_t + eps)
        float rms_gradients = sqrt(mean_squared_gradients[used_idx] + eps);
        // -RMS[Delta_x]_{t-1}/ RMS[g]_t * g_t
        float detla_x = -rms_delta_x / rms_gradients * grads[used_idx];
        // Accumulate Updates
        // E[Delta_x^2]_t  = decay_rate * E[Delta_x^2]_{t-1} + ( 1 - decay_rate ) * Delta_x_t^2
        mean_squared_delta_x[used_idx] =
                decay_rate * mean_squared_delta_x[used_idx] + (1 - decay_rate) * detla_x * detla_x;

        // Update weights
        // NOTE: we are doing gradient ascent
        (*(param->getWeightsPtr()))[used_idx] -= detla_x;
    }
}


void AdaBelief::adjustWeights(std::vector<float> &grads,
                         std::set<unsigned int> &used_idxs,
                         boost::shared_ptr<Param> param) {
    // Referrence: https://juntang-zhuang.github.io/adabelief/
    // Referrence: https://arxiv.org/pdf/2010.07468.pdf
    // Adam use one base iterator
    iteration_count += 1;
    for (auto &used_idx: used_idxs) {
        // Update biased first moment estimate
        float g_t = grads[used_idx];
        // m_t = beta_1 * m_{t-1} + ( 1 - beta_1 ) * g_t
        float m_t = first_moment_vector[used_idx] * beta_1 + (1.0 - beta_1) * g_t;
        first_moment_vector[used_idx] = m_t;

        // Update biased second raw moment estimate
        // v_t = beta_2 * v_{t-1} + ( 1 - beta_2 ) * (g_t - m_t)^2
        float s_t = beta_2 * second_moment_vector[used_idx]
                     + (1.0 - beta_2) * std::pow(g_t - m_t, 2);
        second_moment_vector[used_idx] = s_t;

        if(this->use_amsgrad){
            // Referrence ON THE CONVERGENCE PROOF OF AMSGRAD https://arxiv.org/pdf/1904.03590.pdf
            // select max between v_t and all v_t-1 before it
            s_t = std::max(s_t, second_moment_max_vector[used_idx]);
            // keep track of v_hat which is the max of v_t some far
            second_moment_max_vector[used_idx] = s_t;
        }
        
        // Compute bias-corrected first moment estimate
        // hat_m_t = m_t / ( 1 - beta_1 ^ t)
        float m_hat = m_t / (1.0 - std::pow(beta_1, iteration_count));
        // Compute bias-corrected second raw moment estimate
        // hat_s_t = s_t / ( 1 - beta_2 ^ t)
        float s_hat = s_t / (1.0 - std::pow(beta_2, iteration_count));

        // Update parameters
        // theta_t = theta_{t-1} - alpha * m_hat / ( sqrt(v_hat) + eps)
        (*(param->getWeightsPtr()))[used_idx] += learning_rate * m_t / (std::sqrt(s_hat) + eps);
    }
}