//
// Created by feiw on 15/04/18.
//

#include "Solver.h"

void Solver::adjustWeights(std::vector<double> &grads,
                           std::set<unsigned int> &used_idxs, boost::shared_ptr<Param> param) {
    adjustWeights(grads, *(param->getWeightsPtr()), used_idxs);
}

Sgd::Sgd(double learning_rate){
    this->learning_rate = learning_rate;
}

void Sgd::adjustWeights(std::vector<double> &grads,
                             std::vector<double> &weights,
                             std::set<unsigned int> &used_idxs) {
    for (auto &used_idx: used_idxs)
        weights[used_idx] += learning_rate * grads[used_idx];
}


Momentum::Momentum(unsigned int length, double learning_rate, double momentum) {
    prev_v = std::vector<double>(length, 0.0);
    this->learning_rate = learning_rate;
    this->momentum = momentum;
}

void Momentum::adjustWeights(std::vector<double> &grads,
                             std::vector<double> &weights,
                             std::set<unsigned int> &used_idxs) {

    for (auto &used_idx: used_idxs) {
        double v = momentum * prev_v[used_idx] + learning_rate * grads[used_idx];
        weights[used_idx] += v;
        prev_v[used_idx] = v;
    }
}

Aadm::Aadm(unsigned int length,
           double learning_rate,
           double beta_1,
           double beta_2,
           double eps) {

    first_moment_vector = std::vector<double>(length, 0.0);
    second_moment_vector = std::vector<double>(length, 0.0);
    this->learning_rate = learning_rate;
    this->beta_1 = beta_1;
    this->beta_2 = beta_2;
    this->eps = eps;
    this->iteration_count = 0;
}


void Aadm::adjustWeights(std::vector<double> &grads,
                         std::vector<double> &weights,
                         std::set<unsigned int> &used_idxs) {

    // Adam use one base iterator
    iteration_count += 1;
    for (auto &used_idx: used_idxs) {
        // Update biased first moment estimate
        // m_t = beta_1 * m_{t-1} + ( 1 - beta_1 ) * g_t
        double m_t = first_moment_vector[used_idx] * beta_1 + (1.0 - beta_1) * grads[used_idx];
        first_moment_vector[used_idx] = m_t;

        // Update biased second raw moment estimate
        // v_t = beta_2 * v_{t-1} + ( 1 - beta_2 ) * g_t^2
        double v_t = beta_2 * second_moment_vector[used_idx]
                     + (1.0 - beta_2) * std::pow(grads[used_idx], 2);
        second_moment_vector[used_idx] = v_t;

        // Compute bias-corrected first moment estimate
        // hat_m_t = m_t / ( 1 - beta_1 ^ t)
        double m_hat = m_t / (1.0 - std::pow(beta_1, iteration_count));
        // Compute bias-corrected second raw moment estimate
        // hat_v_t = v_t / ( 1 - beta_2 ^ t)
        double v_hat = v_t / (1.0 - std::pow(beta_2, iteration_count));

        // Update parameters
        // theta_t = theta_{t-1} - alpha * m_hat / ( sqrt(v_hat) + eps)
        weights[used_idx] += learning_rate * m_hat / (std::sqrt(v_hat) + eps);
    }
}

// AMSgrad sucks
AMSgrad::AMSgrad(unsigned int length,
                 double learning_rate,
                 double beta_1,
                 double beta_2,
                 double eps) {

    first_moment_vector = std::vector<double>(length, 0.0);
    second_moment_vector = std::vector<double>(length, 0.0);
    second_moment_max_vector = std::vector<double>(length, 0.0);
    this->learning_rate = learning_rate;
    this->beta_1 = beta_1;
    this->beta_2 = beta_2;
    this->eps = eps;
    this->iteration_count = 0;
}


void AMSgrad::adjustWeights(std::vector<double> &grads,
                            std::vector<double> &weights,
                            std::set<unsigned int> &used_idxs) {
    // Adam use one base iterator
    iteration_count += 1;
    for (auto &used_idx: used_idxs) {
        // Update biased first moment estimate
        // m_t = beta_1 * m_{t-1} + ( 1 - beta_1 ) * g_t
        double m_t = first_moment_vector[used_idx] * beta_1 + (1.0 - beta_1) * grads[used_idx];
        first_moment_vector[used_idx] = m_t;

        // Update biased second raw moment estimate
        // v_t = beta_2 * v_{t-1} + ( 1 - beta_2 ) * g_t^2
        double v_t = beta_2 * second_moment_vector[used_idx]
                     + (1.0 - beta_2) * std::pow(grads[used_idx], 2);
        second_moment_vector[used_idx] = v_t;

        double v_hat = std::max(v_t, second_moment_max_vector[used_idx]);
        // keep track of v_hat which is the max of v_t some far
        second_moment_max_vector[used_idx] = v_hat;

        // Compute bias-corrected first moment estimate
        // hat_m_t = m_t / ( 1 - beta_1 ^ t)
        double m_hat = m_t / (1.0 - std::pow(beta_1, iteration_count));
        // Compute bias-corrected second raw moment estimate
        // hat_v_t = v_t / ( 1 - beta_2 ^ t)
        v_hat = v_t / (1.0 - std::pow(beta_2, iteration_count));

        // Update parameters
        // theta_t = theta_{t-1} - alpha * m_hat / ( sqrt(v_hat) + eps)
        weights[used_idx] += learning_rate * m_hat / (sqrt(v_hat) + eps);
    }
}

Adadelta::Adadelta(unsigned int length,
                   double learning_rate,
                   double decay_rate,
                   double eps) {

    mean_squared_delta_x = std::vector<double>(length, 0.0);
    mean_squared_gradients = std::vector<double>(length, 0.0);
    this->learning_rate = learning_rate;
    this->decay_rate = decay_rate;
    this->eps = eps;
    this->iteration_count = 0;
}


void Adadelta::adjustWeights(std::vector<double> &grads,
                             std::vector<double> &weights,
                             std::set<unsigned int> &used_idxs) {

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
        double rms_delta_x = sqrt(mean_squared_delta_x[used_idx] + eps);
        // RMS[g]_t = sqrt(E[g^2]_t + eps)
        double rms_gradients = sqrt(mean_squared_gradients[used_idx] + eps);
        // -RMS[Delta_x]_{t-1}/ RMS[g]_t * g_t
        double detla_x = -rms_delta_x / rms_gradients * grads[used_idx];
        // Accumulate Updates
        // E[Delta_x^2]_t  = decay_rate * E[Delta_x^2]_{t-1} + ( 1 - decay_rate ) * Delta_x_t^2
        mean_squared_delta_x[used_idx] =
                decay_rate * mean_squared_delta_x[used_idx] + (1 - decay_rate) * detla_x * detla_x;

        // Update weights
        // NOTE: we are doing gradient ascent
        weights[used_idx] -= detla_x;
    }
}