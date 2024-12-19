//
// Created by feiw on 15/04/18.
//

#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include "Param.h"

class Solver {
public:
    virtual void adjustWeights(std::vector<float> &grads,
                               std::set<unsigned int> &used_idxs,
                               boost::shared_ptr<Param> param) = 0;

    void setLearningRate(const float lr) { this->learning_rate = lr; };

    virtual ~Solver() = default;
protected:
    float learning_rate;
};

class Sgd: public Solver {
public:
    Sgd(float learning_rate);
    void adjustWeights(std::vector<float> &grads,
                       std::set<unsigned int> &used_idxs,
                       boost::shared_ptr<Param> param) override;
};

class Momentum : public Solver {
public:
    Momentum(unsigned int length, float learning_rate, float momentum);

    void adjustWeights(std::vector<float> &grads,
                       std::set<unsigned int> &used_idxs,
                       boost::shared_ptr<Param> param) override;

private:
    std::vector<float> prev_v;
    float momentum;

};

class Adam : public Solver {
public:
    Adam(unsigned int length,
         float learning_rate,
         float beta_1,
         float beta_2,
         float eps,
         bool use_amsgrad);

    void adjustWeights(std::vector<float> &grads,
                       std::set<unsigned int> &used_idxs,
                       boost::shared_ptr<Param> param) override;

protected:
    float beta_1;
    float beta_2;
    float eps;
    int iteration_count;
    bool use_amsgrad;
    std::vector<float> first_moment_vector;
    std::vector<float> second_moment_vector;
    std::vector<float> second_moment_max_vector;
};

// simple version of Adam W
// where eta is always 1.0
class AdamW : public Adam {
public:
    AdamW(unsigned int length,
         float learning_rate,
         float beta_1,
         float beta_2,
         float eps,
         bool use_amsgrad,
         float w
         ) : Adam(length, learning_rate, beta_1, beta_2, eps, use_amsgrad) {
        this->w = w;
    };

    void adjustWeights(std::vector<float> &grads,
                       std::set<unsigned int> &used_idxs,
                       boost::shared_ptr<Param> param) override;
    
private:
    float w;
};

// Ada Belief 
// https://arxiv.org/abs/2010.07468
class AdaBelief : public Adam {
public:
    AdaBelief(unsigned int length,
         float learning_rate,
         float beta_1,
         float beta_2,
         float eps,
         bool use_amsgrad) : Adam(length, learning_rate, beta_1, beta_2, eps, use_amsgrad) {};

    void adjustWeights(std::vector<float> &grads,
                       std::set<unsigned int> &used_idxs,
                       boost::shared_ptr<Param> param) override;
};

class AdaDelta : public Solver {
public:
    AdaDelta(unsigned int length,
             float learning_rate,
             float decay_rate,
             float eps);

    void adjustWeights(std::vector<float> &grads,
                       std::set<unsigned int> &used_idxs,
                       boost::shared_ptr<Param> param) override;

private:
    float decay_rate;
    float eps;
    int iteration_count;
    std::vector<float> mean_squared_delta_x;
    std::vector<float> mean_squared_gradients;
};


#endif //SOLVER_H
