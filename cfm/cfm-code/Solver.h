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
    virtual void adjustWeights(std::vector<double> &grads,
                               std::set<unsigned int> &used_idxs,
                               boost::shared_ptr<Param> param);

    virtual void adjustWeights(std::vector<double> &grads, std::vector<double> &weights,
                               std::set<unsigned int> &used_idxs)=0;

    void setLearningRate(const double lr) { this->learning_rate = lr; };

    virtual ~Solver() = default;
protected:
    double learning_rate;
};

class Sgd: public Solver {
public:
    Sgd(double learning_rate);
    void adjustWeights(std::vector<double> &grads,
                       std::vector<double> &weights,
                       std::set<unsigned int> &used_idxs) override;
};

class Momentum : public Solver {
public:
    Momentum(unsigned int length, double learning_rate, double momentum);

    void adjustWeights(std::vector<double> &grads,
                       std::vector<double> &weights,
                       std::set<unsigned int> &used_idxs) override;

private:
    std::vector<double> prev_v;
    double momentum;

};

class Adam : public Solver {
public:
    Adam(unsigned int length,
         double learning_rate,
         double beta_1,
         double beta_2,
         double eps);

    void adjustWeights(std::vector<double> &grads,
                       std::vector<double> &weights,
                       std::set<unsigned int> &used_idxs) override;

protected:
    double beta_1;
    double beta_2;
    double eps;
    int iteration_count;
    std::vector<double> first_moment_vector;
    std::vector<double> second_moment_vector;

};

// simple version of Adam W
// where eta is always 1.0
class AdamW : public Adam {
public:
    AdamW(unsigned int length,
         double learning_rate,
         double beta_1,
         double beta_2,
         double eps,
          double w) : Adam(length, learning_rate, beta_1, beta_2, eps) {
        this->w = w;
    };

    void adjustWeights(std::vector<double> &grads,
                       std::vector<double> &weights,
                       std::set<unsigned int> &used_idxs) override;

    void add_bais_index(int index) {
        bias_indices.insert(index);
    }

    void add_bais_indices(std::vector<int> & indices) {
        for(const auto & index : indices)
            bias_indices.insert(index);
    }
    
private:
    double w;
    std::set<int> bias_indices;
};

class Adadelta : public Solver {
public:
    Adadelta(unsigned int length,
             double learning_rate,
             double decay_rate,
             double eps);

    void adjustWeights(std::vector<double> &grads,
                       std::vector<double> &weights,
                       std::set<unsigned int> &used_idxs) override;

private:
    double decay_rate;
    double eps;
    int iteration_count;
    std::vector<double> mean_squared_delta_x;
    std::vector<double> mean_squared_gradients;
};


class AMSgrad : public Solver {
public:
    AMSgrad(unsigned int length,
            double learning_rate,
            double beta_1,
            double beta_2,
            double eps);

    void adjustWeights(std::vector<double> &grads,
                       std::vector<double> &weights,
                       std::set<unsigned int> &used_idxs) override;

private:
    double beta_1;
    double beta_2;
    double eps;
    int iteration_count;
    std::vector<double> first_moment_vector;
    std::vector<double> second_moment_vector;
    std::vector<double> second_moment_max_vector;
};


#endif //SOLVER_H
