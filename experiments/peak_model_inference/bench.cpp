#include "qalgorithms_peak_inference.h"

#include <chrono>
#include <cmath>
#include <iostream>

using namespace qAlgorithms;

namespace
{
    volatile double resultSink = 0;
}

int main()
{
    constexpr size_t length = 49, center = 24;
    for (size_t traces : {1u, 4u, 16u})
    {
        std::vector<std::vector<float>> raw(traces, std::vector<float>(length));
        std::vector<const float *> pointers;
        for (size_t j = 0; j < traces; ++j)
        {
            for (size_t i = 0; i < length; ++i)
            {
                const double u = .1 * (static_cast<double>(i) - center);
                raw[j][i] = static_cast<float>(std::exp(7 + .02 * j + .6 * u
                    - (2 + .02 * j) * std::cosh(u) - .2 * std::sinh(u) + .01 * std::sin(i * 1.4 + j)));
            }
            pointers.push_back(raw[j].data());
        }
        for (size_t trial = 0; trial < 9; ++trial)
        {
            const auto start = std::chrono::steady_clock::now();
            for (size_t pass = 0; pass < 5; ++pass)
            {
                const auto inference = peak_model::inferHyperbolicWindow(pointers, length, center, center, .1);
                if (inference.status != peak_model::InferenceStatus::available)
                {
                    return 1;
                }
                resultSink = inference.covariance.back();
            }
            const auto duration = std::chrono::steady_clock::now() - start;
            const double nanoseconds = std::chrono::duration<double, std::nano>(duration).count() / 5;
            std::cout << traces << ',' << trial << ',' << nanoseconds << '\n';
        }
    }
    return 0;
}
