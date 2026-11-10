#include "qalgorithms_peak_search.h"

#include <chrono>
#include <cmath>
#include <iostream>

using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

namespace
{
    volatile double resultSink = 0;
}

int main()
{
    constexpr size_t length = 49;
    HyperbolicSearchOptions options;
    options.selection = WindowSelection::prespecified;
    for (size_t traces : {1u, 4u, 16u})
    {
        std::vector<std::vector<float>> raw(traces, std::vector<float>(length));
        std::vector<const float *> pointers;
        for (size_t j = 0; j < traces; ++j)
        {
            for (size_t i = 0; i < length; ++i)
            {
                const double u = .1 * (static_cast<double>(i) - 24);
                raw[j][i] = static_cast<float>(std::exp(7 + .02 * j + .5 * u
                    - 2 * std::cosh(u) - .2 * std::sinh(u) + .01 * std::sin(i * 1.4 + j)));
            }
            pointers.push_back(raw[j].data());
        }
        const auto plan = makeUniformHyperbolicSearchPlan(length, .1, 6, 12, 4, options);
        for (size_t trial = 0; trial < 9; ++trial)
        {
            const auto start = std::chrono::steady_clock::now();
            const auto result = plan.evaluate(pointers);
            const auto elapsed = std::chrono::steady_clock::now() - start;
            resultSink = static_cast<double>(result.candidateIndices.size());
            std::cout << traces << ',' << trial << ',' << plan.windows().size() << ','
                      << std::chrono::duration<double, std::micro>(elapsed).count() << '\n';
        }
    }
    return 0;
}
