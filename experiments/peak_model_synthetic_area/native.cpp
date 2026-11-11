#include "qalgorithms_peak_inference.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace qAlgorithms;
using namespace qAlgorithms::peak_model;

int main()
{
    try
    {
        using Legacy = Model<PeakModelKind::PiecewiseQuadratic>;
        const double missing = std::numeric_limits<double>::quiet_NaN();
        std::cout << std::setprecision(17);
        size_t scale;
        double step;
        while (std::cin >> scale >> step)
        {
            if ((scale != 4 && scale != 8) || !(step > 0 && std::isfinite(step)))
            {
                throw std::runtime_error("Invalid simulation geometry");
            }
            const size_t n = 2 * scale + 1;
            std::vector<float> raw(n), logs(n), predicted(n);
            for (size_t i = 0; i < n; ++i)
            {
                if (!(std::cin >> raw[i]) || !(raw[i] > 0 && std::isfinite(raw[i])))
                {
                    throw std::runtime_error("Invalid simulation intensity");
                }
                logs[i] = static_cast<float>(std::log(static_cast<double>(raw[i])));
            }
            std::vector<RegCoeffs> fits;
            Legacy::runningFit(logs.data(), n, &fits);
            auto legacy = fits.back();
            if (legacy.x0 != scale || legacy.scale != scale)
            {
                throw std::runtime_error("Missing central fit");
            }
            double rss = 0;
            for (size_t i = 0; i < n; ++i)
            {
                const double residual = logs[i] - Legacy::logValue(&legacy, static_cast<double>(i) - scale);
                rss += residual * residual;
            }
            Legacy::refitAmplitude(raw.data(), {0, static_cast<int32_t>(n)}, predicted.data(), &legacy);
            double legacyArea = missing, legacySE = missing;
            if (legacy.b2 < 0 && legacy.b3 < 0)
            {
                legacyArea = Legacy::integral(&legacy, step, rss / (n - 4), &legacySE);
            }
            std::cout << legacy.b0 << ' ' << legacy.b1 / step << ' ' << legacy.b2 / (step * step) << ' '
                      << legacy.b3 / (step * step) << ' ' << legacyArea << ' ' << legacySE;

            HyperbolicModel::Geometry geometry(scale, step);
            const auto inference = inferHyperbolicWindow({raw.data()}, n, scale, scale, geometry);
            HyperbolicModel::Fit selected;
            HyperbolicModel::runningFit({logs.data()}, n, geometry, [&](const auto &fit)
            {
                if (fit.center == scale && fit.scale == scale)
                {
                    selected = fit;
                }
            });
            const auto properties = evaluateHyperbolicWindow(selected, 0, geometry, raw.data(), n);
            const double hyperSE = inference.status == InferenceStatus::available
                ? inference.properties.front().areaUncertainty.value : missing;
            const auto &shape = selected.locals.front();
            std::cout << ' ' << properties.correctedIntercept << ' ' << selected.sharedSlope << ' ' << shape.curvature
                      << ' ' << shape.skew << ' ' << properties.areaU.value << ' ' << hyperSE << '\n';
        }
        if (!std::cin.eof() || !std::cout)
        {
            throw std::runtime_error("Failed to complete native fit comparison");
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
