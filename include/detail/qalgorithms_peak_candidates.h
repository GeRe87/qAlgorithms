#ifndef QALGORITHMS_PEAK_CANDIDATES_H
#define QALGORITHMS_PEAK_CANDIDATES_H

namespace qAlgorithms
{
    namespace detail
    {
        // Each model owns its candidate type and validity rules. The
        // collector sees every decision in input order, before selection.
        template<typename Candidate, typename Inputs, typename Validate, typename Collect>
        void findValidRegressions(const Inputs &inputs, Validate validate, Collect collect)
        {
            for (const auto &input : inputs)
            {
                Candidate candidate;
                const bool valid = validate(input, &candidate);
                collect(candidate, valid);
            }
        }
    }
}

#endif
