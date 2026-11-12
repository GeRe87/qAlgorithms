#include "qalgorithms_superspectrum_tracks.h"

#include <sstream>
#include <stdexcept>

namespace
{
    void require(bool condition)
    {
        if (!condition)
            throw std::runtime_error("Superspectrum track association failed");
    }

    qAlgorithms::SuperspectrumObservation peak(double mz, uint32_t scan, double height, size_t distance)
    {
        return {mz, mz - .01, mz + .01, .001, {{scan, 1, height, distance}}};
    }
}

int main()
{
    std::ostringstream output;
    qAlgorithms::SuperspectrumTracks tracks(output, 5);
    const auto initial = tracks.append({peak(100, 4, 20, 2), peak(100.006, 4, 60, 2)}, 0);
    require(initial.size() == 2 && initial[0] != initial[1]);
    tracks.append({}, 1);
    tracks.append({}, 2);
    tracks.append({}, 3);
    // The last possible overlapping window can reconnect to the closest track.
    const auto resumed = tracks.append({peak(100.005, 4, 61, 2)}, 4);
    require(resumed.size() == 1 && resumed[0] == initial[1]);
    // A mass match alone cannot connect a window without a shared scan.
    const auto distinct = tracks.append({peak(100.005, 5, 80, 2)}, 5);
    require(distinct.size() == 1 && distinct[0] != initial[0] && distinct[0] != initial[1]);
    tracks.finish();
    std::istringstream rows(output.str());
    std::string line;
    size_t count = 0;
    bool singleton = false, earlierHeight = false;
    while (std::getline(rows, line))
    {
        ++count;
        singleton |= line == "1,1,100,,4,1,20,0";
        if (line.rfind("2,2,", 0) == 0)
            earlierHeight = line.substr(line.find(",4,1,")) == ",4,1,60,0";
    }
    require(count == 4 && singleton && earlierHeight);
}
