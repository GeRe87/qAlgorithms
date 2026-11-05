#include "../src/qalgorithms_qbin.cpp" // NOLINT
#include "common_test_utils.hpp"
#include <limits>

using namespace qAlgorithms;

static void test_create_eic_counts_measured_points()
{
    // Deliberately unsorted scans with a gap: cumulative DF counts measured
    // centroids, regardless of scan number or the order they were inserted.
    const uint16_t input_scans[] = {5, 2, 7, 4, 6};
    const uint32_t sorted_scans[] = {2, 4, 5, 6, 7};
    const std::vector<float> retention_times = {10, 10.5F, 11, 11.5F, 12, 12.5F, 13, 13.5F};
    std::vector<CentroidPeak> centroids(5);
    Bin bin;
    for (size_t i = 0; i < centroids.size(); ++i)
    {
        centroids[i].number_MS1 = input_scans[i];
        centroids[i].ID = 100 + input_scans[i];
        centroids[i].mz = 100;
        centroids[i].area = 20;
        centroids[i].height = 10;
        centroids[i].DQSC = 0.9F;
        bin.pointsInBin.push_back(&centroids[i]);
    }

    const EIC eic = bin.createEIC(&retention_times);
    assert(eic.df.size() == 5, "The EIC must retain all five measured centroids\n", NULL);
    for (size_t i = 0; i < 5; ++i)
    {
        assert(eic.df[i] == i + 1,
               "Cumulative measured-point count at index %zu is %u instead of %zu\n",
               i, eic.df[i], i + 1);
        assert(eic.scanNumbers[i] == sorted_scans[i],
               "EIC must be ordered by scan number\n", NULL);
        assert(eic.cenID[i] == 100 + sorted_scans[i],
               "Centroid identity changed while ordering the EIC\n", NULL);
        assert(eic.RT[i] == retention_times[sorted_scans[i]],
               "EIC retention time must match its source scan\n", NULL);
    }
    assert(sumOfCumulative(eic.df.data(), 0, 5) == 5,
           "The complete EIC must contain five measured points\n", NULL);
    assert(sumOfCumulative(eic.df.data(), 1, 3) == 3,
           "A three-centroid subrange must contain three measured points\n", NULL);
}

static void test_create_eic_area_domain()
{
    const std::vector<float> times = {0, 1, 2, 3, 4};
    std::vector<CentroidPeak> points(5);
    Bin bin;
    for (size_t i = 0; i < points.size(); ++i)
    {
        points[i].number_MS1 = uint16_t(i);
        points[i].mz = 100;
        points[i].area = 0.125F;
        bin.pointsInBin.push_back(&points[i]);
    }
    const EIC eic = bin.createEIC(&times);
    verify(eic.ints_area.size() == points.size());
    for (float area : eic.ints_area)
        verify(area == 0.125F);
    for (float invalidArea : {0.0F, -1.0F, std::numeric_limits<float>::infinity(),
                              std::numeric_limits<float>::quiet_NaN()})
    {
        points[0].area = invalidArea;
        bool rejected = false;
        try
        {
            (void)bin.createEIC(&times);
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
        }
        verify(rejected);
    }
}

int main(void)
{
    test_create_eic_counts_measured_points();
    test_create_eic_area_domain();
    return 0;
}
