#include "qalgorithms_read_file.h"
#include <cstddef>
#include <cstdint>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#pragma GCC diagnostic ignored "-Wfloat-conversion"

#include "../src/qalgorithms_read_file.cpp" // NOLINT
#include "common_test_utils.hpp"

using namespace qAlgorithms;

static void test_spec_decode_nocomp_double()
{
    const double arr_D[10] = {0.135761618738162, 2.75537138578063, 1.72069261613833, 0.0218284366498446, 1.65262662060337, 0.637948708991408, 0.0385504743656586, 0.275107929479815, 1.71114153324745, 0.152519432142112};
    // reference string encoded using python script (see end of file)
    const char encoded[] = "jiZEAKNgwT+WGDInAAsGQGOGDPv0h/s/KLQiljFalj/72H+cKHH6P7SXNGkTauQ//qRKReO8oz/uG/9JXpvRP2MAwvHVYPs/ugWHusGFwz8=";
    const size_t encode_len = 108; // 109 above includes null terminator

    std::vector<float> result;
    assert(decodeSpectrum(encoded, sizeof(encoded) - 1, &result, 10, false, true) == 0,
           "Uncompressed double decoding failed");

    assert(result.size() == 10, "Incorrect number of Elements decoded\n", NULL);
    for (size_t i = 0; i < 10; i++)
    {
        assert(flt_equal(arr_D[i], result[i], FLT_EPSILON),
               "Element %d was decoded incorrectly (difference: %f)",
               i, abs(arr_D[i] - result[i]));
    }

    // also check if encoding works
    std::vector<char> encode = encode_base64_dbl(arr_D, 10);
    assert(encode.size() == encode_len, "Incorrect length of encoding", NULL);
    for (size_t i = 0; i < encode_len; i++)
    {
        assert(encode[i] == encoded[i], "Non-Matching encoding result", NULL);
    }
}

static void test_spec_decode_iscomp_double()
{
    const double arr_D[10] = {0.135761618738162, 2.75537138578063, 1.72069261613833, 0.0218284366498446, 1.65262662060337, 0.637948708991408, 0.0385504743656586, 0.275107929479815, 1.71114153324745, 0.152519432142112};

    const uint8_t *arr_char = (uint8_t *)arr_D;
    const size_t arr_char_len = 10 * sizeof(double);
    std::vector<char> intermediate;

    compress_and_encode(arr_char, arr_char_len, &intermediate);

    std::vector<float> result;
    assert(decodeSpectrum(intermediate.data(), intermediate.size(), &result, 10, true, true) == 0,
           "Compressed double decoding failed");

    assert(result.size() == 10, "Incorrect number of Elements decoded\n", NULL);
    for (size_t i = 0; i < 10; i++)
    {
        assert(flt_equal(arr_D[i], result[i], FLT_EPSILON),
               "Element %d was decoded incorrectly (difference: %f)",
               i, abs(arr_D[i] - result[i]));
    }
}

static void test_spec_decode_iscomp_float()
{
    const float arr_F[10] = {0.135761618738162, 2.75537138578063, 1.72069261613833, 0.0218284366498446, 1.65262662060337, 0.637948708991408, 0.0385504743656586, 0.275107929479815, 1.71114153324745, 0.152519432142112};

    const uint8_t *arr_char = (uint8_t *)arr_F;
    const size_t arr_char_len = 10 * sizeof(float);
    std::vector<char> intermediate;

    compress_and_encode(arr_char, arr_char_len, &intermediate);

    std::vector<float> result;
    assert(decodeSpectrum(intermediate.data(), intermediate.size(), &result, 10, true, false) == 0,
           "Compressed float decoding failed");

    assert(result.size() == 10, "Incorrect number of Elements decoded\n", NULL);
    for (size_t i = 0; i < 10; i++)
    {
        assert(flt_equal(arr_F[i], result[i], FLT_EPSILON),
               "Element %d was decoded incorrectly (difference: %f)",
               i, abs(arr_F[i] - result[i]));
    }
}

template <class T>
static std::string encoded_values(const std::vector<T> &values, const bool compressed)
{
    const size_t byte_count = values.size() * sizeof(T);
    std::vector<char> encoded;
    if (compressed)
        compress_and_encode(reinterpret_cast<const uint8_t *>(values.data()), byte_count, &encoded);
    else
    {
        encoded.resize(simdutf::base64_length_from_binary(byte_count));
        const size_t written = simdutf::binary_to_base64(reinterpret_cast<const char *>(values.data()), byte_count, encoded.data());
        assert(written == encoded.size(), "Could not encode test values");
    }
    return {encoded.begin(), encoded.end()};
}

template <class T>
static std::string binary_array(const std::vector<T> &values, const bool compressed, const bool mz)
{
    const std::string encoded = encoded_values(values, compressed);
    return "<binaryDataArray encodedLength='" + std::to_string(encoded.size()) + "'>"
           "<cvParam accession='" + std::string(mz ? "MS:1000514" : "MS:1000515") + "'/>"
           "<cvParam accession='" + std::string(sizeof(T) == sizeof(double) ? "MS:1000523" : "MS:1000521") + "'/>"
           "<cvParam accession='" + std::string(compressed ? "MS:1000574" : "MS:1000576") + "'/>"
           "<binary>" + encoded + "</binary></binaryDataArray>";
}

static std::string test_spectrum(const std::string &arrays, const std::string &rt = "1",
                                 const std::string &polarity = "MS:1000130", const std::string &level = "1",
                                 const std::string &unit = "UO:0000010", const std::string &length = "3")
{
    return "<spectrum id='test' defaultArrayLength='" + length + "'>"
           "<cvParam accession='MS:1000511' value='" + level + "'/>"
           "<cvParam accession='MS:1000128'/><cvParam accession='" + polarity + "'/>"
           "<scanList count='1'><scan><cvParam accession='MS:1000016' value='" + rt +
           "' unitAccession='" + unit + "'/></scan></scanList>"
           "<binaryDataArrayList count='2'>" + arrays + "</binaryDataArrayList></spectrum>";
}

static std::string test_document(const std::string &spectra, const std::string &count = "1", const bool indexed = false)
{
    const std::string body = "<mzML xmlns='http://psi.hupo.org/ms/mzml'><run id='test'>"
                             "<spectrumList count='" + count + "'>" + spectra + "</spectrumList></run></mzML>";
    return indexed ? "<indexedmzML>" + body + "</indexedmzML>" : body;
}

struct TemporaryMzml
{
    std::filesystem::path directory;
    std::filesystem::path path;

    explicit TemporaryMzml(const std::string &xml)
    {
        static size_t sequence = 0;
        const auto stamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        directory = std::filesystem::temp_directory_path() /
                    ("qalgorithms-parser-test-" + std::to_string(stamp) + "-" + std::to_string(sequence++));
        assert(std::filesystem::create_directory(directory), "Could not create fixture directory");
        path = directory / "input.mzML";
        std::ofstream output(path, std::ios::binary);
        output << xml;
        assert(output.good(), "Could not write mzML fixture");
    }

    ~TemporaryMzml()
    {
        std::error_code ignored;
        std::filesystem::remove_all(directory, ignored);
    }
};

static void test_independent_array_encodings_and_roots()
{
    const std::vector<double> mz64 = {100.25, 101.5, 102.75};
    const std::vector<float> intensity32 = {10, 20, 30};
    const std::vector<float> mz32 = {200.25F, 201.5F, 202.75F};
    const std::vector<double> intensity64 = {40, 50, 60};
    // Intensity precedes m/z, and each array has its own precision/compression.
    const std::string first = test_spectrum(binary_array(intensity32, true, false) + binary_array(mz64, false, true));
    const std::string second = test_spectrum(binary_array(mz32, true, true) + binary_array(intensity64, false, false), "2");

    for (const bool indexed : {false, true})
    {
        TemporaryMzml fixture(test_document(first + second, "2", indexed));
        XML_File file(fixture.path.c_str(), mzML);
        assert(!file.defective && file.number_spectra == 2, "Valid indexed or plain mzML was rejected");
        std::vector<float> mz;
        std::vector<float> intensity;
        assert(get_spectrum(&file, &mz, &intensity, 0) == 0, "Mixed-encoding spectrum was rejected");
        assert(mz == std::vector<float>({100.25F, 101.5F, 102.75F}) && intensity == intensity32,
               "Array types, widths or compression were mixed up");
        assert(get_spectrum(&file, &mz, &intensity, 1) == 0, "Per-spectrum encoding was not respected");
        assert(mz == mz32 && intensity == std::vector<float>({40, 50, 60}), "Second spectrum decoded incorrectly");
        file.free_linknodes();
        file.free_linknodes(); // Manual cleanup and the destructor must both be safe.
    }

    TemporaryMzml single(test_document(first));
    XML_File one_scan(single.path.c_str(), mzML);
    assert(!one_scan.defective && one_scan.polarityMode == Polarities::positive,
           "A single-spectrum mzML must not fail the polarity scan");
}

static void test_selected_retention_times()
{
    const std::string arrays = binary_array(std::vector<double>{100, 101, 102}, false, true) +
                               binary_array(std::vector<float>{10, 20, 30}, false, false);
    const std::string spectra = test_spectrum(arrays, "1") +
                                test_spectrum(arrays, "2", "MS:1000129") +
                                test_spectrum(arrays, "3", "MS:1000130", "2") +
                                test_spectrum(arrays, "0.1", "MS:1000130", "1", "UO:0000031");
    TemporaryMzml fixture(test_document(spectra, "4"));
    XML_File file(fixture.path.c_str(), mzML);
    assert(!file.defective && file.polarityMode == Polarities::mixed, "Mixed-polarity fixture was rejected");
    const auto selected = filter_spectra(&file, true, Polarities::positive, true);
    assert(selected == std::vector<uint32_t>({0, 3}), "MS2 and negative scans were not filtered");
    std::vector<float> times;
    get_spectra_RT(&file, &selected, &times);
    assert(times == std::vector<float>({1, 6}), "RTs must use the selected indices and convert minutes to seconds");
    const std::vector<uint32_t> invalid_order = {3, 0};
    get_spectra_RT(&file, &invalid_order, &times);
    assert(times.empty(), "Non-increasing selected RTs must be rejected");
    const std::vector<uint32_t> invalid_index = {99};
    get_spectra_RT(&file, &invalid_index, &times);
    assert(times.empty(), "An invalid RT index must be rejected without an assertion");
}

static void test_referenced_array_parameters()
{
    const std::vector<double> values = {100, 101, 102};
    const std::string arrays = "<binaryDataArray><referenceableParamGroupRef ref='mass'/>"
                               "<binary>" + encoded_values(values, true) + "</binary></binaryDataArray>" +
                               binary_array(std::vector<float>{10, 20, 30}, false, false);
    std::string xml = test_document(test_spectrum(arrays));
    const std::string group = "<referenceableParamGroupList count='1'><referenceableParamGroup id='mass'>"
                              "<cvParam accession='MS:1000514'/><cvParam accession='MS:1000523'/>"
                              "<cvParam accession='MS:1000574'/></referenceableParamGroup></referenceableParamGroupList>";
    xml.insert(xml.find("<run"), group);
    TemporaryMzml fixture(xml);
    XML_File file(fixture.path.c_str(), mzML);
    std::vector<float> mz;
    std::vector<float> intensity;
    assert(!file.defective && get_spectrum(&file, &mz, &intensity, 0) == 0,
           "Referenced binary-array parameters were not resolved");
    assert(mz == std::vector<float>({100, 101, 102}), "Referenced array decoded incorrectly");
}

static void test_bad_documents_fail_cleanly()
{
    const std::string arrays = binary_array(std::vector<double>{100, 101, 102}, false, true) +
                               binary_array(std::vector<float>{10, 20, 30}, false, false);
    const std::string valid_scan = test_spectrum(arrays);
    const std::vector<std::string> invalid_documents = {
        "<notMzML/>", "<mzML>", "<mzML/>", "<indexedmzML/>",
        test_document("", "0"), test_document(valid_scan, "2"),
        test_document(valid_scan, "1junk"), test_document(valid_scan, "184467440737095516160"),
        test_document(test_spectrum(arrays, "1", "missing-polarity")),
        test_document(test_spectrum(arrays, "1", "MS:1000130", "invalid-level")),
        test_document(test_spectrum(arrays, "not-a-time")),
        test_document(test_spectrum(arrays, "1", "MS:1000130", "1", "unknown-unit")),
        test_document(test_spectrum(arrays, "1", "MS:1000130", "1", "UO:0000010", "3bad"))};
    for (const std::string &xml : invalid_documents)
    {
        TemporaryMzml fixture(xml);
        XML_File file(fixture.path.c_str(), mzML);
        assert(file.defective, "Malformed mzML was not rejected");
        file.free_linknodes();
        file.free_linknodes();
    }
}

static void test_bad_arrays_fail_cleanly()
{
    const std::string mz_array = binary_array(std::vector<double>{100, 101, 102}, false, true);
    const std::string intensity_array = binary_array(std::vector<float>{10, 20, 30}, false, false);
    const std::vector<std::string> invalid_arrays = {
        mz_array, mz_array + mz_array + intensity_array,
        "<binaryDataArray><cvParam accession='MS:1000514'/><binary>AA==</binary></binaryDataArray>" + intensity_array,
        binary_array(std::vector<double>{100, 101}, false, true) + intensity_array,
        mz_array + binary_array(std::vector<float>{10, 20}, true, false)};
    for (const std::string &arrays : invalid_arrays)
    {
        TemporaryMzml fixture(test_document(test_spectrum(arrays)));
        XML_File file(fixture.path.c_str(), mzML);
        assert(!file.defective, "Binary payloads must be checked when decoding, not eagerly in the constructor");
        std::vector<float> mz = {1};
        std::vector<float> intensity = {1};
        assert(get_spectrum(&file, &mz, &intensity, 0) != 0, "Malformed array was accepted");
        assert(mz.empty() && intensity.empty(), "Partial results must be cleared after an array error");
    }

    TemporaryMzml fixture(test_document(test_spectrum("", "1", "MS:1000130", "1", "UO:0000010", "0")));
    XML_File empty_scan(fixture.path.c_str(), mzML);
    std::vector<float> mz;
    std::vector<float> intensity;
    assert(!empty_scan.defective && get_spectrum(&empty_scan, &mz, &intensity, 0) == 0 && mz.empty() && intensity.empty(),
           "An empty scan should be represented by empty arrays");
}

static void test_decoder_errors_and_uncompressed_float()
{
    std::vector<float> output;
    const std::vector<float> values = {10, 20, 30};
    const std::string plain = encoded_values(values, false);
    assert(decodeSpectrum(plain.data(), plain.size(), &output, values.size(), false, false) == 0 && output == values,
           "Uncompressed float decoding failed");
    assert(decodeSpectrum("!!!!", 4, &output, 3, false, false) != 0 && output.empty(),
           "Invalid base64 was accepted");
    assert(decode_base64("!!!!", 4).empty(), "Invalid base64 must return an empty error result");
    assert(decodeSpectrum(plain.data(), plain.size(), &output, 4, false, false) != 0 && output.empty(),
           "A truncated uncompressed array was accepted");
    assert(decodeSpectrum(plain.data(), plain.size(), &output, 2, false, false) != 0 && output.empty(),
           "An oversized uncompressed array was accepted");
    assert(decodeSpectrum(plain.data(), plain.size(), &output, 3, true, false) != 0 && output.empty(),
           "Corrupt zlib input was accepted");
    const std::string compressed = encoded_values(values, true);
    assert(decodeSpectrum(compressed.data(), compressed.size(), &output, 4, true, false) != 0 && output.empty(),
           "A decompressed array with the wrong length was accepted");
    assert(decodeSpectrum(plain.data(), plain.size(), &output, std::numeric_limits<size_t>::max(), false, true) != 0,
           "An overflowing declared length was accepted");
    const std::string nonfinite = encoded_values(std::vector<double>{std::numeric_limits<double>::infinity()}, false);
    assert(decodeSpectrum(nonfinite.data(), nonfinite.size(), &output, 1, false, true) != 0 && output.empty(),
           "Non-finite input was accepted");

    const std::vector<double> closeMasses{100, 100.00000001, 100.00000003};
    const std::string precise = encoded_values(closeMasses, true);
    std::vector<double> decodedMasses;
    assert(decodeSpectrum(precise.data(), precise.size(), &decodedMasses, closeMasses.size(), true, true) == 0 &&
               decodedMasses == closeMasses,
           "The profile decoder must retain distinct m/z values below float resolution");
}

int main()
{
    test_spec_decode_nocomp_double();
    test_spec_decode_iscomp_double();
    test_spec_decode_iscomp_float();
    test_independent_array_encodings_and_roots();
    test_selected_retention_times();
    test_referenced_array_parameters();
    test_bad_documents_fail_cleanly();
    test_bad_arrays_fail_cleanly();
    test_decoder_errors_and_uncompressed_float();
    return 0;
}

// python code used to generate test string (LLM code):

// import base64
// import struct

// doubles = [1.5, 2.7, 3.14, 4.0, 5.999]

// # Convert doubles to bytes
// # 'd' format character represents a double (8 bytes each)
// bytes_data = struct.pack(f'{len(doubles)}d', *doubles)

// # Encode to base64
// base64_string = base64.b64encode(bytes_data).decode('utf-8')

// print(base64_string)
