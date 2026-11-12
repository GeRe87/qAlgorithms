#include "qalgorithms_read_file.h"
#include "qalgorithms_datatypes.h"

#include "../external/pugixml/pugixml.hpp"
#include "../external/simdutf/simdutf.h" // use a fast base64 decode function that makes proper use of SIMD
#include "qalgorithms_utils.h"

// #include <filesystem>
#include <cassert>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

// CMake selects both the header and the matching compression library.
#ifdef QALGORITHMS_USE_ZLIB_NG
    #include <zlib-ng.h>
#else
    #include <zlib.h>
inline int32_t zng_uncompress(uint8_t *dest, size_t *destLen, const uint8_t *source, size_t sourceLen)
{
    uLongf destLen2 = *destLen;
    uint32_t z_ret = uncompress(dest, &destLen2, source, sourceLen);
    *destLen = destLen2;
    return z_ret;
}
inline size_t zng_compressBound(size_t sourceLen)
{
    return compressBound(sourceLen);
}
inline int32_t zng_compress(uint8_t *dest, size_t *destLen, const uint8_t *source, size_t sourceLen)
{
    uLongf destLen2 = *destLen;
    uint32_t z_ret = compress(dest, &destLen2, source, sourceLen);
    *destLen = destLen2;
    return z_ret;
}
#endif
// The standard zlib API is the portable default; native zlib-ng is optional.

namespace qAlgorithms
{
    static bool isCentroided_fun(const XML_File *file);

    static bool read_size(const char *text, size_t *value)
    {
        if (!text || !*text)
            return false;
        const char *end = text + std::strlen(text);
        const auto parsed = std::from_chars(text, end, *value);
        return parsed.ec == std::errc() && parsed.ptr == end;
    }

    static pugi::xml_node spectrum_node(const XML_File *file, const size_t index)
    {
        if (!file || file->defective || !file->linknodes || index >= file->linknodes->size())
            return {};
        return file->linknodes->at(index);
    }

    // CV names are descriptive text; the accession is the stable identifier.
    // mzML also permits these parameters to be supplied by a shared group.
    static pugi::xml_node cv_param(const pugi::xml_node &node, const char *accession)
    {
        pugi::xml_node param = node.find_child_by_attribute("cvParam", "accession", accession);
        if (param)
            return param;
        pugi::xml_node root = node.root().child("mzML");
        if (!root)
            root = node.root().child("indexedmzML").child("mzML");
        for (pugi::xml_node ref : node.children("referenceableParamGroupRef"))
        {
            const pugi::xml_node group = root.child("referenceableParamGroupList")
                                             .find_child_by_attribute("referenceableParamGroup", "id",
                                                                      ref.attribute("ref").value());
            param = group.find_child_by_attribute("cvParam", "accession", accession);
            if (param)
                return param;
        }
        return {};
    }

    static bool binary_encoding(const pugi::xml_node &array, bool *compressed, bool *f64)
    {
        const bool is_compressed = bool(cv_param(array, "MS:1000574"));
        const bool is_uncompressed = bool(cv_param(array, "MS:1000576"));
        const bool is_double = bool(cv_param(array, "MS:1000523"));
        const bool is_float = bool(cv_param(array, "MS:1000521"));
        if (is_compressed == is_uncompressed || is_double == is_float ||
            cv_param(array, "MS:1002312") || cv_param(array, "MS:1002313") || cv_param(array, "MS:1002314"))
            return false;
        *compressed = is_compressed;
        *f64 = is_double;
        return true;
    }

    static Polarities get_polarity_mode(const XML_File *file)
    {
        Polarities mode = Polarities::unknown_polarity;
        for (size_t specNum = 0; specNum < file->linknodes->size(); ++specNum)
        {
            if (spectrum_ms_level(file, specNum) != 1 ||
                (!file->isCentroided && !spectrum_is_profile(file, specNum)))
                continue;
            const Polarities polarity = spectrum_polarity(file, specNum);
            if (polarity == Polarities::unknown_polarity)
                return polarity;
            if (mode == Polarities::unknown_polarity)
                mode = polarity;
            else if (polarity != mode)
                return Polarities::mixed;
        }
        return mode;
    };

    XML_File::XML_File(const path_char *file, const SourceFileType type)
    {
        filetype = type;
        if (filetype != mzML)
        {
            defective = true;
            (void)fprintf(stderr, "Error: only mzML input is supported.\n");
            return;
        }

        pugi::xml_parse_result loading_result = mzml_base_document.load_file(file);

        if (loading_result)
        {
            mzml_root_node = mzml_base_document.document_element();
        }
        else
        {
            (void)fprintf(stderr, "Error: .mzML file could not be opened. Error description:\n%s\n",
                          loading_result.description());
            defective = true;
            return;
        }

        if (std::strcmp(mzml_root_node.name(), "indexedmzML") == 0)
            mzml_root_node = mzml_root_node.child("mzML");
        if (std::strcmp(mzml_root_node.name(), "mzML") != 0)
        {
            defective = true;
            (void)fprintf(stderr, "Error: expected an mzML document, optionally wrapped in indexedmzML.\n");
            return;
        }

        pugi::xml_node spec_list = mzml_root_node.child("run").child("spectrumList");
        size_t declared_count = 0;
        if (!spec_list || !read_size(spec_list.attribute("count").value(), &declared_count) ||
            declared_count > std::numeric_limits<uint32_t>::max())
        {
            defective = true;
            (void)fprintf(stderr, "Error: missing spectrumList or invalid spectrum count in mzML.\n");
            return;
        }
        if (declared_count == 0)
        {
            defective = true;
            (void)fprintf(stderr, "Error: the supplied mzML file does not contain any spectra.\n");
            return;
        }

        // Do not allocate an input-supplied count before checking the actual nodes.
        linknodes = new std::vector<pugi::xml_node>;
        for (pugi::xml_node child : spec_list.children("spectrum"))
        {
            linknodes->push_back(child);
        }
        if (linknodes->size() != declared_count)
        {
            defective = true;
            (void)fprintf(stderr, "Error: spectrumList count does not match the number of spectra.\n");
            return;
        }
        number_spectra = static_cast<uint32_t>(declared_count);

        for (size_t index = 0; index < number_spectra; ++index)
        {
            const uint32_t level = spectrum_ms_level(this, index);
            const pugi::xml_node spec = linknodes->at(index);
            size_t array_length = 0;
            if (level == 0 || !read_size(spec.attribute("defaultArrayLength").value(), &array_length))
            {
                defective = true;
                (void)fprintf(stderr, "Error: spectrum %zu has an invalid MS level or array length.\n", index);
                return;
            }
            if (level == 1)
            {
                const bool profile = bool(cv_param(spec, "MS:1000128"));
                const bool centroid = bool(cv_param(spec, "MS:1000127"));
                if (profile == centroid || spectrum_polarity(this, index) == Polarities::unknown_polarity ||
                    !std::isfinite(spectrum_rt(this, index)))
                {
                    defective = true;
                    (void)fprintf(stderr, "Error: MS1 spectrum %zu needs a valid spectrum type, polarity and scan start time.\n", index);
                    return;
                }
            }
        }

        isCentroided = isCentroided_fun(this);
        polarityMode = get_polarity_mode(this);
        if (polarityMode == Polarities::unknown_polarity)
        {
            defective = true;
            (void)fprintf(stderr, "Error: no MS1 spectra with a known polarity were found.\n");
            return;
        }
        zlib_compression = spectrum_is_compressed(this, 0);
        precision_f64 = spectrum_is_float64(this, 0);
    };

    XML_File::~XML_File()
    {
        delete linknodes;
    }

    void XML_File::free_linknodes()
    {
        delete linknodes;
        linknodes = nullptr;
        defective = true;
    };

    template <typename T>
    static int32_t decodeSpectrum(const char *const binaryData,
                                  const size_t encodedLength,
                                  std::vector<T> *spectrum,
                                  const size_t finalSize,
                                  bool compression, bool f64)
    {
        spectrum->clear();
        const size_t element_size = f64 ? sizeof(double) : sizeof(float);
        if (finalSize > spectrum->max_size() || finalSize > std::numeric_limits<size_t>::max() / element_size)
            return 1;
        if (finalSize == 0)
            return encodedLength == 0 ? 0 : 1;
        const size_t expected_bytes = finalSize * element_size;
        std::vector<char> decoded_string = decode_base64(binaryData, encodedLength);
        if (decoded_string.empty())
            return 1;
        if (compression)
        {
            std::vector<char> buffer(expected_bytes);
            size_t buffer_size = buffer.size();
            const int32_t status = zng_uncompress(reinterpret_cast<uint8_t *>(buffer.data()), &buffer_size,
                                                   reinterpret_cast<const uint8_t *>(decoded_string.data()),
                                                   decoded_string.size());
            if (status != Z_OK || buffer_size != expected_bytes)
                return 2;
            decoded_string.swap(buffer);
        }
        else if (decoded_string.size() != expected_bytes)
        {
            return 3;
        }

        spectrum->resize(finalSize);
        for (size_t index = 0; index < finalSize; ++index)
        {
            // memcpy avoids alignment and aliasing assumptions about the byte buffer.
            double value = 0;
            if (f64)
                std::memcpy(&value, decoded_string.data() + index * element_size, sizeof(double));
            else
            {
                float value_float = 0;
                std::memcpy(&value_float, decoded_string.data() + index * element_size, sizeof(float));
                value = value_float;
            }
            if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<T>::max())
            {
                spectrum->clear();
                return 4;
            }
            spectrum->at(index) = static_cast<T>(value);
        }
        return 0;
    }

    template <typename T>
    static int32_t readSpectrum(const XML_File *file, std::vector<T> *spectrum_mz,
        std::vector<T> *spectrum_int, size_t index)
    {
        spectrum_mz->clear();
        spectrum_int->clear();
        const pugi::xml_node spec = spectrum_node(file, index);
        size_t spectrum_size = 0;
        if (!spec || !read_size(spec.attribute("defaultArrayLength").value(), &spectrum_size))
        {
            (void)fprintf(stderr, "Error: spectrum %zu has an invalid array length or index.\n", index);
            return 1;
        }
        if (spectrum_size == 0)
            return 0;

        pugi::xml_node mz_array;
        pugi::xml_node intensity_array;
        for (pugi::xml_node array : spec.child("binaryDataArrayList").children("binaryDataArray"))
        {
            const bool is_mz = bool(cv_param(array, "MS:1000514"));
            const bool is_intensity = bool(cv_param(array, "MS:1000515"));
            if ((is_mz && is_intensity) || (is_mz && mz_array) || (is_intensity && intensity_array))
            {
                (void)fprintf(stderr, "Error: spectrum %zu has ambiguous m/z or intensity arrays.\n", index);
                return 1;
            }
            if (is_mz)
                mz_array = array;
            if (is_intensity)
                intensity_array = array;
        }

        const pugi::xml_node arrays[] = {mz_array, intensity_array};
        std::vector<T> *outputs[] = {spectrum_mz, spectrum_int};
        for (size_t array_index = 0; array_index < 2; ++array_index)
        {
            const pugi::xml_node array = arrays[array_index];
            bool compressed = false;
            bool f64 = false;
            size_t array_length = spectrum_size;
            const pugi::xml_attribute length_attribute = array.attribute("arrayLength");
            if (!array || !array.child("binary") || !binary_encoding(array, &compressed, &f64) ||
                (length_attribute && (!read_size(length_attribute.value(), &array_length) || array_length != spectrum_size)))
            {
                (void)fprintf(stderr, "Error: spectrum %zu has missing or unsupported %s array metadata; "
                                      "32/64-bit floats with zlib or no compression are supported.\n",
                              index, array_index == 0 ? "m/z" : "intensity");
                spectrum_mz->clear();
                spectrum_int->clear();
                return 1;
            }
            const char *binary_data = array.child("binary").child_value();
            if (decodeSpectrum(binary_data, std::strlen(binary_data), outputs[array_index], spectrum_size, compressed, f64) != 0)
            {
                (void)fprintf(stderr, "Error: spectrum %zu has invalid %s binary data, decompression or decoded length.\n",
                              index, array_index == 0 ? "m/z" : "intensity");
                spectrum_mz->clear();
                spectrum_int->clear();
                return 2;
            }
        }
        return 0;
    }

    int32_t get_spectrum(const XML_File *file, std::vector<float> *spectrum_mz,
        std::vector<float> *spectrum_int, size_t index)
    {
        return readSpectrum(file, spectrum_mz, spectrum_int, index);
    }

    int32_t get_spectrum(const XML_File *file, std::vector<double> *spectrum_mz,
        std::vector<double> *spectrum_int, size_t index)
    {
        return readSpectrum(file, spectrum_mz, spectrum_int, index);
    }

    void get_spectra_RT(const XML_File *file, const std::vector<uint32_t> *indices, std::vector<float> *const RTs)
    {
        RTs->clear();
        RTs->reserve(indices->size());
        for (uint32_t index : *indices)
        {
            const float rt = spectrum_rt(file, index);
            if (!std::isfinite(rt) || (!RTs->empty() && rt <= RTs->back()))
            {
                (void)fprintf(stderr, "Error: selected spectrum %u has an invalid or non-increasing retention time.\n", index);
                RTs->clear();
                return;
            }
            RTs->push_back(rt);
        }
    };

    std::vector<uint32_t> filter_spectra(const XML_File *file,
                                         const bool ms1,
                                         const Polarities polarity,
                                         const bool profile_mode)
    {
        // return a vector of all indices that are relevant to the query. Properties are checked in order of regularity.
        if (!file || file->defective || !file->linknodes)
            return {};
        const size_t numSpectra = file->number_spectra;
        std::vector<uint32_t> indices;
        indices.reserve(numSpectra);

        for (uint32_t specNum = 0; specNum < numSpectra; specNum++)
        {
            bool isProfile = spectrum_is_profile(file, specNum);
            if (isProfile != profile_mode) // this does not allow for processing of partially centroided data - change?
                continue;

            Polarities polarity_spec = spectrum_polarity(file, specNum);
            if (polarity_spec != polarity)
                continue;

            uint32_t level = spectrum_ms_level(file, specNum);
            bool isMS1 = 1 == level;
            if (isMS1 != ms1)
                continue; // only ms1 or msn data can be retrieved at once.

            indices.push_back(specNum);
        }
        indices.shrink_to_fit();
        return indices;
    }

    static bool isCentroided_fun(const XML_File *file)
    {
        size_t centroided = 0;
        size_t profile = 0;

        for (size_t specNum = 0; specNum < file->number_spectra; ++specNum)
        {
            uint32_t level = spectrum_ms_level(file, specNum);
            if (level != 1)
                continue;

            if (spectrum_is_profile(file, specNum))
            {
                profile += 1;
            }
            else
            {
                centroided += 1;
            }
        }

        return centroided > 0 && profile == 0;
    }

    bool spectrum_is_compressed(const XML_File *file, const size_t specNum)
    {
        const pugi::xml_node array = spectrum_node(file, specNum).child("binaryDataArrayList").child("binaryDataArray");
        return bool(cv_param(array, "MS:1000574"));
    }

    bool spectrum_is_float64(const XML_File *file, const size_t specNum)
    {
        const pugi::xml_node array = spectrum_node(file, specNum).child("binaryDataArrayList").child("binaryDataArray");
        return bool(cv_param(array, "MS:1000523"));
    }

    bool spectrum_is_profile(const XML_File *file, const size_t specNum)
    {
        return bool(cv_param(spectrum_node(file, specNum), "MS:1000128"));
    }

    uint32_t spectrum_ms_level(const XML_File *file, const size_t specNum)
    {
        const pugi::xml_node param = cv_param(spectrum_node(file, specNum), "MS:1000511");
        size_t level = 0;
        if (!read_size(param.attribute("value").value(), &level) || level > std::numeric_limits<uint32_t>::max())
            return 0;
        return static_cast<uint32_t>(level);
    }

    Polarities spectrum_polarity(const XML_File *file, const size_t specNum)
    {
        const pugi::xml_node spec = spectrum_node(file, specNum);
        const bool positive_scan = bool(cv_param(spec, "MS:1000130"));
        const bool negative_scan = bool(cv_param(spec, "MS:1000129"));
        if (positive_scan && !negative_scan)
            return Polarities::positive;
        if (negative_scan && !positive_scan)
            return Polarities::negative;

        return Polarities::unknown_polarity;
    }

    float spectrum_rt(const XML_File *file, const size_t specNum)
    {
        const pugi::xml_node scan = spectrum_node(file, specNum).child("scanList").child("scan");
        const pugi::xml_node rt_node = cv_param(scan, "MS:1000016");
        const char *value_text = rt_node.attribute("value").value();
        char *value_end = nullptr;
        errno = 0;
        const float value = std::strtof(value_text, &value_end);
        if (value_end == value_text || *value_end != '\0' || errno == ERANGE || !std::isfinite(value) || value < 0)
            return std::numeric_limits<float>::quiet_NaN();

        const char *unit_accession = rt_node.attribute("unitAccession").value();
        const char *unit_name = rt_node.attribute("unitName").value();
        const bool minutes = *unit_accession ? std::strcmp(unit_accession, "UO:0000031") == 0
                                            : std::strcmp(unit_name, "minute") == 0;
        const bool seconds = *unit_accession ? std::strcmp(unit_accession, "UO:0000010") == 0
                                            : std::strcmp(unit_name, "second") == 0;
        if (!minutes && !seconds)
            return std::numeric_limits<float>::quiet_NaN();
        const float rt = value * (minutes ? 60.0F : 1.0F);
        return std::isfinite(rt) ? rt : std::numeric_limits<float>::quiet_NaN();
    }

    // Decodes a Base64 string into a string with binary data using the simdutf library subset chosen by '--with-base64'
    // (https://github.com/simdutf/simdutf/tree/master?tab=readme-ov-file#single-header-version-with-limited-features).
    std::vector<char> decode_base64(const char *encoded_string, const size_t encoded_length)
    {
        if (!encoded_string || encoded_length == 0)
            return {};
        const size_t decoded_length = simdutf::maximal_binary_length_from_base64(encoded_string, encoded_length);
        if (decoded_length == 0)
            return {};
        std::vector<char> output(decoded_length);
        simdutf::result simd_res = simdutf::base64_to_binary(encoded_string,
                                                             encoded_length,
                                                             output.data(), simdutf::base64_default,
                                                             simdutf::last_chunk_handling_options::strict);

        if (simd_res.error != 0) // [[unlikely]]
        {
            return {}; // error message is handled one function above
        }
        output.resize(simd_res.count);
        return output;
    };

    std::vector<char> encode_base64_dbl(const double *input_dbl, const size_t in_len_dbl)
    {
        const size_t length = in_len_dbl * sizeof(double);
        const char *input = (const char *)input_dbl;
        std::vector<char> buffer(simdutf::base64_length_from_binary(length));
        size_t written = simdutf::binary_to_base64(input, length, buffer.data());
        assert(written == buffer.size());
        return buffer;
    }

    void compress_zlib(const std::vector<char> *uncompressed_string, std::vector<char> *output_string)
    {
        output_string->resize(zng_compressBound(uncompressed_string->size()));
        size_t outSize = output_string->size();
        const int32_t status = zng_compress(reinterpret_cast<uint8_t *>(output_string->data()), &outSize,
                                            reinterpret_cast<const uint8_t *>(uncompressed_string->data()),
                                            uncompressed_string->size());
        if (status != Z_OK)
        {
            output_string->clear();
            return;
        }
        output_string->resize(outSize);
    };

    void decompress_inPlace(std::vector<char> *in_out, const size_t decompressedSize)
    {
        std::vector<char> buffer(decompressedSize);
        size_t outSize = decompressedSize;
        const int32_t status = zng_uncompress(reinterpret_cast<uint8_t *>(buffer.data()), &outSize,
                                              reinterpret_cast<const uint8_t *>(in_out->data()), in_out->size());
        if (status != Z_OK || outSize != decompressedSize)
        {
            in_out->clear();
            return;
        }
        in_out->swap(buffer);
    }

    void compress_and_encode(const uint8_t *input_arr,
                             const size_t length,
                             std::vector<char> *output_string)
    {
        // compress using zlib-ng
        const size_t compressionBound = zng_compressBound(length);
        std::vector<char> buffer(compressionBound);
        size_t outSize = compressionBound;
        const int32_t status = zng_compress(reinterpret_cast<uint8_t *>(buffer.data()), &outSize, input_arr, length);
        if (status != Z_OK || outSize > buffer.size())
        {
            output_string->clear();
            return;
        }

        // encode as base64
        size_t base64_size = simdutf::base64_length_from_binary(outSize);
        output_string->resize(base64_size);
        const size_t written = simdutf::binary_to_base64(buffer.data(), outSize, output_string->data());
        assert(written == base64_size);
    }

} // namespace qAlgorithms
