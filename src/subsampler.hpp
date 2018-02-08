/*!
 * @file subsampler.hpp
 *
 * @brief Subsampler class header file
 */

#pragma once

#include <stdlib.h>
#include <vector>
#include <memory>

namespace bioparser {
    template<class T>
    class Parser;
}

namespace rampler {

class Sequence;

class Subsampler;
std::unique_ptr<Subsampler> createSubsampler(const std::string& sequences_path,
    uint32_t reference_length);

class Subsampler {
public:
    ~Subsampler();

    void initialize();

    void subsample(std::vector<std::unique_ptr<Sequence>>& dst, uint32_t coverage);

    friend std::unique_ptr<Subsampler> createSubsampler(
        const std::string& sequences_path, uint32_t reference_length);
private:
    Subsampler(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
        uint32_t reference_length);
    Subsampler(const Subsampler&) = delete;
    const Subsampler& operator=(const Subsampler&) = delete;

    std::unique_ptr<bioparser::Parser<Sequence>> sparser_;
    uint32_t reference_length_;
    uint64_t sequences_length_;
};

}
