/*!
 * @file subsampler.cpp
 *
 * @brief Subsampler class source file
 */

#include <random>

#include "sequence.hpp"
#include "subsampler.hpp"

#include "bioparser/bioparser.hpp"

namespace rampler {

constexpr uint32_t kChunkSize = 1024 * 1024 * 1024; // ~ 1GB

template<class T>
void shrinkToFit(std::vector<std::unique_ptr<T>>& src, uint64_t begin) {

    uint64_t i = begin;
    for (uint64_t j = begin; i < src.size(); ++i) {
        if (src[i] != nullptr) {
            continue;
        }

        j = std::max(j, i);
        while (j < src.size() && src[j] == nullptr) {
            ++j;
        }

        if (j >= src.size()) {
            break;
        } else if (i != j) {
            src[i].swap(src[j]);
        }
    }
    if (i < src.size()) {
        src.resize(i);
    }
}

std::unique_ptr<Subsampler> createSubsampler(const std::string& sequences_path,
    uint32_t reference_length) {

    std::unique_ptr<bioparser::Parser<Sequence>> sparser = nullptr;

    auto extension = sequences_path.substr(std::min(sequences_path.rfind('.'),
        sequences_path.size()));
    if (extension == ".fasta" || extension == ".fa") {
        sparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            sequences_path);
    } else if (extension == ".fastq" || extension == ".fq") {
        sparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            sequences_path);
    } else {
        fprintf(stderr, "[rampler::createSubsampler] error: "
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fa, .fastq, .fq)!\n", sequences_path.c_str());
        exit(1);
    }

    if (reference_length == 0) {
        fprintf(stderr, "[rampler::createSubsampler] error: "
            "invalid reference length\n");
        exit(1);
    }

    return std::unique_ptr<Subsampler>(new Subsampler(std::move(sparser),
        reference_length));
}

Subsampler::Subsampler(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
    uint32_t reference_length)
        : sparser_(std::move(sparser)), reference_length_(reference_length),
        sequences_length_(0) {
}

Subsampler::~Subsampler() {
}

void Subsampler::initialize() {

    if (sequences_length_ != 0) {
        fprintf(stderr, "[rampler::Subsampler::initialize] warning: "
            "object already initialized!\n");
        return;
    }

    sparser_->reset();
    while (true) {
        std::vector<std::unique_ptr<Sequence>> sequences;
        auto status = sparser_->parse_objects(sequences, kChunkSize);

        for (const auto& it: sequences) {
            sequences_length_ += it->data().size();
        }

        if (!status) {
            break;
        }
    }
}

void Subsampler::subsample(std::vector<std::unique_ptr<Sequence>>& dst,
    uint32_t coverage) {

    if (coverage * reference_length_ > sequences_length_) {
        fprintf(stderr, "[rampler::Subsampler::subsample] warning: "
            "insufficient data for coverage %u!\n", coverage);
        return;
    }

    std::random_device r;
    std::mt19937 generator(r());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    double ratio = (coverage * reference_length_) / static_cast<double>(
        sequences_length_);

    sparser_->reset();
    while (true) {
        uint64_t l = dst.size();
        auto status = sparser_->parse_objects(dst, kChunkSize);

        for (uint64_t i = l; i < dst.size(); ++i) {
            if (distribution(generator) > ratio) {
                dst[i].reset();
            }
        }

        shrinkToFit(dst, l);

        if (!status) {
            break;
        }
    }
}

}
