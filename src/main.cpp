#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "sequence.hpp"
#include "subsampler.hpp"

#include "bioparser/bioparser.hpp"

static struct option options[] = {
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    std::vector<std::string> input_parameters;

    char argument;
    while ((argument = getopt_long(argc, argv, "h", options, nullptr)) != -1) {
        switch (argument) {
            case 'h':
            default:
                help();
                exit(1);
        }
    }

    for (int32_t i = optind; i < argc; ++i) {
        input_parameters.emplace_back(argv[i]);
    }

    if (input_parameters.size() < 3) {
        fprintf(stderr, "[rampler::] error: missing input parameter(s)!\n");
        help();
        exit(1);
    }

    auto subsampler = rampler::createSubsampler(input_parameters[0],
        atoi(input_parameters[1].c_str()));
    subsampler->initialize();

    std::string base_name = input_parameters[0].substr(input_parameters[0].find_last_of("/\\") + 1,
        input_parameters[0].rfind('.'));
    std::string extension = input_parameters[0].substr(input_parameters[0].rfind('.'));

    for (uint32_t i = 2; i < input_parameters.size(); ++i) {
        uint32_t coverage = atoi(input_parameters[i].c_str());

        std::vector<std::unique_ptr<rampler::Sequence>> sequences;
        subsampler->subsample(sequences, coverage);

        if (sequences.empty()) {
            continue;
        }

        std::string out_path = base_name + "_" + input_parameters[i] + "x" +
            extension;
        auto out = fopen(out_path.c_str(), "w");

        if (sequences[0]->quality().empty()) {
            for (const auto& it: sequences) {
                fprintf(out, ">%s\n%s\n", it->name().c_str(), it->data().c_str());
            }
        } else {
            for (const auto& it: sequences) {
                fprintf(out, "@%s\n%s\n+\n%s\n", it->name().c_str(), it->data().c_str(),
                    it->quality().c_str());
            }
        }

        fclose(out);
    }

    return 0;
}

void help() {
    printf(
        "usage: rampler [options ...] <sequences> <reference length> <coverage> [<coverage> ...]\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format containing sequences to be\n"
        "        subsampled\n"
        "    <reference length>\n"
        "        integer denoting length of the reference genome (or assembly)\n"
        "        from which the sequences originate\n"
        "    <coverage>\n"
        "        integer denoting desired coverage of the subsampled sequences\n"
        "\n"
        "    options:\n"
        "        -h, --help\n"
        "            prints out the help\n");
}
