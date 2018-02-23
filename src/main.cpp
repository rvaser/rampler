#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "sequence.hpp"
#include "sampler.hpp"

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
                exit(0);
        }
    }

    if (optind == argc) {
        help();
        exit(1);
    }

    for (int32_t i = optind; i < argc; ++i) {
        input_parameters.emplace_back(argv[i]);
    }

    bool do_subsample = false, do_split = false;
    if (input_parameters[0] == "subsample") {
        do_subsample = true;
    } else if (input_parameters[0] == "split") {
        do_split = true;
    } else {
        fprintf(stderr, "[rampler::] error: unkown mode %s!\n", input_parameters[0].c_str());
        exit(1);
    }

    if ((do_subsample && input_parameters.size() < 3) ||
        (do_split && input_parameters.size() < 2)) {

        fprintf(stderr, "[rampler::] error: missing input parameter(s)!\n");
        help();
        exit(1);
    }

    auto sampler = rampler::createSampler(input_parameters[1]);
    sampler->initialize();

    if (do_split) {
        sampler->split(atoi(input_parameters[2].c_str()));
    } else if (do_subsample) {
        uint32_t reference_length = atoi(input_parameters[2].c_str());
        for (uint32_t i = 3; i < input_parameters.size(); ++i) {
            sampler->subsample(reference_length, atoi(input_parameters[i].c_str()));
        }
    }

    return 0;
}

void help() {
    printf(
        "usage: rampler [options ...] <mode>\n"
        "\n"
        "    <mode>\n"
        "        subsample <sequences> <reference length> <coverage> [<coverage> ...]\n"
        "\n"
        "            <sequences>\n"
        "                input file in FASTA/FASTQ format containing sequences\n"
        "                to be subsampled\n"
        "            <reference length>\n"
        "                integer denoting length of the reference genome (or\n"
        "                assembly) from which the sequences originate\n"
        "            <coverage>\n"
        "                integer denoting desired coverage of the subsampled\n"
        "                sequences\n"
        "\n"
        "        split <sequences> <chunk size>\n"
        "\n"
        "            <sequences>\n"
        "                input file in FASTA/FASTQ format containing sequences\n"
        "                which will be split into smaller chunks\n"
        "            <chunk size>\n"
        "                integer denoting the desired chunk size in bytes\n"
        "\n"
        "    options:\n"
        "        -h, --help\n"
        "            prints out the help\n");
}
