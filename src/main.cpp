#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "bioparser/bioparser.hpp"

static struct option options[] = {
    {"coverage", required_argument, 0, 'c'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    std::vector<std::string> input_parameters;
    uint32_t coverage = 30;

    char argument;
    while ((argument = getopt_long(argc, argv, "c:h", options, nullptr)) != -1) {
        switch (argument) {
            case 'c':
                coverage = atoi(optarg);
                break;
            case 'h':
            default:
                help();
                exit(1);
        }
    }

    for (int32_t i = optind; i < argc; ++i) {
        input_parameters.emplace_back(argv[i]);
    }

    if (input_parameters.size() < 2) {
        fprintf(stderr, "[rampler::] error: missing input parameter(s)!\n");
        help();
        exit(1);
    }

    return 0;
}

void help() {
    printf(
        "usage: rampler [options ...] <reference length> <sequences>\n"
        "\n"
        "    <reference length>\n"
        "        integer denoting length of the reference genome (or assembly)\n"
        "        from which the sequences originate\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format containing seqeunces to be\n"
        "        subsampled\n"
        "\n"
        "    options:\n"
        "        -c, --coverage <int>\n"
        "            default: 30\n"
        "            desired coverage of the subsampled sequences\n"
        "        -h, --help\n"
        "            prints out the help\n");
}
