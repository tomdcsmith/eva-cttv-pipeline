import argparse
from collections import defaultdict
import sys


# Script that takes two tab separated files,
# one with RCV to ontology URL mappings and one with RCV to trait name mappings,
# and outputs a file with a trait name per line, with the URLs connected through the RCVs


def main():
    parser = ArgParser(sys.argv)

    rcv_to_urls = get_rcv_to_urls(parser.rcv_2_url_filepath)

    rcv_to_trait = get_rcv_to_trait(parser.rcv_to_trait_filepath)

    trait_to_urls = get_trait_to_url(rcv_to_trait, rcv_to_urls)

    with open(parser.outfile_path, "wt") as f:
        for trait, urls in trait_to_urls.items():
            url_string = '\t'.join(urls)
            out_string = '\t'.join([trait, url_string])
            f.write(out_string + "\n")


class ArgParser:
    def __init__(self, argv):
        parser = argparse.ArgumentParser()

        parser.add_argument("-u", dest="rcv_2_url_filepath", help="Tab separated file containing Clinvar RCV accessions in first column and urls in second column", required=True)
        parser.add_argument("-t", dest="rcv_2_trait_filepath", help="Tab separated file containing Clinvar RCV accessions in first column and traits in second column", required=True)
        parser.add_argument("-o", dest="outfile_path", required=True)

        args = parser.parse_args(args=argv[1:])

        self.rcv_2_url_filepath = args.rcv_2_url_filepath
        self.rcv_2_trait_filepath = args.rcv_2_trait_filepath
        self.outfile_path = args.outfile_path


def get_rcv_to_urls(file_path):
    rcv_to_urls = {}
    with open(file_path, "rt") as f:
        for line in f:
            line_list = line.rstrip().split("\t")
            rcv = line_list[0]
            urls_str = line_list[1]
            urls_set = set(urls_str.split(","))
            rcv_to_urls[rcv] = urls_set
    return rcv_to_urls


def get_rcv_to_trait(file_path):
    rcv_to_trait = {}
    with open(file_path, "rt") as f:
        for line in f:
            line_list = line.rstrip().split("\t")
            rcv = line_list[0]
            trait = line_list[1]
            rcv_to_trait[rcv] = trait.lower()
    return rcv_to_trait


def get_trait_to_url(rcv_to_trait, rcv_to_urls):
    rcv_to_trait = rcv_to_trait
    rcv_to_urls = rcv_to_urls

    trait_to_urls = defaultdict(set)
    for rcv, trait in rcv_to_trait.items():
        if rcv not in rcv_to_urls:
            continue

        new_urls = rcv_to_urls[rcv]
        trait_to_urls[trait].update(new_urls)

    return trait_to_urls


if __name__ == "__main__":
    main()
