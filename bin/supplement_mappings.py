import argparse
import sys


# To supplement older trait to ontology mappings with older mappings
# Both old and current mappings files should have trait in first column


def main():
    parser = ArgParser(sys.argv)

    with open(parser.outfile_path, "wt") as outfile:
        curr_traits = process_current_mappings(parser.curr_file_path, outfile)
        add_prev_mappings(parser.prev_file_path, outfile, curr_traits)


class ArgParser:
    def __init__(self, argv):
        parser = argparse.ArgumentParser()

        parser.add_argument("-l", dest="old_maps_filepath", required=True)
        parser.add_argument("-c", dest="curr_maps_filepath", required=True)
        parser.add_argument("-o", dest="outfile_path", required=True)

        args = parser.parse_args(args=argv[1:])

        self.old_maps_filepath = args.old_maps_filepath
        self.curr_maps_filepath = args.curr_maps_filepath
        self.outfile_path = args.outfile_path


def process_current_mappings(file_path, outfile_handle):
    curr_traits = set()
    with open(file_path, "rt") as f:
        for line in f:
            line_list = line.rstrip().split("\t")
            trait = line_list[0].lower()
            curr_traits.add(trait)

            outfile_handle.write(line)

    return curr_traits


def add_prev_mappings(file_path, outfile_handle, curr_traits):
    with open(file_path, "rt") as f:
        for line in f:
            line_list = line.rstrip().split("\t")
            if len(line_list) < 2:
                continue

            trait = line_list[0].lower()
            if trait in curr_traits:
                continue

            outfile_handle.write(line)


if __name__ == '__main__':
    main()