import argparse
import gzip
import sys


def main():
    parser = ArgParser(sys.argv)

    infile_path = parser.infile_path
    outfile_path = parser.outfile_path

    with open_file(infile_path, "rt") as varsum_file:
        with open_file(outfile_path, "wt") as outfile:
            for line in varsum_file:
                line_list = line.rstrip().split("\t")
                if skip_varsum_line(line_list):
                    continue

                output_lines = get_output_lines(line_list)
                for output_line in output_lines:
                    outfile.write(output_line + "\n")


def get_output_lines(line_list):
    start = line_list[19]
    stop = line_list[20]

    if int(stop) - int(start) > 50000:
        return []

    chrom = line_list[18]
    ref = line_list[21] if line_list[21] != "na" else "-"
    alt = line_list[22] if line_list[22] != "na" else "-"
    strand = "+"
    type = line_list[1]
    svtype = get_svtype(type)
    rcvs = line_list[11].split(";")
    rs = "rs" + line_list[9] if line_list[9] != "-1" else "-1"
    nsv = line_list[10] if line_list[10] != "-" else "-1"
    ncbi_geneid = line_list[3]

    output_lines = []
    for rcv in rcvs:
        output_line = build_output_line(chrom, start, stop, ref, alt, strand, svtype, rcv, rs, nsv, ncbi_geneid, type)
        output_lines.append(output_line)

    return output_lines


def build_output_line(chrom, start, stop, ref, alt, strand, svtype, rcv, rs, nsv, ncbi_geneid, type):
    output_line_list = [chrom, start, stop, ref, alt, strand, svtype, rs, rcv, ncbi_geneid, nsv, type]
    output_line = "\t".join(output_line_list)

    return output_line


def get_svtype(type):
    type_to_svtype_dict = {"deletion": "DEL", "insertion": "INS", "duplication": "DUP"}
    if type.lower() in type_to_svtype_dict:
        return type_to_svtype_dict[type]
    else:
        return "INS"


def skip_varsum_line(line_list):
    if line_list[0].startswith("#"):
        return True

    clin_sig = line_list[5]
    if "pathogenic" not in clin_sig.lower():
        return True

    assembly = line_list[12]
    if assembly.lower() != "grch38":
        return True

    return False


def open_file(file_path, mode):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


class ArgParser:
    def __init__(self, argv):
        description = """
                Script for extracting the coordinates, variant IDs, and NCBI gene ID for ClinVar
                records from a ClinVar variant_summary file.Output is a tab separated file with one
                ClinVar record per line, and containing the columns: chromosome, start pos,
                stop pos, reference and alternate alleles, strand, RS ID, RCV ID, NCBI gene ID, NSV ID.
                If a column doesn't have a value then a value of '-1' is used.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="infile_path", required=True,
                            help="path to variant summary file from ClinVar")
        parser.add_argument("-o", dest="outfile_path", required=True,
                            help="Tab separated output file, with one ClinVar record data per " +
                                 "line. '-1' represents an empty value.")

        args = parser.parse_args(args=argv[1:])

        self.infile_path = args.infile_path
        self.outfile_path = args.outfile_path


if __name__ == '__main__':
    main()
