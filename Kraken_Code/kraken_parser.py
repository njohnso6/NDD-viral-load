#!/usr/bin/env python

import argparse
from ete3 import NCBITaxa


class KrakenData:
    def __init__(self, report_line):
        self.tax_level = self.determine_tax_level(report_line.split()[3])
        self.num_reads = int(report_line.split()[1])
        self.num_reads_direct = int(report_line.split()[2])
        self.name = report_line.split()[5:]
        self.name = ' '.join(self.name).rstrip()
        self.ncbi_tax_id = int(report_line.split()[4])

    def determine_tax_level(self, tax_coding):
        if tax_coding == 'U':
            return 'Unclassified'
        elif tax_coding == 'K':
            return 'Kingdom'
        elif tax_coding == 'D':
            return 'Domain'
        elif tax_coding == 'P':
            return 'Phylum'
        elif tax_coding == 'C':
            return 'Class'
        elif tax_coding == 'O':
            return 'Order'
        elif tax_coding == 'F':
            return 'Family'
        elif tax_coding == 'G':
            return 'Genus'
        elif tax_coding == 'S':
            return 'Species'
        elif tax_coding == '-':
            return 'Other'


def taxid_to_lineage_string(taxid):
    tax_order = ['kingdom', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    outstr = ''
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    for level in tax_order:
        for tid in names:
            rank = ncbi.get_rank([tid])
            if rank[tid] == 'superkingdom':
                rank[tid] = 'domain'
            if rank[tid] == level:
                outstr += level[0] + '_' + names[tid] + ';'
    return outstr[:-1]


def determine_unassigned_rank(taxid):
    """
    Given a taxid, will use ete3 to look at all its descendants. Based on what it finds, will infer what taxonomic
    level the taxid should be at. Useful for things that have 'no rank' according to NCBI.
    :param taxid: NCBI taxid, should be an integer
    :return: string that says what taxonomy level we're at, one of the options from tax_order
    """
    tax_order = ['kingdom', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    ncbi = NCBITaxa()
    descendants = ncbi.get_descendant_taxa(taxid, intermediate_nodes=True)
    lowest_rank = 900
    for descendant in descendants:
        rank = ncbi.get_rank([descendant])
        if rank[descendant] in tax_order:
            rank_number = tax_order.index(rank[descendant])
            if rank_number < lowest_rank:
                lowest_rank = rank_number
    return tax_order[lowest_rank - 1]


if __name__ == '__main__':
    taxonomy_order = ['Unclassified', 'Kingdom', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Other']
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file',
                        type=str,
                        required=True,
                        nargs='+',
                        help='Path to one or more kraken output file(s), generated by kraken-report. Separate'
                             ' each file with a space.')
    parser.add_argument('-l', '--level',
                        type=str,
                        default='Family',
                        choices=taxonomy_order,
                        help='Taxonomy level you want a report for. Defaults to Family level.')
    parser.add_argument('-t', '--taxonomy',
                        type=str,
                        default='Bacteria',
                        help='Subset of data to look at (i.e. only look at families within Firmicutes). Defaults to'
                             ' examining all bacteria.')
    parser.add_argument('-o', '--output_file',
                        type=str,
                        required=True,
                        help='Path to output file. Will be written in CSV format.')
    parser.add_argument('-f', '--full_taxonomy',
                        default=False,
                        action='store_true',
                        help='If enabled, will output full taxonomy as headers instead of only current genus/species/'
                             'whatever else.')
    args = parser.parse_args()

    # Keep all of our data in one giant list that contains dictionaries - don't think we'll ever run into
    # large enough datasets that this is a problem.
    output_list_of_dicts = list()

    for input_file in args.input_file:
        # Dictionary will store read counts for each family/genus/whatever found, as well as the sample that we're on.
        output_dict = dict()
        output_dict['Sample'] = input_file
        # Flags! Only start writing output once we've seen our taxonomy of interest.
        write_output = False
        tax_level = None
        try:
            with open(input_file) as f:
                lines = f.readlines()
        except FileNotFoundError:
            print('ERROR: Could not find file {}: Please make sure the path to the file is correct. '
                  'Exiting...'.format(input_file))
            quit()
        check_tax_level = True
        for line in lines:
            x = KrakenData(line)
            if x.name == args.taxonomy:  # Check if we've hit desired taxonomy. If yes, set our write output flag
                # and go to next loop iteration.
                write_output = True
                tax_level = taxonomy_order.index(x.tax_level)
                check_tax_level = False
                # continue
            # If we've started writing, check that we haven't escaped our tax level. If we have, stop.
            if tax_level is not None and check_tax_level:
                if taxonomy_order.index(x.tax_level) <= tax_level:
                    break

            # If we're at our desired level, make sure we write our output.
            if x.tax_level.upper() == args.level.upper() and write_output and x.num_reads > 0:
                if args.full_taxonomy:
                    full_taxonomy = taxid_to_lineage_string(x.ncbi_tax_id)
                    output_dict[full_taxonomy] = str(x.num_reads)
                else:
                    output_dict[x.name] = str(x.num_reads)
            # Sometimes we get groups or subphylums or superorders or something that cause all sorts of problems.
            # ETE3 classifies them as no rank, but they still get reads assigned to them that end up being missed.
            # If this is the case, we need to check their taxonomic level to see if they're above our tax level,
            # and if so output their number of directly assigned reads.
            elif x.tax_level == 'Other' and write_output:
                if args.full_taxonomy:
                    full_taxonomy = taxid_to_lineage_string(x.ncbi_tax_id)
                    # Only want to write output if the level isn't beyond specified tax level.
                    # Since we can't do it directly, just check that the identifiers for levels beyond desired
                    # level aren't there (i.e. no _g denoting genus found when we only want to go to family level).
                    good_to_write = True
                    for i in range(taxonomy_order.index(args.level), len(taxonomy_order) - 1):
                        if taxonomy_order[i][0].lower() + '_' in full_taxonomy:
                            good_to_write = False
                    if good_to_write:
                        if ',' in x.name:
                            x.name = x.name.replace(',', ':')
                        # Final check: sometimes these taxonomic inbetween groups are truly an extra level (so
                        # between a level between family and genus, for example), and so for a genus level report we
                        # only want to report number of reads that map to this level and no further(num_reads_direct).
                        # Other times, these essentially take the place of a level (so it may not be called a genus, but
                        # is the only level between family and species, and therefore everything mapping to it or
                        # anything below it should be reported(num_reads).
                        # Need to figure out which of these two scenarios is happening.
                        # Use descendant nodes to figure out what level we're equivalent to. If equivalent to desired
                        # level, we're taking it's place, anything higher is just another level.
                        inferred_tax_level = determine_unassigned_rank(x.ncbi_tax_id)
                        if inferred_tax_level.capitalize() == args.level and x.num_reads != 0:
                            output_dict[full_taxonomy + ';' + x.name] = str(x.num_reads)
                        elif x.num_reads_direct != 0:
                            output_dict[full_taxonomy + ';' + x.name] = str(x.num_reads_direct)
            # Other case is that we're at a tax level higher than desired. Here, report unassigned reads.
            elif write_output and taxonomy_order.index(args.level) > taxonomy_order.index(x.tax_level) \
                    and x.num_reads_direct > 0:
                if args.full_taxonomy:
                    full_taxonomy = taxid_to_lineage_string(x.ncbi_tax_id)
                    output_dict[full_taxonomy + '_unassigned'] = str(x.num_reads_direct)
                else:
                    output_dict[x.name + '_unassigned'] = str(x.num_reads_direct)

            check_tax_level = True
        output_list_of_dicts.append(output_dict)

    # The above gets us to a point where we know what each sample has - now need to write it all to a nice output file.
    # Make sure to account for data that's not in all files!
    # To do so, do some brute force all-vs-all checking.
    for sample_one in output_list_of_dicts:
        for sample_two in output_list_of_dicts:
            for item in sample_one:
                if item not in sample_two:
                    sample_two[item] = 0

    output_header = 'Sample,'
    for item in output_list_of_dicts[0]:
        if item != 'Sample':
            output_header += item + ','
    output_header = output_header[:-1]
    with open(args.output_file, 'w') as outfile:
        outfile.write(output_header + '\n')
        for sample in output_list_of_dicts:
            output_line = ''
            for item in output_header.split(','):
                output_line += str(sample[item]) + ','
            output_line = output_line[:-1]
            outfile.write(output_line + '\n')

