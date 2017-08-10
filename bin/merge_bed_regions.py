#!/usr/bin/env python

'''
Greg McInes
Stanford University
gmcinnes@stanford.edu
'''

'''
This script takes in a bed file and collapses all the regions sharing the same tag into a single region.  
bedtools only joins overlapping regions, not based on gene name, as far as I can tell.

For example:

Input:
chr1 1000 2000 GENE_A
chr1 3000 4000 GENE_A

Output:
chr1 1000 4000 GENE_A
'''

import argparse

class MergeBed(object):
    def __init__(self, file, buffer=0, debug=False):
        self.debug = debug
        self.buffer = buffer
        self.run(file)

    def run(self, file):
        genes = {}
        gene_order = []
        with open(file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split()
                if fields[3] in genes:
                    genes[fields[3]].update_coordinates(int(fields[1]), int(fields[2]))
                else:
                    gene = Gene(fields[0], fields[1], fields[2], fields[3])
                    genes[gene.name] = gene
                    gene_order.append(gene.name)

        for gene in gene_order:
            genes[gene].print_bed(buffer=self.buffer)


class Gene(object):
    def __init__(self, chr, start, stop, name):
        self.chr = chr
        self.start = int(start)
        self.stop = int(stop)
        self.name = name

    def update_coordinates(self, new_start, new_stop):
        if new_start < self.start:
            self.start = new_start
        if new_stop > self.stop:
            self.stop = new_stop

    def print_bed(self, buffer=0):
        print("%s\t%s\t%s\t%s" % (self.chr, max(self.start-buffer, 0), self.stop+buffer, self.name))



"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This is a script I wrote')
    parser.add_argument("bed_file", help="Bed file to merge genes")
    parser.add_argument("-b", "--buffer", type=int, help="Add a buffer upstream and downstream of the merged regions")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    MergeBed(options.bed_file, options.buffer, options.debug)

