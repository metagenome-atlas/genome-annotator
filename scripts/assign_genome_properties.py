#!/usr/bin/env python

import sys
sys.path.append("/Users/silas/Documents/GitHub/pygenprop")

from pygenprop.database_file_parser import parse_genome_property_file
from pygenprop.assignment_file_parser import parse_interproscan_file
from pygenprop.results import GenomePropertiesResults
import pandas as pd
import os


def asign_properties(interpro_assignments,genprop_flat_file):
    """
        Function to assign genomeproperties to one ore more interproscan annotation files.

        params:
            assignment_file_paths: one ore more interproscan annotation files
            genprop_flat_file: path to genomeProperties.txt flatfile

        returns a GenomePropertiesResults object.
        Takes sample name form file name.

    """


    if type(interpro_assignments) == str:
        interpro_assignments=[interpro_assignments]


    genprop_tree = parse_genome_property_file(open(genprop_flat_file))

    assignments = []

    for path in interpro_assignments:
        with open(path) as assignment_file:
            file_assignment_results = parse_interproscan_file(assignment_file)

            file_assignment_results.sample_name = os.path.splitext(os.path.split(path)[-1])[0]
            assignments.append(file_assignment_results)

    result = GenomePropertiesResults(*assignments, genome_properties_tree=genprop_tree)

    return result


if __name__ == "__main__":


    try:
        results= asign_properties(snakemake.input.interpro,snakemake.input.genprop_flat_file)

        results.property_results.to_csv(snakemake.output[0],sep='\t')

    except NameError:

        import argparse

        p = argparse.ArgumentParser()

        p.add_argument("-g","--genprop-flat-file",required=True)
        p.add_argument("-o","--output-file",required=True)
        p.add_argument("interpro_assignments",nargs='+')
        args = p.parse_args()

        results= asign_properties(args.interpro_assignments,args.genprop_flat_file)

        results.property_results.to_csv(args.output_file,sep='\t')
