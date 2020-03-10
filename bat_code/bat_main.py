#
# Python 2.7
#
# v. 1.0
#
# BAT
# Bfactor Alignment Tool
#
# Joe Bemister
# 07/17/2015
#
#
# last modified:
# 07/17/2015
#
# ################
# DESCRIPTION
# ################
#
# Program reads in a pdb file with bfactor values (originally intended for
# Proflex flexibility indeces) and aligns them based on a DALI alignment file.


import sys
import os

#importing own python files from same dir
import bat_optparse
import bat_read_align_file
import bat_print_align
import bat_read_pdb

#load menu
bat_optparse

#Writing to File at the same time (From Sebastian Raschka's HETHER)
#######################################################

if '.txt' in bat_optparse.salign:
    output_file = bat_optparse.salign[:-4] + '_bat.csv'
else:
    output_file = bat_optparse.salign + '_bat.csv'

print('Results written to {}'.format(output_file))

#######################################################






def main():
    # opens, reads, and closes the input files. returns dictionary of c alphas.
    d_target = bat_read_pdb.read(bat_optparse.starget)
    d_query = bat_read_pdb.read(bat_optparse.squery)

    # opens, reads, and closes the input file. returns a dictionary query alignment
    # and a pdb code for the target and query.
    d_query_align, target, query = bat_read_align_file.read(bat_optparse.salign)

    # prints the alignment to the output file
    bat_print_align.print_align(d_target, d_query, d_query_align, target, query, output_file)

main()
