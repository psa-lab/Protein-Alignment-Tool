# Module that provides the input menu
# for BAT
#
#
# Joe Bemister
#
# 03/31/2015
#
#
# modified:
#
#
#


import optparse
import sys
import os

#Initialize
parser = optparse.OptionParser()
#log = os.path.join(os.getcwd(), 'log.txt')
#align = os.path.join(os.getcwd(), 'bat.csv')

#Feed options

parser.add_option("-t", help="Location of the target pdb file", dest="target", action="store", metavar="<file>")
parser.add_option("-q", help="Location of the pdb file to align to target", dest="query", action="store", metavar="<file>")
parser.add_option("-a", help="Location of the file with the desired alignment", dest="align", action="store", metavar="<file>")


#Tell optparse to parse arguments
#opts contains all values received via command line
#every other arguments that are not recognized by our parser

(opts, args) = parser.parse_args()


#checking for mandatory options (based on Sebastian Raschka's HETHER)
if opts.target is None:
    print("\n\nERROR: target file is missing\n\n")
    parser.print_help()
    print("\n\n")
    sys.exit()
elif opts.query is None:
    print("\n\nERROR: query file is missing\n\n")
    parser.print_help()
    print("\n\n")
    sys.exit()
elif opts.align is None:
    print("\n\nERROR: alignment file is missing\n\n")
    parser.print_help()
    print("\n\n")
    sys.exit()


#access variables by opts.<dest>
starget = opts.target
squery = opts.query
salign = opts.align
