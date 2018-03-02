# Module that reads output alignment file from DALI
# for BAT
#
#
# Joe Bemister
#
# 07/30/2015
#
#
# modified:
#
#
#
import sys


def read(sFilename):
    '''Opens a DALI alignment file and returns an alignment dictionary for the query
    as well as a pdb code for the target and the query protein.'''
    try:
        objFile = open(sFilename, 'r')
    except IOError:
        print "\n\n    ERROR: could not open {}\n\n".format(sFilename)
        sys.exit()

    d_query_align = {}

    for line in objFile:
        if 'Position 1' in line:
            target = line[-5:-1]
        elif 'Position 2' in line:
            query = line[-5:-1]
        elif '<=>' in line:
            index1, index2 = line.find('<=>'), line.rfind('<=>')
            for i in range(int(line[index1 - 11: index1 - 7]), int(line[index1 - 5: index1 - 1]) + 1):
                d_query_align[i] = i + int(line[index2 + 8: index2 + 12]) - int(line[index1 - 11: index1 - 7])

    objFile.close()

    return d_query_align, target, query
