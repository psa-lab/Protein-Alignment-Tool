# Module that reads a pdb file
# for BAT
#
#
# Joe Bemister
#
# 08/03/2015
#
#
# modified:
#
#
#

def read(sFilename):
    '''Opens a pdb file and returns a dictionary of chain ID, residue ID,
    and the bfactor associated with a residue number.'''
    try:
        objFile = open(sFilename, 'r')
    except IOError:
        print("\n\n    ERROR: could not open {}\n\n".format(sFilename))
        sys.exit()


    d_calpha = {}

    for line in objFile:
        #would include specified chain here
        if line[:4] == 'ATOM' and line[13:15] == 'CA':
            d_calpha[int(line[22:26])] = [line[21], line[17:20], float(line[60:66])]
        else:
            continue

    objFile.close()

    return d_calpha
