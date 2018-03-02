# Module that outputs the alignment onto a csv
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


def print_align(d_target, d_query, d_query_align, target, query, ofile):
    '''Outputs the alignment to the output_file. Takes the form of:

    ,Chain ID and Residue Number,<chain ID><residue number>,...
    <pdb code>,Residue ID,<residue ID>,...
    ,Flexibility Index,<flexibility index>,...
                                                                    '''
    try:
        outFile = open(ofile, 'w')
    except IOError:
        print "\n\n    ERROR: could not open {}\n\n".format(ofile)
        sys.exit()

    output_1, output_2, output_3 = ',Chain ID and Residue Number', target + ',Residue ID', ',Flexibility Index'
    for key in sorted(d_target.keys()):
        output_1 += ',' + d_target[key][0] + str(key)
        output_2 += ',' + d_target[key][1]
        output_3 += ',' + str(d_target[key][2])
    #print('{}\n{}\n{}\n'.format(output_1, output_2, output_3))
    outFile.write('{}\n{}\n{}\n\n'.format(output_1, output_2, output_3))

    output_a, output_b, output_c = ',Chain ID and Residue Number', query + ',Residue ID', ',Flexibility Index'
    count = 0
    while count < len(output_1.split(',')) - 2:
        if count in d_query_align.keys():
            output_a += ',' + d_query[d_query_align[count]][0] + str(d_query_align[count])
            output_b += ',' + d_query[d_query_align[count]][1]
            output_c += ',' + str(d_query[d_query_align[count]][2])
        else:
            output_a += ',-'
            output_b += ',-'
            output_c += ',-'
        count += 1
    #print('{}\n{}\n{}\n'.format(output_a, output_b, output_c))
    outFile.write('{}\n{}\n{}\n\n'.format(output_a, output_b, output_c))
    outFile.close()

