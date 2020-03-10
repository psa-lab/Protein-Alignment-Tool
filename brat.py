#
# Python 2.7 or 3.3+
#
# v. 1.0.1
#
# BRAT
# Binding Residue Alignment Tool
#
# Joe Bemister
# 2016-06-27
#
#
# last modified:
# Joe Bemister
# 2016-08-30
#
# ################
# DESCRIPTION
# ################
#
# Program reads in a csv or pdb file with important query residues
# and/or takes ligand info from the user and determines binding
# residues based on distance, aligns them based on a DALI alignment
# file, and emphasizes important residues from the query (and
# residues from the target aligned with the query residues) within a
# given radius of the target ligand.


import sys
import os
import time
import optparse

#importing own python files from same dir
#import brat_optparse
#import brat_print_output
#import brat_read
#import brat_pdb_functions


######################################################
#from brat_optparse

#Initialize
parser = optparse.OptionParser()

#Feed options

parser.add_option("-q", help="Location of the query pdb file", dest="query", action="store", metavar="<file>")
parser.add_option("-t", help="Location of the target pdb file to align to query", dest="target", action="store", metavar="<file>")
parser.add_option("-a", help="Location of the file with the desired alignment", dest="align", action="store", metavar="<file>")
parser.add_option("-b", help="Location of the binding residue file (csv/txt or pdb)", dest="res_file", action="store", metavar="<file>")
parser.add_option("-l", help="Case sensitive ligand code, chain, and number in file. Format as: <lig_code>,<chain>,<number>" + \
"   If no chain is presented for the ligand then format as: <lig_code>,,<number>", dest="lig_ID", action="store")
parser.add_option("-s", help="Span of residue numbers in the query file for the ligand. Format as: <chain>,<beg_number>-<end_number>,...", dest="lig_span", action="store")
parser.add_option("-r", help="The radius around the ligand in which all residues count as binding residues (default = 5)", dest='radius', action="store")


#Tell optparse to parse arguments
#opts contains all values reveived via command line
#every other arguments that are not recognized by our parser

(opts, args) = parser.parse_args()


#checking for mandatory options
if opts.query is None:
    print("\n\nERROR: query file is missing\n\n")
    parser.print_help()
    print("\n\n")
    sys.exit()
elif opts.target is None:
    print("\n\nERROR: target file is missing\n\n")
    parser.print_help()
    print("\n\n")
    sys.exit()
elif opts.align is None:
    print("\n\nERROR: alignment file is missing\n\n")
    parser.print_help()
    print("\n\n")
    sys.exit()
elif opts.lig_ID is None and opts.lig_span is None and opts.res_file is None:
    print("\n\nERROR: no input for binding residues\n\n")
    parser.print_help()
    print("\n\n")
    sys.exit()
elif opts.lig_ID and not opts.radius or opts.lig_span and not opts.radius:
    print("\n\nUser has requested BRAT's radius identification of binding residues"
        "\nbut did not specify a radius. The default 4.5 Angstrom radius will"
        "\nbe used.")
    opts.radius = 4.5


#access variables by opts.<dest>
squery = opts.query
starget = opts.target
salign = opts.align
sres_file = opts.res_file
slig_ID = opts.lig_ID
slig_span = opts.lig_span
#checking for usable radius input
try:
    iradius = float(opts.radius)
    if iradius <= 0:
        print("\n\nERROR: the radius must be a positive integer.\n\n")
        sys.exit()
except ValueError:
    print("\n\nERROR: the radius must be a positive integer.\n\n")
    sys.exit()
#if radius is not needed and not specified
except TypeError:
    iradius = False


######################################################
#from brat_read

def missing_query_align(d_target_align):
    '''takes the dictionary of an alignment the target residues to the query residues. This function creates
    a list of lists of target residues that are present outside of the alignment that cause gaps in the query.
    '''
    prev, l_q_align, l_begin = 0, [], []
    for val in sorted(d_target_align.values()):
        if val - prev > 1:
            tmp_list = []
            l_begin.append(prev)
            for i in range(prev, val):
                tmp_list.append(i)
            l_q_align.append(tmp_list)
        prev = val
    l_q_align.append(l_begin)
    return l_q_align

def read_pdb_CA(sFilename, incl_chain = False):
    '''Opens a pdb file and returns a list of the residues present in
    the pdb.'''
    try:
        objFile = open(sFilename, 'r')
    except IOError:
        print("\n\n    ERROR: could not open {}\n\n".format(sFilename))
        sys.exit()


    l_CA_res, l_res, chain, prev_mut = [], [], ' ', []
    count = 0
    # dictionary of possible mutated side chains
    d_special_residues = {'TPO':['THR', 'PHOSPHOTHREONINE'], 'PTR':['TYR', 'PHOSPHOTYROSINE'], 'SEP':['SER', 'PHOSPHOSERINE'],
        '2ML':['LEU', '2-METHYLLEUCINE'], 'M3L':['LYS', 'N-TRIMETHYLLYSINE'], 'MLY':['LYS', 'N-DIMETHYL-LYSINE'],
        'ACE':['XXX', 'ACETYL GROUP'], 'LMQ':['GLN', '(3S)-3-METHYL-L-GLUTAMINE'], 'MEQ':['GLN', 'N5-METHYLGLUTAMINE'],
        'NCY':['CYS', 'N-METHYLCYSTEINE'], 'HIC':['HIS', 'METHYL-HISTIDINE'], 'AGM':['ARG', '5-METHYLARGININE'],
        '2MR':['ARG', 'N3,N4-DIMETHYLARGININE'], 'SME':['MET', 'METHIONINE SULFOXIDE']}

    for line in objFile:
        # to handle alternative side chain locations
        if line[:4] == 'ATOM' and line[16] == 'B' or line[:4] == 'ATOM' and line[16] == 'C' or line[:4] == 'ATOM' and line[16] == '2':
            continue
        # handle mutated residues. This assumes that any entry of these is a part of the protein and not a ligand
        # even if it is a HETATM entry.
        if line[:4] == 'ATOM' and line[17:20] in d_special_residues or line[:6] == 'HETATM' and line[17:20] in d_special_residues:
            if line[17:20] not in prev_mut:
                print('\nIn file: {}\nBRAT found the modified residue {} in chain {} and used {} for the alignment.'.format(
                    sFilename, line[17:20], line[21], d_special_residues[line[17:20]][0]))
                prev_mut.append(line[17:20])
            line = 'ATOM  ' + line[6:17] + d_special_residues[line[17:20]][0] + line[20:]
        # skip lines that aren't ATOM or HETATM entries
        if not incl_chain and line[:4] != 'ATOM' or incl_chain and line[:4] != 'ATOM' and line[:6] != 'HETATM':
            continue
        # build up list of residues
        if line[21] != chain:
            if l_res != []:
                l_CA_res.append([chain, l_res])
                l_res = []
            chain = line[21]
        if line[13:15] == 'CA':
            if incl_chain:
                l_res.append(line[17:28].replace(' ', ''))
            else:
                l_res.append((line[17:20] + line[22:28]).replace(' ', ''))
        else:
            continue

    l_CA_res.append([chain, l_res])

    objFile.close()

    return l_CA_res

def read_pdb_coords(sFilename, l_lig_ID, range = False):
    '''Opens a pdb file and returns a list of the residues and a
    list of the heavy atoms in the ligand.'''
    try:
        objFile = open(sFilename, 'r')
    except IOError:
        print("\n\n    ERROR: could not open {}\n\n".format(sFilename))
        sys.exit()

    l_atom_res, l_lig_atom, prev_atom = [], [], ' '

    # list of ligand heavy atoms if only want to take heavy atoms in account for getting radius
    #l_heavy_atoms = ["O", "C", "P", "S", "N", "CA", "MG", "MN", "Fe", "CU", "ZN", "NA", "K"]

    # dictionary of possible mutated side chains
    d_special_residues = {'TPO':['THR', 'PHOSPHOTHREONINE'], 'PTR':['TYR', 'PHOSPHOTYROSINE'], 'SEP':['SER', 'PHOSPHOSERINE'],
        '2ML':['LEU', '2-METHYLLEUCINE'], 'M3L':['LYS', 'N-TRIMETHYLLYSINE'], 'MLY':['LYS', 'N-DIMETHYL-LYSINE'],
        'ACE':['XXX', 'ACETYL GROUP'], 'LMQ':['GLN', '(3S)-3-METHYL-L-GLUTAMINE'], 'MEQ':['GLN', 'N5-METHYLGLUTAMINE'],
        'NCY':['CYS', 'N-METHYLCYSTEINE'], 'HIC':['HIS', 'METHYL-HISTIDINE'], 'AGM':['ARG', '5-METHYLARGININE'],
        '2MR':['ARG', 'N3,N4-DIMETHYLARGININE'], 'SME':['MET', 'METHIONINE SULFOXIDE']}

    for line in objFile:
        # to handle alternate side chain locations would fail if 4 possible alternate side chains
        if line[:4] == 'ATOM' and line[16] == 'B' or line[:4] == 'ATOM' and line[16] == 'C' or line[:4] == 'ATOM' and line[16] == '2':
            continue
        # handle mutated residues. This assumes that any entry of these is apart of the protein and not a ligand
        # even if it is a HETATM entry.
        if line[:4] == 'ATOM' and line[17:20] in d_special_residues or line[:6] == 'HETATM' and line[17:20] in d_special_residues:
            line = line[:17] + d_special_residues[line[17:20]][0] + line[20:]
        # gathering the atoms based on which type of ligand (single molecule or span) was input by user
        if range:
            # this, for now, excludes residues with a number in the span even if different chain
            if line[:4] == 'ATOM' and int(line[22:26]) not in l_lig_ID[1]:
                l_atom_res.append([[float(line[31:38]), float(line[39:46]), float(line[47:54])], line[17:26].replace(' ','')])
            elif line[:4] == 'ATOM' and line[21] == l_lig_ID[0] and int(line[22:26]) in l_lig_ID[1] or \
                line[:6] == 'HETATM' and line[21] == l_lig_ID[0] and int(line[22:26]) in l_lig_ID[1]:
                l_lig_atom.append([float(line[31:38]), float(line[39:46]), float(line[47:54])])
        else:
            if line[:4] == 'ATOM':
                #could include more info in dict if needed. such as chain ID or bfactor
                l_atom_res.append([[float(line[31:38]), float(line[39:46]), float(line[47:54])], line[17:26].replace(' ','')])
            elif line[:6] == 'HETATM' and line[17:26].replace(' ', '') == ''.join(l_lig_ID):
            #elif line[:6] == 'HETATM' and line[17:21].replace(' ','') == l_lig_ID[0]:
                #print(line[17:26] + ' ' + line[17:26].replace(' ', '') + ' ' + ''.join(l_lig_ID))
                l_lig_atom.append([float(line[31:38]), float(line[39:46]), float(line[47:54])])

    objFile.close()

    return l_atom_res, l_lig_atom

def read_bind_res_csv(sFilename, combined_list):
    '''Opens a csv file and returns two lists of the items in the csv.
    One is to note that they came from the file. The other is meant to
    be appended and contain residues from another source.
    '''

    try:
        objFile = open(sFilename, 'r')
    except IOError:
        print("\n\n    ERROR: could not open {}\n\n".format(sFilename))
        sys.exit()


    user_bind_res, tmp_list, chain = [], [], ''

    for line in objFile:
        if line == '\n':
            continue
        for item in line.split(','):
            if item[-1] == '\n':
                item = item[:-1]
            # deal with a faulty binding residue entry
            if not item[:3].isalpha() or item[3].isdigit() or len(item) > 4 and not item[4:].isdigit():
                #print("\nbinding residue {} is not of proper format.".format(item))
                #print("It will not be considered a binding residue for the alignment.\n")
                #continue
                print("\nbinding residue {} is not of proper format.".format(item))
                print("This can cause an issue with the alignment. Please fix this and resubmit.\n")
                sys.exit()
            item = item.strip()

            # deal with binding residues without a chain. this function assumes there is
            # there is only one chain in the entire file when no chain is given.
            # I have decided to call a binding residue following format "ASP200" to be
            # considered a faulty entry. If changed then decomment below.
            #if item[3].isdigit():
            #    chain = ' '
            #    tmp_list.append(item[:3] + ' ' + item[3:])
            #    combined_list.append(item[:3] + ' ' + item[3:])
            #    continue

            # deal with binding residues without a chain. this function assumes there is
            # only one chain in the entire file when no chain is given.
            if item[3] == ' ':
                chain = ' '
                #tmp_list.append(item)
                #combined_list.append(item[:3] + item[4:])
                #continue

            # if a residue is not of the same chain. adds tmp_list to user_bind_res
            # and restarts tmp_list to only contain residues of the next chain.
            if chain != ' ' and item[3] != chain:
                if tmp_list == []:
                    chain = item[3]
                else:
                    user_bind_res.append([chain, tmp_list])
                    tmp_list, chain = [], item[3]
            tmp_list.append(item)
            combined_list.append(item)

    user_bind_res.append([chain, tmp_list])
    #print(combined_list)

    objFile.close()

    return user_bind_res, combined_list

def read_bind_res_pdb(sFilename, combined_list):
    '''Opens a pdb file and returns two lists of the items in the csv.
    One is to note that they came from the file. The other is meant to
    be appended and contain residues from another source.
    '''

    res_list = read_pdb_CA(sFilename, incl_chain=True)
    for item in res_list:
        for bind_res in item[1]:
            #print(bind_res)
            combined_list.append(bind_res)

    return res_list, combined_list

def read_DALI_alignment(sFilename):
    '''Opens a DALI alignment file and returns an alignment dictionary for the query
    as well as a pdb code for the query and the target protein. The dictionary looks
    as follows:

    {<query_listpos>:<target_listpos>, <query_listpos>:<target_listpos>, ...}     '''

    try:
        objFile = open(sFilename, 'r')
    except IOError:
        print("\n\n    ERROR: could not open {}\n\n".format(sFilename))
        sys.exit()

    l_target_align, l_query_align = [], []

    # read 1st two lines with query and target titles
    line = objFile.readline()
    query = line[-5:-1]
    line = objFile.readline()
    target = line[-5:-1]

    # get to the actual alignment
    while ("# Structural equivalences" not in line):
        line = objFile.readline()

    # get information from the actual alignment
    chain_a, chain_b, d_target_align = ' ', ' ', {}
    for line in objFile:
        if '<=>' not in line:
            continue
        if line[11] != chain_a or line[18] != chain_b:
            if d_target_align != {}:
                l_target_align.append([chain_a, chain_b, d_target_align])
                # responsible for '-'s in query output
                l_query_align.append([chain_a, chain_b, missing_query_align(d_target_align)])

                d_target_align = {}
            chain_a, chain_b = line[11], line[18]
        index1, index2 = line.find('<=>'), line.rfind('<=>')
        for i in range(int(line[index1 - 11: index1 - 7]) - 1, int(line[index1 - 5: index1 - 1])):
            d_target_align[i] = i + int(line[index1 + 4: index1 + 8]) - int(line[index1 - 11: index1 - 7])

    l_target_align.append([chain_a, chain_b, d_target_align])
    # responsible for '-'s in query output
    l_query_align.append([chain_a, chain_b, missing_query_align(d_target_align)])


    objFile.close()
    return l_target_align, l_query_align, query, target


######################################################
#from brat_pdb_functions

def grab_radius(l_query_res, l_atom_coord, l_res_win_rad, radius):
    '''
    finds atoms within a radius of the atom coordinates.
    Requires a list of the query residues as such...:

    [[[x_coord, y_coord, z_coord], <res_ID><chain><res_num>], ...]

    a list as...:

    [x_coord, y_coord, z_coord]

    a previously created list to return formatted as such...:

    [<res_ID><res_num>, <res_ID><res_num>, ...]

    and a float radius.
    '''

    #l_res_win_rad_new = l_res_win_rad

    for atom in l_query_res:
        fDist =  (sum([(l_atom_coord[i]-atom[0][i])**2 for i in range(3)]))**0.5
        round(fDist,11)
        if fDist <= radius and atom[1] not in l_res_win_rad:
            l_res_win_rad.append(atom[1])

    return l_res_win_rad

def get_bind_res(file, lig_ID, radius, total_bind_res, b_main_span = False):
    '''
    Grabs the residues of a protein surrounding the ligand. returns a list of these
    residues as such:

    [<res_ID><res_num>, <res_ID><res_num>, ...]

    requires str of an unopened file, str of ID of ligand in file, and list to add
    residues to.
    '''

    ### CHECKING INPUT
    if not os.path.isfile(file):
        print("ERROR: {} does not exists.\n".format(file))
        raise FileError
    if file.split(".")[-1].lower() != "pdb":
        print("ERROR: {} is not a pdb file\n".format(file))
        raise FileError

    brat_bind_res, tmp_list, chain_brat_bind_res, chain = [], [], [], ''

    # gather binding residues regardless of chain
    # implementation of span or single ligand atom collection
    if b_main_span:
        l_query_res, l_heavy_atom = read_pdb_coords(file, lig_ID, range=True)
    else:
        l_query_res, l_heavy_atom = read_pdb_coords(file, lig_ID)
    for l_atom in l_heavy_atom:
        brat_bind_res = grab_radius(l_query_res, l_atom, brat_bind_res, radius)
    #print("brat_bind_res: {}".format(brat_bind_res))

    # create new list of binding residues sorted by chain
    for res in brat_bind_res:
        if res[3].isdigit() or chain == ' ':
            chain = ' '
            tmp_list.append(res[:3] + ' ' + res[3:])
            continue
        elif res[3] != chain:
            if tmp_list == []:
                chain = res[3]
            else:
                chain_brat_bind_res.append([chain, tmp_list])
                tmp_list, chain = [], res[3]
        tmp_list.append(res)
    chain_brat_bind_res.append([chain, tmp_list])
    #print("chain_bind_res: {}".format(chain_brat_bind_res))

    # add binding residues to total_bind_res
    for res in brat_bind_res:
        if res not in total_bind_res:
            total_bind_res.append(res)


    return chain_brat_bind_res, total_bind_res


######################################################
#from brat_print_output

def add_html_spaces(s_num):
    ''' takes a number (str) and adds html spaces using "&nbsp;" until len == len(1000).
        primarily left this one line function in to prevent print_align from being
        even harder to understand.
    '''
    # Don't leave spaces anymore so will probably have to adjust this to until len = len("1000")
    while len(s_num.replace('&nbsp;', ' ')) < len("1000"):
        s_num = '&nbsp;' + s_num
    return s_num
    #return s_num.replace(' ', '&nbsp;')

def print_html_header(query, target):
    ''' takes the target and query title from the alignment file and returns an html header
        including a title with the target and query title.
    '''
    header = '<!doctype html>\n<html>\n<head>\n<meta charset="utf-8">\n' + \
                '<meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes">\n' + \
                '<title>BRAT alignment of {} and {}</title>\n</head>\n'.format(query, target)
    return header

#def equal_title_length(query, target):
#    if target > query:
#        while(len(query) < len(target)):
#            query += ' '
#    else:
#        while(len(target) < len(query)):
#            target += ' '
#    return query, target


def print_align(l_query, l_target, d_target_align, l_query_align, query, target, query_chain, target_chain, bind_res, ofile, query_f, target_f):
    ''' Outputs the alignment to the output_file. While doing this it bolds the binding
        residues included in bind_res.

        Input:
            l_query (list): includes the residues included in the query pdb file
            l_target (list): includes residues included in the target pdb file
            d_target_align (dictionary): includes the indexes {<que_ind>:<tar_ind>,...} of the query and target lists
                that should be aligned.
            l_query_align (list): includes the indexes of the query where target residues are missing from the
                alignment as such:
                    [[<miss_ind_1a>, <miss_in_2>,...],[<miss_ind_1b>,...],...[<miss_ind_1a>, <miss_ind_1b>, ...]]
            query (string): the title of the query protein
            target (string): the title of the target protein
            bind_res (list): list of the residues defined as being involved in binding
            ofile (file open for writing): file opened in brat_main using "with open() as " so no need to close
            query_f (string): path of the query file
            target_f (string): path of the target file

        The output contains lines of 50 residues formatted as such:

        <query_ID>' '<starting_res_num><start_res_ID>...<end_res_ID><ending_res_num>
        <target_ID>' '<starting_res_num><start_res_equiv_ID>...<end_res_equiv_ID><ending_res_num>

        <query_ID>' '<starting_res_num><start_res_ID>...<end_res_ID><ending_res_num>
        <target_ID>' '<starting_res_num><start_res_equiv_ID>...<end_res_equiv_ID><ending_res_num>
'''

# for testing in input is correct
#    count = 0
#    print("l_query\n")
#    print(l_query)
#    for item in l_query:
#        if count == 50:
#            break
#        print(item)
#        count += 1
#    count = 0
#    print("l_target\n")
#    for item in l_target:
#        if count == 50:
#            break
#        print(item)
#        count += 1
#    count = 0
#    print("d_target_align\n")
#    for k,v in d_target_align.items():
#        if count == 50:
#            break
#        print("k: " + str(k) + " v: " + str(v))
#        count += 1
#    count = 0
#    print("l_query_align\n")
#    for item in l_query_align:
#        if count == 50:
#            break
#        print(item)
#        count += 1
#    count = 0
#    print("bind_res: {}".format(bind_res))



    d_res_ID = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', \
        'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', \
        'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'XXX':'X'}
    bind_res_dict = {'chain': query_chain}

    ofile.write(print_html_header(query, target))
    ofile.write('<body>\n<div style="font-family:courier">\n<p>\n<br>Alignment for {} chain {} and {} chain {}<br>\n'.format(query, query_chain, target, target_chain) + '{} = {}<br>\n{} = {}<br><br>\n'.format(query, query_f, target, target_f))

    # build output and write to file
    output_1, output_2, count, prev_res, prev_bind_res = query + ' ', target + ' ', 0, '', False

#    for res_q in sorted(l_query, key=lambda string: int(string[3:])):
#        ind = sorted(l_query, key=lambda string: int(string[3:])).index(res_q)
    for res_q in l_query:
        ind = l_query.index(res_q)
        key_error = False

# Block responsible for ending and starting lines if query residues is present in new line start position.
################################################################

        # starting new lines after 50 residues
        if count == 50:
            if prev_bind_res == False:
                #ofile.write(output_1 + ' ' + add_html_spaces(sorted(l_query, key=lambda string: int(string[3:]))[ind - 1][3:]) + '<br>\n')
                ofile.write(output_1 + ' ' + add_html_spaces(l_query[ind - 1][3:]) + '<br>\n')
                # might not accurately count if last residue in line has no match in query
                #ofile.write(output_2 + ' ' + add_html_spaces(l_target[d_target_align[min(d_target_align.keys())] - 1][3:])  + '<br>\n<br>\n')
                ofile.write(output_2 + ' ' + add_html_spaces(prev_res[3:])  + '<br>\n<br>\n')
                output_1, output_2, count = query + ' ', target + ' ', 0
            else:
                #ofile.write(output_1 + '</b> ' + add_html_spaces(sorted(l_query, key=lambda string: int(string[3:]))[ind - 1][3:]) + '<br>\n')
                ofile.write(output_1 + '</b> ' + add_html_spaces(l_query[ind - 1][3:]) + '<br>\n')
                # may not accurately count if last residue in line has no match in query
                ofile.write(output_2 + '</b> ' + add_html_spaces(prev_res[3:])  + '<br>\n<br>\n')
                output_1, output_2, count, prev_bind_res = query + ' ', target + ' ', 0, False

        # adding leading resnum when necessary.
        if count == 0:
            output_1 += add_html_spaces(res_q[3:]) + ' '
            try:
                output_2 += add_html_spaces(l_target[d_target_align[ind]][3:]) + ' '
            except KeyError:
                output_2 += add_html_spaces(l_target[d_target_align[min(d_target_align.keys())]][3:]) + ' '
                key_error = True

###################################################################################################
# End of block responsible for ending and starting lines if query residues is present in new line start position.


# The block that outputs the residues. It takes into account binding residues and missing residues.
###################################################################################################

    # The subblock that handles cases where the query residue is present in the position.
    ############################################################
        # if current residue is considered a binding residue without a binding residue directly before
        if res_q.replace(' ', '') in bind_res and prev_bind_res == False:
            prev_bind_res = True
            output_1 += '<b>' + d_res_ID[res_q[:3]]
            try:
                output_2 += '<b>' + d_res_ID[l_target[d_target_align[ind]][:3]]
                bind_res_dict[res_q[:3] + query_chain + res_q[3:]] = l_target[d_target_align[ind]][:3] + target_chain + l_target[d_target_align[ind]][3:]
                prev_res = l_target[d_target_align[ind]]
            # outputs a - if nothing in d_target_align with the key as well as if the residue is not in d_res_ID
            except KeyError:
                output_2 += '<b>-'
                bind_res_dict[res_q[:3] + query_chain + res_q[3:]] = '-'
                key_error = True
        # if current residue is not considered a binding residue but there is a binding residue directly before
        elif res_q.replace(' ', '') not in bind_res and prev_bind_res:
            prev_bind_res = False
            output_1 += '</b>' + d_res_ID[res_q[:3]]
            try:
                output_2 += '</b>' + d_res_ID[l_target[d_target_align[ind]][:3]]
                prev_res = l_target[d_target_align[ind]]
            except KeyError:
                output_2 += '</b>-'
                key_error = True
            except IndexError:
                print(l_target)
                print(d_target_align[ind])
        # if current residue is considered a binding residue and there is a binding residue directly before
        elif res_q.replace(' ','') in bind_res and prev_bind_res:
            output_1 += d_res_ID[res_q[:3]]
            try:
                output_2 += d_res_ID[l_target[d_target_align[ind]][:3]]
                bind_res_dict[res_q[:3] + query_chain + res_q[3:]] = l_target[d_target_align[ind]][:3] + target_chain + l_target[d_target_align[ind]][3:]
                prev_res = l_target[d_target_align[ind]]
            except KeyError:
                output_2 += '-'
                key_error = True
        # if current residue is not considered a binding residue and there is not a binding residue directly before
        else:
            output_1 += d_res_ID[res_q[:3]]
            try:
                output_2 += d_res_ID[l_target[d_target_align[ind]][:3]]
                prev_res = l_target[d_target_align[ind]]
            except KeyError:
                output_2 += '-'
                key_error = True

    ############################################################
    # End of the subblock that handles cases where the query residue is present in the position.


    # The subblock that handles cases where the query residue is not present in the position.
    ############################################################

        # in the event ind not in d_target_align.
        try:
            if d_target_align[ind] in l_query_align[-1]:
                count += 1
                for l in l_query_align[:-1]:
                    if d_target_align[ind] in l:
                        for l_ in l[1:]:
                            # end and start lines if the query residue is not in present position.
                            # not yet tested in all situations
                            if count == 50:
                                #ofile.write(output_1 + ' ' + add_html_spaces(sorted(l_query, key=lambda string: int(string[3:]))[ind - 1][3:]) + '<br>\n')
                                ofile.write(output_1 + ' ' + add_html_spaces(l_query[ind - 1][3:]) + '<br>\n')
                                # might not accurately count if last residue in line has no match in query
                                ofile.write(output_2 + ' ' + add_html_spaces(prev_res[3:])  + '<br>\n<br>\n')
                                output_1, output_2, count = query + ' ', target + ' ', 0
                                #output_1 += add_html_spaces(sorted(l_query, key=lambda string: int(string[3:]))[ind + 1][3:]) + ' '
                                output_1 += add_html_spaces(l_query[ind + 1][3:]) + ' '
                                output_2 += add_html_spaces(l_target[l_][3:]) + ' '

                            output_1 += '-'
                            output_2 += d_res_ID[l_target[l_][:3]]
                            prev_res = l_target[l_]
                            count += 1
                        count -= 1
        # catching the event that ind not in d_target_align.
        except KeyError:
            key_error = True

    ############################################################
    # End of the subblock that handles cases where the query residue is not present in the position.

###################################################################################################
# End of the block that outputs the residues.


        # deleting the key,value just used in order to allow for minimum to be used on dictionary
        # for a new line of output that doesn't have a target residue present in the leading position.
        if key_error == False:
            del d_target_align[ind]

        count += 1

    # outputing the final line since the last line will never be written from the for res_q loop above.
    # this is because the last iteration through cannot possibly have a count of 50.
    ofile.write(output_1 + ' ' + add_html_spaces(res_q[3:]) + '<br>\n')
    # might not accurately count if last residue in line has no match in query
    ofile.write(output_2 + ' ' + add_html_spaces(prev_res[3:])  + '<br></p>')

    return bind_res_dict


def print_res(res_list, l_targ_bind_res, ofile, user = False):
    '''Outputs the residue list to the ALREADY OPEN brat_res.csv output file. Takes the form of:

    BRAT_input: <residue_ID><chainID><residue_number>, <residue_ID><chainID><residue_number>, ...

    or

    USER_input: <residue_ID><chainID><residue_number>, <residue_ID><chainID><residue_number>, ...

    depending on if the input was provided by the user or by BRAT.          '''
    #print(res_list)
    if user:
        out1 = "USER_bind_res, "
    else:
        out1 = "BRAT_bind_res, "
    out2 = "corresponding_target_res, "

    #print("sorted res_list: {}".format(sorted(res_list)))
    #print("l_targ_bind_res: {}".format(l_targ_bind_res))

    for chain in sorted(res_list):
        for corr_chain in l_targ_bind_res:
            if corr_chain['chain'] != chain[0]:
                continue
            for res in sorted(chain[1], key=lambda string: int(string[4:])):
                out1 += res + ','
                out2 += corr_chain[res] + ','

    ofile.write(out1[:-1] + '\n')
    ofile.write(out2[:-1] + '\n\n')


######################################################
#from brat_main

#load menu
#brat_optparse

# establishing output files
if '.txt' in salign:
    output_file1 = salign[:-4] + '_brat_res.csv'
    output_file2 = salign[:-4] + '_brat_aligned.html'
else:
    output_file1 = salign + '_brat_res.csv'
    output_file2 = salign + '_brat_aligned.html'




def main():
    # opens, reads, and closes the input file. returns a dictionary query alignment
    # and a pdb code for the query and target.
    l_target_align, l_query_align, query, target = read_DALI_alignment(salign)
    #print(d_target_align)

    # opens, reads, and closes the input files. returns dictionary of c alphas.
    l_query = read_pdb_CA(squery)
    l_target = read_pdb_CA(starget)

    # The next few lines determine the binding residues via user input or BRAT
    # input. These residues are placed into an individual list based on input
    # method. Another list total_bind_res is the union of the two individual
    # lists.
    total_bind_res = []

    # Identifies (if requested by -l) binding residues and adds user inputed binding
    # residues and/or brat residues to the output file.
    if sres_file:
        index = sres_file.rfind('.')
        if sres_file[index:] == '.csv' or sres_file[index:] == '.txt':
            user_bind_res, total_bind_res = read_bind_res_csv(sres_file, total_bind_res)
            #print("user_bind_res: {}".format(user_bind_res))
        elif sres_file[index:] == '.pdb':
            user_bind_res, total_bind_res = read_bind_res_pdb(sres_file, total_bind_res)
            #print("user_bind_res: {}".format(user_bind_res))
        else:
            print('\n\n Binding residue file is not of proper format.\n')
            sys.exit()
    if slig_span:
        index = slig_span.find(',')
        index2 = slig_span.find('-')
        lig_ID = [slig_span[:index],[i for i in range(int(slig_span[index+1:index2]),int(slig_span[index2+1:]) + 1)]]
        brat_bind_res, total_bind_res = get_bind_res(squery, lig_ID, iradius, total_bind_res, b_main_span = True)
        #print("brat_bind_res: {}".format(brat_bind_res))
    if slig_ID:
        lig_ID = slig_ID.split(',')
        brat_bind_res, total_bind_res = get_bind_res(squery, lig_ID, iradius, total_bind_res)
        #print("brat_bind_res: {}".format(brat_bind_res))


    # creates a list containing a list for the binding residues of each chain
    l_bind_res, tmp_list, chain = [], [], ''
    #print("sorted total_bind_res: {}".format(sorted(total_bind_res, key=lambda string: string[3:])))
    for item in sorted(total_bind_res, key=lambda string: string[3:]):
        #print("res: {}, digit: {}".format(item, item[3].isdigit()))
        if item[3].isdigit():
            chain = ' '
            tmp_list.append(item)
            continue
        if chain != ' ' and item[3] != chain:
            if tmp_list == []:
                chain = item[3]
            else:
                l_bind_res.append([chain, sorted(tmp_list, key=lambda string: int(string[3:]))])
                tmp_list, chain = [], item[3]
        tmp_list.append(item[:3] + item[4:])
        #print("tmp_list: {}".format(tmp_list))
    l_bind_res.append([chain, sorted(tmp_list, key=lambda string: int(string[3:]))])


    #checking all values before output
    #print("l_target_align:\n{}\n".format(l_target_align))
    #print("l_query_align:\n{}\n".format(l_query_align))
    #print("l_query:\n{}\n".format(l_query))
    #print("l_target:\n{}\n".format(l_target))
    #print("l_bind_res:\n{}\n".format(l_bind_res))

    l_targ_bind_res = []

    # takes each chain alignment and finds the corresponding query align (for '-'s), query residues, target residues, and binding residues
    # then outputs them to the file
    with open(output_file2, 'w') as ofile2:
        for chain_targ_align in l_target_align:
            #print("chain_targ_align:\n{}\n".format(chain_targ_align))
            for chain_quer_align in l_query_align:
                if chain_targ_align[0] == chain_quer_align[0] and chain_targ_align[1] == chain_quer_align[1]:
                    #print("chain_quer_align:\n{}\n".format(chain_quer_align))
                    for chain_quer_res in l_query:
                        if chain_targ_align[0] == chain_quer_res[0]:
                            #print("chain_quer_res:\n{}\n".format(chain_quer_res))
                            for chain_targ_res in l_target:
                                if chain_targ_align[1] == chain_targ_res[0]:
                                    #print("chain_targ_res:\n{}\n".format(chain_targ_res))
                                    for chain_bind_res in l_bind_res:
                                        if chain_targ_align[0] == chain_bind_res[0]:
                                            #print("chain_bind_res:\n{}\n".format(chain_bind_res))
                                            targ_bind_res = print_align(chain_quer_res[1], chain_targ_res[1], chain_targ_align[2], chain_quer_align[2], query, target, chain_targ_align[0], chain_targ_align[1], chain_bind_res[1], ofile2, squery, starget)
                                            l_targ_bind_res.append(targ_bind_res)
                                            break
                                    if chain_targ_align[0] != chain_bind_res[0]:
                                        targ_bind_res = print_align(chain_quer_res[1], chain_targ_res[1], chain_targ_align[2], chain_quer_align[2], query, target, chain_targ_align[0], chain_targ_align[1], [], ofile2, squery, starget)
                                    break
                            break
                    break

        # output BRAT version, date and time stamp, lab info, and closing statements for the html
        ofile2.write('\n\n<p><br>BRAT v1.0.1<br>\n' + time.strftime("%B %d %Y %H:%M:%S") + '<br>\n' + \
            'Protein Structural Analysis & Design Laboratory<br>\n' + \
            'www.kuhnlab.bmb.msu.edu</p>\n</div>\n\n</body>\n</html>')

    #print("l_targ_bind_res: {}".format(l_targ_bind_res))

    # printing binding residue csv output
    with open(output_file1, 'w') as ofile1:
        ofile1.write('Binding residues from the alignment of {} and {}\n\n'.format(query, target))
        #print("user_bind_res: \n{}\n".format(user_bind_res))
        #print("brat_bind_res: \n{}\n".format(brat_bind_res))
        #print("total_bind_res: \n{}\n".format(l_bind_res))
        #print("l_targ_bind_res: \n{}\n".format(l_targ_bind_res))
        if sres_file:
            print_res(user_bind_res, l_targ_bind_res, ofile1, user=True)
        if slig_span or slig_ID:
            print_res(brat_bind_res, l_targ_bind_res, ofile1)

        # printing self-documentation of binding residue identification to brat_res.csv output file
        if sres_file and slig_span:
            ofile1.write('\nUSER_bind_res was defined from {}\n'.format(salign) +
                'BRAT_bind_res was defined as the residues within {}'.format(iradius) +
                'Angstrom of residue span #:{} of chain:{} in {}\n'.format(slig_span[index + 1:], slig_span[:index], squery))
        elif sres_file and slig_ID:
            ofile1.write('\nUSER_bind_res was defined from {}\n'.format(salign) +
                'BRAT_bind_res was defined as the residues within {} '.format(iradius) +
                'Angstrom of ligand ID:{} #:{} of chain:{} in {}\n'.format(lig_ID[0], lig_ID[2], lig_ID[1], squery))
        elif sres_file:
            ofile1.write('\nUSER_bind_res was defined from {}\n'.format(salign))
        elif slig_span:
            ofile1.write('BRAT_bind_res was defined as the residues within {} '.format(iradius) +
                'Angstrom of residue span #:{} of chain:{} in {}\n'.format(slig_span[index + 1:], slig_span[:index], squery))
        else:
            ofile1.write('BRAT_bind_res was defined as the residues within {} '.format(iradius) +
                'Angstrom of ligand ID:{} #:{} of chain:{} in {}\n'.format(lig_ID[0], lig_ID[2], lig_ID[1], squery))

        # printing BRAT version and date and time stamp to brat_res.csv output file
        ofile1.write('\nBRAT v1.0.1\n' + time.strftime("%B %d %Y  %H:%M:%S") + '\nProtein Structural Analysis & Design Laboratory\nwww.kuhnlab.bmb.msu.edu')


main()

print('\n\nBinding residues written to {}'.format(output_file1))
print('Results written to {}'.format(output_file2))
