'''
Remove controls in the microarray and only keep valid genes.
'''

import sys

def remove_control(in_file, out_file):

    in_fh = open(in_file,'r')
    out_fh = open(out_file, 'w')
    out_fh.write(in_fh.readline()) # write header
    for line in in_fh:
        toks = line.strip().split('\t')
        gid = toks.pop(0).split('_')[0]
        if gid == 'Pae' or gid == 'ig' or gid.startswith('AFFX'):
            continue #skip controls/etc
        else:
            out_fh.write(gid + '\t' + '\t'.join(toks) + '\n')   
    out_fh.close()
    in_fh.close()

remove_control(sys.argv[1],sys.argv[2])