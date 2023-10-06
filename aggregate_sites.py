import argparse
from collections import defaultdict
import gzip
import bz2
import os
import warnings

parser = argparse.ArgumentParser(description='')
parser.add_argument('--input_dir','-i', action='store',type=str,required=True,
                    help='Input directory with all the files')
parser.add_argument('--chromosome','-c', action='store',type=str,required=True,help='')
args = parser.parse_args()

if os.path.isdir(os.path.abspath(args.input_dir)):
    files= []
    for (dirpath, dirnames, filenames) in os.walk(os.path.abspath(args.input_dir)):
        for filename in filenames:
            files.append(dirpath+'/'+filename)
else:    
    raise Exception("Input must be a directory")

def openfile(file):
    '''
    Opens a file
    '''
    if file.endswith('.gz'):
        opened_file = gzip.open(file,'rt')
    elif file.endswith('bz') or file.endswith('bz2'):
        opened_file = bz2.open(file,'rt')
    else:
        opened_file = open(file,'rt')
    return opened_file

final_dict_meth= defaultdict(int)
final_dict_all= defaultdict(int)

for file in files:
    warnings.warn("Processing {} file".format(file))
    with openfile(file) as inp:
        for line in inp:
            line=line.rstrip().split('\t')
            if line[0] == args.chromosome and int(line[9]) > 0:
                final_dict_meth[(line[0],line[1],line[2],line[5])] += int(round(int(line[9])*(float(line[10])/100)))
                final_dict_all[(line[0],line[1],line[2],line[5])] += int(line[9])
    
for key, val in final_dict_all.items():
    meth= final_dict_meth[key]
    print('\t'.join(key)+'\t'+str(val)+'\t'+str(meth)+'\t'+str(meth/val))


