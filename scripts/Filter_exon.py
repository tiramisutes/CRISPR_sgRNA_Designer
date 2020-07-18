#!/usr/bin/env python
# encoding: utf-8

"""
@version: v1.0
@author: Zhongping Xu
@license: Apache Licence
@contact: 1099808298@qq.com
@site: http://tiramisutes.github.io/
@software: python2.7
@time: 2020/7/17 18:46
@description：Filter fasta for run BE
"""

import os, sys, argparse
import regex as re
#import re
from Bio import SeqIO

#################################################################################
################################ General parameters #############################
#################################################################################
flank = 20


### pattern for sgRNA
pattern = '[ATCG]{21}GG|CC[ATCG]{21}'
patternF = '[ATCG]{21}GG'
patternR = 'CC[ATCG]{21}'

#################################################################################
################################### Function ####################################
#################################################################################

def filter_fasta(input, module, window_st = 14, window_ed = 16):
    all_results = []
    for seq_record in SeqIO.parse(input, "fasta"):
        #sgRNAs = re.findall( pattern, str(seq_record.seq), re.I)
        sgRNAs = re.findall(pattern, str(seq_record.seq), overlapped=True)
        seq_record.id = seq_record.id.replace("::","_").replace(":","_")
        if sgRNAs:
            #print seq_record.id, "; ", seq_record.seq, ";", str(sgRNAs)
            sgRNAF = []
            sgRNAR = []
            BE = []
            for sgRNA in sgRNAs:
                if re.match(patternF, sgRNA):
                    seed_st = flank - window_ed
                    seed_ed = seed_st + (window_ed - window_st + 1)
                    seed = sgRNA[seed_st:seed_ed]
                    sgRNAF.append(seed)
                    if module == "ABE" and "A" in seed:
                        BE.append("yes")
                    elif module == "ABE" and "A" not in seed:
                        BE.append("no")
                    elif module == "CBE" and "C" in seed:
                        BE.append("yes")
                    elif module == "CBE" and "C" not in seed:
                        BE.append("no")
                    else:
                        pass
                elif re.match(patternR, sgRNA):
                    seed_st = window_ed
                    seed_ed = seed_st + (window_ed - window_st + 1)
                    seed = sgRNA[seed_st:seed_ed]
                    sgRNAR.append(seed)
                    if module == "ABE" and "T" in seed:
                        BE.append("yes")
                    elif module == "ABE" and "T" not in seed:
                        BE.append("no")
                    elif module == "CBE" and "G" in seed:
                        BE.append("yes")
                    elif module == "CBE" and "G" not in seed:
                        BE.append("no")
                    else:
                        pass
            if not "yes" in BE:
                #print seq_record.id, "; ", seq_record.seq, ";", str(sgRNAs), ";", sgRNAF, ";", sgRNAR, ";", BE, ";", "delete"
                pass
            else:
                #print seq_record.id, "; ", seq_record.seq, ";", str(sgRNAs), ";", sgRNAF, ";", sgRNAR, ";", BE, ";", "pass"
                results = seq_record.id + "^" + str(seq_record.seq)
                #print results
                all_results.append(results)
                #filter_output.write(results + "\n")
    return(all_results)


def save_to_file(data, output):
    with open(output, 'w') as filehandle:
        filehandle.writelines("%s\n" % line for line in data)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, help = "The input fasta file, can have path in file name")
    parser.add_argument("-o", "--output", type = str, help = "The output file to run for BE_Designer, can have path in file name")
    parser.add_argument("--module", type = str, help = "The BE editor type, can be ABE or CBE", default = "ABE")
    parser.add_argument("--window_st", type = int, help = "start position of base editing window (14 for ABE and 13 for CBE) (default = 14)", default = 14)
    parser.add_argument("--window_ed", type = int, help = "end position of base editing window (16 for ABE and 17 for CBE) (default = 16)", default = 16)

    return parser.parse_args()


def main():

    args = parse_args()

    input_fasta = args.input
    output_results = args.output
    module = args.module
    window_st = args.window_st
    window_ed = args.window_ed

    filter_input_fasta = filter_fasta(input = input_fasta, module = module, window_st = window_st, window_ed = window_ed)
    save_to_file(data = filter_input_fasta, output = output_results)


if __name__ == '__main__':
    main()