#!/usr/bin/env python3

import math
import pysam
import tqdm
import numpy as np
import re
import pandas as pd
import os
import argparse
import gzip

def find1(stro, ch):
  # 0.100 seconds for 1MB str
  npbuf = np.frombuffer(bytes(stro,'utf-8'), dtype=np.uint8) # Reinterpret str as a char buffer
  return np.where(npbuf == ord(ch))[0]

def convert_to_coordinate(seq, Ml, Mm, which=u"T"):
    assert (len(Ml) == len(Mm))

    Mm = Mm + 1
    result = np.zeros((len(seq))) + np.nan
    n_which = 0
    
    cum = np.cumsum(Mm)
    cum -= 1
    
    if which != "N":
        pos = find1(seq, which)
    else:
        pos = np.arange(len(result))

    if len(cum) != 0:
        if cum[-1] > len(pos) - 1:
            cum = cum[:np.argmax(cum > len(pos) - 1)]
            Ml = Ml[:len(cum)]

        result[pos[cum]] = np.array(Ml) / 255

    return result

def smooth(ser, sc):
    return np.array(pd.Series(ser).rolling(sc, min_periods=1, center=True).mean())

def process_bam(bam, verbose=False, fill_nan=False, res=1, maxi=None, chs=None, no_seq=False, remove_less_than=None, tqdm_do=False):
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(bam, "r")
    pysam.set_verbosity(save)

    Read = {}
    
    monitor = lambda x : x
    if tqdm_do:
        monitor = lambda x: tqdm.tqdm(x)
    for ir, read in monitor(enumerate(samfile)):
        if verbose:
            print(read)

        seq = read.get_forward_sequence() # return the actual sequence of the mapped read
        attr = {}

        if no_seq:
            attr["mapped_strand"] = "+"
        else:
            if read.is_reverse:
                attr["mapped_strand"] = "-"
            else:
                attr["mapped_strand"] = "+"

            attr["mapped_chrom"] = read.reference_name
            pos = read.get_reference_positions()

            try:
                attr["mapped_start"] = pos[0]
                attr["mapped_end"] = pos[-1]
                attr["seq"] = seq
                attr["cigar"] = read.cigarstring # For InDel handling
            except:
                pass

            if chs is not None and attr["mapped_chrom"] not in chs:
                continue

        for tag in ["Ml","ML"]:
            if read.has_tag(tag):
                Ml = read.get_tag(tag)
                break

        for tag in ["Mm","MM"]:
            if read.has_tag(tag):
                Mmt = read.get_tag(tag).split(";")[:-1]
                break

        Mm = {}
        base_ref = {}
        skipped_base = {}
        for Smm in Mmt:
            mod = Smm[2:3].upper()
            try:
                shift = np.fromstring(Smm[Smm.index(",")+1:], dtype=np.int32, sep=',')
            except:
                print("Strange read skipping")
                print(read)
                print(Smm)
                continue

            Mm[mod] = shift
            base_ref[mod] = Smm[:1]
            skipped_base[mod] = Smm[3:4]

        if not Mm:
            continue  # Avoid downstream NameError for undefined 'mod'
        
        Nn = {}
        start = 0
        for mod in Mm.keys(): # Nn[mod] not working
            val = convert_to_coordinate(seq, Ml[start:start+len(Mm[mod])], Mm[mod], which=base_ref[mod])
            start += len(Mm[mod])

            # if attr["mapped_strand"] == "-":
            #     val = val[::-1]

            if skipped_base[mod] == ".":
                val[np.isnan(val) & (np.array(list(seq)) == base_ref[mod])] = 0

            if res != 1:
                n = res * (len(val)//res)
                val = val[:n].reshape(-1, res)
                npts =  np.nansum(~np.isnan(val), axis=1, dtype=np.float16)
                val= np.nanmean(val, axis=1, dtype=np.float16)
            else:
                npts =  np.array(~np.isnan(val), dtype=int)
            Nn[mod] = val
            Nn[mod+"_npoints"] = npts
        # skip=False

        if remove_less_than is not None:
            if Nn == {}:
                continue
            else:
                for m in Mm.keys():
                    if type(remove_less_than) == float:
                        th = remove_less_than
                    else:
                        th = remove_less_than.get(m,0)
                    if np.nanmean(Nn[m]) <= th:
                        continue

        Read[read.query_name] = [attr,Nn]

    return Read



#####################   Main  ##########################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_bam', type=str)
    parser.add_argument('--out_folder', type=str)
    parser.add_argument('--name_files', type=str)

    args = parser.parse_args()
   
    args.out_folder = args.out_folder.rstrip('/') + '/'

    # Create output folder if it does not exist
    os.makedirs(args.out_folder, exist_ok=True)

    reads = process_bam(args.in_bam, verbose=False, fill_nan=False, res=1, maxi=None, chs=None, tqdm_do=False)
    reads_keys = [*reads]
    num_out_file = re.search(r'mod_splitted_(.*?).bam', args.in_bam).group(1)

    # Write .fa file
    fa_filename = args.out_folder + args.name_files + "_" + str(num_out_file) + ".fa.gz"
    with gzip.open(fa_filename, "wt") as ofile:  # 'wt' mode for text mode writing
        for k in range(len(reads_keys)):
            read = reads[reads_keys[k]][0]
            header = (f">/read_{reads_keys[k]} {{'mapped_chrom': '{read['mapped_chrom']}', 'mapped_strand': '{read['mapped_strand']}', 'mapped_start': {read['mapped_start']}, 'mapped_end': {read['mapped_end']}}}")
            ofile.write(header + "\n" + read["seq"] + "\n")

    # Write .fa_ratio file
    fa_ratio_filename = args.out_folder + args.name_files + "_" + str(num_out_file) + ".fa_ratio.gz"
    with gzip.open(fa_ratio_filename, "wt") as ofile:  # 'wt' mode for text mode writing
        for k in range(len(reads_keys)):
            # print(reads[reads_keys[k]][1])
            ofile.write(">/read_" + reads_keys[k] + "\n" + ' '.join(["{0:0.2f}".format(l) for l in reads[reads_keys[k]][1]['B']]) + "\n") # Key could be defined better

    # Write .fa_cigar file
    fa_cigar_filename = args.out_folder + args.name_files + "_" + str(num_out_file) + ".fa_cigar.gz"
    with gzip.open(fa_cigar_filename, "wt") as ofile:  # 'wt' mode for text mode writing
        for k in range(len(reads_keys)):
            # print(reads[reads_keys[k]][1])
            ofile.write(">/read_" + reads_keys[k] + "\n" + str(reads[reads_keys[k]][0]["cigar"]) + "\n")

