#!/usr/bin/env python3
"""
Author : petros <petros@localhost>
Date   : 2021-09-07
Purpose: Global alignment of DNA sequences
"""

import argparse
import os
import sys
from typing import NamedTuple
import numpy
import pandas
from Bio import SeqIO


class Args(NamedTuple):
    """ Command-line arguments """
    DNA1: str
    DNA2: str
    Out : str
    Match: int
    Mismatch: int
    Gap: int


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Alignment of DNA sequences',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('DNA1',
                        metavar='str',
                        help='First DNA sequence (Fasta file path or sequence)'
                        )

    parser.add_argument('DNA2',
                        metavar='str',
                        help='Second DNA sequence (Fasta file path or sequence)'
                        )
    parser.add_argument('-o',
                        '--out',
                        metavar='str',
                        help="output txt file name",
                        type=str,
                        default='alignment_output.txt')

    parser.add_argument('-m',
                        '--mat',
                        metavar='int',
                        help='positive score for base pair match',
                        type=int,
                        default=2
                        )

    parser.add_argument('-n',
                        '--mis',
                        metavar='int',
                        help='negative score for mismatch',
                        type=int,
                        default=-1
                        )
    parser.add_argument('-g',
                        '--gap',
                        metavar='int',
                        help='negative score for the insertion of a gap',
                        type=int,
                        default=-2
                        )

    args = parser.parse_args()
    
        
    
    if os.path.isfile(args.DNA1):
        args.DNA1 = str(SeqIO.read(args.DNA1, format="fasta").seq)
    else:
        if not set(args.DNA1.upper()).issubset("ACGT"):
            parser.error("Invalid characters in DNA sequence")

    if os.path.isfile(args.DNA2):
        args.DNA2 = str(SeqIO.read(args.DNA2, format="fasta").seq)
    else:
        if not set(args.DNA2.upper()).issubset("ACGT"):
            parser.error("Invalid characters in DNA sequence")

    if args.gap >= -1:
        parser.error(f'-g {args.gap} gap penalty must be a negative integer')
    if args.mat < 1:
        parser.error(f'-m {args.mat} match must be a positive integer')
    if args.mis >= 0:
        parser.error(f'-n {args.mis} mismatch penalty must be a negative integer')

    return Args(DNA1=args.DNA1, DNA2=args.DNA2, Match=args.mat, Mismatch=args.mis, Gap=args.gap, Out=args.out)


# --------------------------------------------------
def main() -> None:
    """ The juicy part """

    args = get_args()
    df = array_score_calculator(array_initialization(args.DNA1, args.DNA2, args.Gap), args.Gap, args.Match, args.Mismatch)
    result = traceback(df, args.Gap, args.Match, args.Mismatch)
    sys.stdout = open(args.Out, 'w')
    print(result)
    sys.stdout.close()


# --------------------------------------------------
def array_initialization(DNA1, DNA2, gap) -> pandas.DataFrame:
    """Initialize empty array"""
    print("initializing score matrix")
    cols = ["0",] + list(DNA1)
    rows = ["0",] + list(DNA2)
    pd = pandas.DataFrame(numpy.zeros((len(DNA2) +1, len(DNA1) +1), dtype=int), columns=cols, index=rows)
    pd.iloc[0] = list(range(0, (len(DNA1) * gap) -1, gap))
    pd.iloc[:, 0] = list(range(0, (len(DNA2) * gap) -1, gap))
    return pd

# --------------------------------------------------
def array_score_calculator(alignment_df, gap, match, mismatch) -> pandas.DataFrame:
    """Computes alignment score for each array position"""
    print("calculating scores for matrix positions")
    seq1 = list(alignment_df.columns.values)
    seq2 = list(alignment_df.index)
    for col in range(1, len(list(alignment_df.columns.values))):
        for row in range(1, len(list(alignment_df.index))):
            calcs = [alignment_df.values[row -1 , col -1] + match_function(seq1[col], seq2[row], match, mismatch),
                     alignment_df.values[row - 1, col] + gap, alignment_df.values[row, col-1] + gap]
            alignment_df.values[row, col] = max(calcs)

    return alignment_df

# --------------------------------------------------
def traceback(alignment_df, gap, match, mismatch) -> str:
    """Traceback function finds the best possible alignment"""
    print("determining best alignment")
    Align1, Align2, Bridge = [], [], []
    seq1 = list(alignment_df.columns.values)
    seq2 = list(alignment_df.index)
    i = len(list(alignment_df.columns.values)) - 1
    j = len(list(alignment_df.index)) - 1
    while i >= 1 or j >=1:
        if i >= 1 and j >= 1 and alignment_df.values[j,i] == alignment_df.values[j-1, i-1] + match_function(seq1[i], seq2[j], match, mismatch):
            Align1.insert(0, seq1[i])
            Align2.insert(0, seq2[j])
            if seq1[i] == seq2[j]:
                Bridge.insert(0, "|")
            else:
                Bridge.insert(0, " ")
            i -= 1
            j -= 1
        elif j >= 1 and alignment_df.values[j,i] == alignment_df.values[j - 1, i] + gap:
            Align1.insert(0, "-")
            Align2.insert(0, seq2[j])
            Bridge.insert(0, " ")
            j -= 1
        else:
            Align1.insert(0, seq1[i])
            Align2.insert(0, "-")
            Bridge.insert(0, " ")
            i -= 1
    print("Done!")
    return '{}\n{}\n{}'.format("".join(Align1), "".join(Bridge), "".join(Align2))
# --------------------------------------------------
def match_function(base1, base2, match, mismatch) -> int:
    """Checks if the two nucleotides are the same"""

    if base1 == base2:
        return match
    else:
        return mismatch

# --------------------------------------------------
if __name__ == '__main__':
    main()




###############################
###############################
#Internal unit tests

def test_match_function() -> None:
    assert match_function("A", "A", 1, -2) == 1

def test_array_init() -> None:
    assert len(array_initialization("ATTA", "ATTA", -2)) == 5

def test_array_score() -> None:
    assert array_score_calculator(array_initialization("ATTA", "ATTA", -2), -2, 2, -1).iloc[4][4] == 8

def test_traceback() -> None:
    assert "||||" in traceback(array_score_calculator(array_initialization("ATTA", "ATTA", -2), -2, 2, -1), gap= -2, match= 2, mismatch= -1)

###############################
###############################
