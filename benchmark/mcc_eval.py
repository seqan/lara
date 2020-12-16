#!/usr/bin/python3

# Author: Jörg Winkler

# Compute Matthews correlation coefficient (MCC) for two structured alignments.
# Usage: ./mcc_eval.py <reference.fasta> <test.fasta>
# Note: You need to have the RNAalifold program in PATH and Biopython installed.

# If you use this script, please cite:
# Winkler, Urgese, Ficarra, Reinert (2021)
# LaRA 2: Parallel and vectorized program for sequence-structure alignment of RNA sequences

from sys import argv
from subprocess import run
from math import sqrt
from Bio.AlignIO import read


# parse an alignment from file and compute a structure string with RNAalifold
def parse_alignment(filename):
    if filename.endswith("aln"):
        ali = read(filename, 'clustal')
    else:
        ali = read(filename, 'fasta')
    proc = run(["RNAalifold", "--noPS", filename], capture_output=True, check=True, text=True)
    structure = proc.stdout.split('\n')[1].split(' ')[0]
    assert len(structure) == len(ali[0])
    return ali, structure


# retrieve paired alignment positions
def parse_structure(structure):
    paired = []
    stack = []
    for idx in range(len(structure)):
        if structure[idx] == '(':
            stack.append(idx)
        elif structure[idx] == ')':
            paired.append((stack.pop(), idx))
    assert len(stack) == 0
    return paired


# transform alignment positions into sequence positions
def to_seq_pos(aligned_seq, basepairs):
    pos = []
    idx = 0
    for char in aligned_seq:
        if char == '-':
            pos.append('-')
        else:
            pos.append(idx)
            idx += 1
    return [(pos[x], pos[y]) for x, y in basepairs if pos[x] != '-' and pos[y] != '-']


if len(argv) != 3:
    print('    Compute Matthews correlation coefficient (MCC) for two structured alignments.')
    print(f'    Usage: {argv[0]} <reference.fasta> <test.fasta>')
    print('    Note:  You need to have the PETfold program in PATH and Biopython installed.')
    exit(1)

# parse alignments and compute a structure string
refAli, refStructure = parse_alignment(argv[1])
testAli, testStructure = parse_alignment(argv[2])

# retrieve paired alignment positions
refPairs = parse_structure(refStructure)
testPairs = parse_structure(testStructure)

# initialise confusion matrix
tp, fp, fn, tn = 0, 0, 0, 0

for seqidx in range(len(refAli)):
    # convert to sequence positions in order to compare matching base pairs
    refSeqPos = to_seq_pos(refAli[seqidx].seq, refPairs)
    testSeqPos = to_seq_pos(testAli[seqidx].seq, testPairs)

    # match base pairs
    l_tp = len([x for x in testSeqPos if x in refSeqPos])      # intersection
    l_fp = len([x for x in testSeqPos if x not in refSeqPos])  # test - ref
    l_fn = len([x for x in refSeqPos if x not in testSeqPos])  # ref - test
    seqLen = len(refAli[seqidx].seq) - refAli[seqidx].seq.count('-')
    l_tn = seqLen * (seqLen - 1) // 2 - l_tp - l_fp - l_fn      # complement(tp+fp+fn)

    # add to the confusion matrix
    tp += l_tp
    fp += l_fp
    fn += l_fn
    tn += l_tn

# Matthews correlation coefficient
if (tp+fp) == 0 or (tp+fn) == 0 or (tn+fp) == 0 or (tn+fn) == 0:
    mcc = tp * tn - fp * fn
else:
    mcc = (tp * tn - fp * fn) / sqrt((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn))
print(f'{mcc:.4f}')
