""" Tests for Sequence_aligner.py """

import os
import platform
import re
from subprocess import getstatusoutput, getoutput

PRG = '../TurtleAlign.py'
RUN = f'python {PRG}' if platform.system() == 'Windows' else PRG
TEST1 = ('./input_data/sequence1.fasta', './input_data/sequence2.fasta')


# --------------------------------------------------
def test_exists() -> None:
    """ Program exists """

    assert os.path.isfile(PRG)
    
# --------------------------------------------------

def test_help_works() -> None:
    """ -h option prints help page """
    rv, out = getstatusoutput('{} -h'.format(RUN))
    assert rv == 0
    assert out.lower().startswith('usage:')
    assert "Alignment" in out
    
# --------------------------------------------------

def test_bad_input() -> None:
    """ fails on bad DNA input """
    rv, out = getstatusoutput('{} ATTAG AFATTA'.format(RUN))
    assert rv != 0
    assert "Invalid characters in DNA" in out

# --------------------------------------------------

def test_bad_gap() -> None:
    """ fails on bad penalty gap """
    rv, out = getstatusoutput('{} -g 5 ATTAG AATTA'.format(RUN))
    assert rv != 0
    assert "gap penalty must" in out
    
# --------------------------------------------------

def test_bad_match() -> None:
    """ fails on bad match score """
    rv, out = getstatusoutput('{} -m -5 ATTAG AATTA'.format(RUN))
    assert rv != 0
    assert "match must be" in out  

# ---------------------------------------------------

def test_bad_mismatch() -> None:
    """ fails on bad mismatch score """
    rv, out = getstatusoutput('{} -n 4 ATTAG AATTA'.format(RUN))
    assert rv != 0
    assert "mismatch penalty must be" in out  

# ---------------------------------------------------

def test_good_run() -> None:
    """ Runs successfully with good fasta file input """
    rv, out = getstatusoutput('{} {} {}'.format(RUN, TEST1[0], TEST1[1]))
    assert rv == 0
    assert 'TGAATTCAGTTA\n|| | || |  |\nTGGA-TC-G--A\n' == open("alignment_output.txt", "r").read()
    os.remove("alignment_output.txt")
    
    
    
    
