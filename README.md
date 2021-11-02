# TurtleAlign

TurtleAlign is a functional Python implementation of the much beloved Needleman–Wunsch Global alignment algorithm. This code is meant for educational purposes,
namely to gain a deeper insight into sequence bioinformatics. If you plan on using this to align actual long gene sequences, bring with you a cup of coffee and plenty of patience.

TurtleAlign is a tool to visualize how two sequences align given a match score, a gap and a mismatch penalty.

TurtleAlign supports command line functionality, argument parsing and outputs the alignment in a txt file. To visualize the resulting txt file, make sure to disable line wrapping.
For more information, run the './TurtleAlign.py -h' command in your terminal

The code comes with two bacterial gene fasta files to allow for illustration (I picked some very short ones on purpose, do not worry).

## Documentation

DNA1: First positional argument, takes a DNA sequence directly in the terminal or a DNA fasta file.

DNA2: Second positional argument,  takes a DNA sequence directly in the terminal or a DNA fasta file.

-o, --out: Output file name

-m, --mat: Positive score for base pair match.  Strictly positive integer value

-n, -mis: Negative score for base pair mismatch. Strictly negative integer value

-g, --gap: Negative score for the insertion of a gap. Strictly inferior to -1

## Tutorial

Open the command line terminal and navigate to the TurtleAlign directory

First, ensure installation of the necessary packages: pip install -r requirements.txt

Run the following command and open the outputed txt file (named alignment_output.txt by default): ./TurtleAlign.py -m 2 -n -1 -g -2 MTuber_atpE.fasta EColi_atpE.fasta

If the txt file looks funny, disable line wrapping in the settings of your text editor and scroll along the sequence from end to end.

## Unit tests

Unit tests have been added for all functions in case new features are added. If features are added and functions changed, make sure to run pytest -xv TurtleAlign.py to verify that
existing functions work as intended. To ensure that the overall code runs and fails gracefully, navigate inside the tests directory and run pytest:

cd tests
pytest -xv TurtleAlign_tests.py

The tests will fail unless you are in the tests directory when you run the pytest command.

Tests can be added or modified accordingly, as new features are added.
