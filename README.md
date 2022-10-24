# idr.mol.feats
Python package to compute bulk molecular features from IDR/protein sequences

Computes a range of molecular features from amino acid sequences of intrinsically disordered protein regions (IDRs). IDRs can either be predicted using state-of-the-art disorder predictors (e.g. SPOT-Disorder), or determined experimentally.

Getting Started

Clone the project from the remote repository: https://github.com/IPritisanac/idrs.mol.feats and start a new branch

Run as: python run_feats.py input_file.txt

Prerequisites

This set of scripts is self-contained and all the dependencies (outside of the standard python packages such as: numpy, scipy, math, etc.) are provided

e.g. run_feats.py is dependent on sequence_features.py and data.py

these are imported in the header:

e.g. from sequence_features import SequenceFeatures


Make sure you have python2.7 or higher versions installed on your system

Description

Classes for calculation of molecular features from protein IDR sequence data

- Parser class (myparser.py) - parses input_file.txt
- input_file.txt -  lists parameters/parameter files
- Alignment class (run_feats.py) - reads in and process IDR alignment sequences
- Features class (sequence_features.py) - contains methods that compute molecular features from sequence files

Contributing

Please read CONTRIBUTING.md for details on our code of conduct, and the process for submitting pull requests to us.

Versioning

For the versions available, see the tags on this repository.

Authors

Alan Moses, Iva Pritisanac & Khaled Elemam
