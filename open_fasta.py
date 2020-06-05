#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 16:15:50 2020

@author: alix de thoisy
@ Institut Pasteur de Guyane
"""
import sys

from Bio import SeqIO
from Bio.Alphabet import Alphabet

def open_fasta (file_path, ab=Alphabet()):
    """
    :param file_path: file of fasta file
    :param *alphabet: optional, alphabet used for sequences
                    "dna" or "protein"
    :type path_aln: String
    :type *alphabet: String

    :return: bioseq, seq_id list of sequence id
    :rtype: list of BioPython sequences, list of String

    :exceptions:
        FileNotFoundError if file could not be opened
        raises FileNotFoundError
    
        other exception:
        raises general Exception
                (used for tracing in the main function)

    .. note:: requires BioPython library
    """
    bioseq = []
    seq_id = []

    try:
        for seq_record in SeqIO.parse(file_path, "fasta", alphabet=ab):
            bioseq.append(seq_record)
            seq_id.append(seq_record.id)
    except FileNotFoundError as not_found_error: 
        raise FileNotFoundError("Error while processing open_fasta function.\
                                 \nFile not found error. {} ".format(not_found_error))
    except err:
        raise Exception("Error while processing open_fasta function.\
                        \nUnable to read fasta file: {}. \nERROR: {} "\
                        .format(file_path, err))


    return bioseq, seq_id