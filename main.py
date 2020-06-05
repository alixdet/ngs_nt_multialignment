#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:12:40 2020

@author: alix de thoisy
@ Institut Pasteur de Guyane
"""

import sys
import time

from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from open_aln import open_aln
from open_fasta import open_fasta
from adding_indels import realign
from prot_to_nt_ma import prot_to_nt_ma


def main(path_aln, path_contig_ma_prot, path_subject_nt, path_contig_nt):
    """
    :param path_aln: path to output of blastx
    :param path_contig_ma_prot: path to fasta file containing the contig 
                                multi-aligned in proteins
    :param path_subject_nt: path to fasta file containing entire subject 
                            sequence in nucleotides
    :param path_contig_nt: path to fasta file containing entire contig 
                            sequence in nucleotides
    :type path_aln: String
    :type path_contig_ma_prot: String
    :type path_subject_nt: String
    :type path_contig_nt: String

    :return: msa_nt
    :rtype: list of BioPython sequences

    .. note:: to have an output, use main_w_output
    """
    #opening the aln file containing the output of blastx
    aln_file = open_aln(path_aln)

    #opening the fasta file containing the contigs
    contig_bioseq_fasta, contig_id_fasta = open_fasta(path_contig_nt, ab=IUPACAmbiguousDNA())
    
    #opening the fasta file containing the subjects
    subject_bioseq_fasta, subject_id_fasta = open_fasta(path_subject_nt, ab=IUPACAmbiguousDNA())

    #running the function to get multiple alignment
    msa_bioseq = realign(aln_file, path_contig_ma_prot)

    """
    from the msa file sequence get all the info required for alignment in the
    several files 
    """
    msa_nt = prot_to_nt_ma(aln_file,
                           subject_bioseq_fasta,
                           contig_bioseq_fasta,
                           msa_bioseq)

    return msa_nt

def main_w_output(path_aln,
                  path_contig_ma_prot,
                  path_subject_nt,
                  path_contig_nt,
                  output_directory):
    """
    :param path_aln: path to output of blastx
    :param path_contig_ma_prot: path to fasta file containing the contig 
                                multi-aligned in proteins
    :param path_subject_nt: path to fasta file containing entire subject 
                            sequence in nucleotides
    :param path_contig_nt: path to fasta file containing entire contig 
                            sequence in nucleotides
    :param output_directory: path of output directory
    :type path_aln: String
    :type path_contig_ma_prot: String
    :type path_subject_nt: String
    :type path_contig_nt: String
    :type output_directory: String

    :print: processing time in seconds
    :output: fasta file named ..._nt_ma.msa with multi-alignment in nucleotids
              in output directory
    :return: nothing

    .. note:: to be run in command line
    """
    t0 = time.time()
    msa_nt = main(path_aln,
                  path_contig_ma_prot,
                  path_subject_nt,
                  path_contig_nt) 
    """
    outputing the result
    """
    aln_name = path_aln.split("/")[-1][:-4]
    output_file  = output_directory + "/" + aln_name + "_nt_ma.msa"
    with open(output_file, "w+") as f_out:
        for seq in msa_nt:
            SeqIO.write(seq, f_out, "fasta")
    
    t1 = time.time()
    print ("Processing time :" + str(t1-t0)) 
    
if __name__ == "__main__":
    main_w_output(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4])

