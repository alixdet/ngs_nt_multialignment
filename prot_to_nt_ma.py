#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 13:49:40 2020

@author: alix de thoisy
@ Institut Pasteur de Guyane
"""
import sys

from Bio import SeqIO

from open_aln import open_aln
from open_fasta import open_fasta

def prot_to_nt_ma(aln, 
                  subject_bioseq_fasta,
                  contig_bioseq_fasta,
                  msa_bioseq):
    """
    :param aln: aln file opened by open_aln function
    :param subject_bioseq_fasta: subjects sequences in nt
    :param contig_bioseq_fasta: contig sequences in nt
    :param msa_bioseq: multi aligned sequences in protein
    :type aln: list of lists
    :type subject_bioseq_fasta: list of BioPython sequences
    :type contig_bioseq_fasta: list of BioPython sequences
    :type msa_bioseq: list of BioPython sequences

    :return: msa_nt
    :rtype: list of BioPython sequences multi aligned in nt

    .. note:: meant to be used in main function
            to have an output, use prot_to_nt_ma_w_output
    """
    """
    extracting the info from aln file, returned in multiple lists 
    by open_aln function
    the lists have the same indexation
    """
    nb_seq_aln = aln[0]
    contig_query, contig_frame = aln[1:3]
    contig_beg, contig_end = aln[3:5]
    subject_acc, subject_sense = aln[6:8]
    subject_beg, subject_end = aln[8:10]

    #initializing the list
    msa_nt = []
    
    #browsing linearly through the multiple alignment in proteins
    for i in range (len(msa_bioseq)):
        """
        case subject
        """
        if msa_bioseq[i].description.split(" ")[-1] == "subject":
            #getting the index on the aln file
            index_in_aln = subject_acc.index(msa_bioseq[i].id)
             
            #getting the index in fasta file and the bioseq
            for j,seq in enumerate(subject_bioseq_fasta):
                # pattern of seq id : 'AVN87675_ORGANISM_Chicken_megrivirus'
                if seq.id.split("_ORGANISM")[0] == msa_bioseq[i].id:
                    index_in_fasta = j
                    break
            else:
                print(msa_bioseq[i].id)
                sys.exit("subject not found in fasta")
            subject_cut = subject_bioseq_fasta[index_in_fasta]
            
            #processing to alignment
            for j in range(len(msa_bioseq[i])):
                if (msa_bioseq[i][j] == '-'):
                    subject_cut = subject_cut[:j*3]+'---'+subject_cut[j*3:]
            msa_nt.append(subject_cut)
            
            #used for the contig
            previous_acc = msa_bioseq[i].id   
        else:
            """
            contig
            """ 
            #simple version
            if msa_bioseq[i].description.split(" ")[-1] == "contig":
                for j in range(nb_seq_aln):
                    if msa_bioseq[i].id == contig_query[j] and \
                        subject_acc[j] == previous_acc:
                            index_in_aln = j
                            break
                beg_nt = contig_beg[index_in_aln]-1
                end_nt = contig_end[index_in_aln]
                frame = True
                if contig_frame[index_in_aln] < 0:
                    frame = False
                
            #contig id's have been modified when met several times
            #used to store the info within the name
            #pattern: >k141_40040__2 _contig_coord_93_653_True
            elif msa_bioseq[i].description.split("_")[-5] == "contig":
                info = msa_bioseq[i].description.split("_")
                beg_nt = int(info[-3])-1
                end_nt = int(info[-2])
                frame = eval(info[-1])          
            
            #if none of the patterns were met
            else:
                #should return an exception ? 
                print("Contig not found in aln")
                print(msa_bioseq[i].id)
                
            #getting the index in fasta file and the bioseq
            for j,seq in enumerate(contig_bioseq_fasta):
                if seq.id == msa_bioseq[i].id:
                    index_in_fasta = j
                    break
                
                #if several occurences of the same query
                #pattern: k141_11076__2
                elif seq.id == msa_bioseq[i].id[:-3]:
                    index_in_fasta = j
                    break
            
            #creating a new biosequence object
            contig_cut = contig_bioseq_fasta[index_in_fasta]     
            #cutting the sequence
            contig_cut.seq = contig_cut.seq[beg_nt:end_nt]
            #changing the attributes
            contig_cut.id = msa_bioseq[i].id
            contig_cut.description= msa_bioseq[i].description
            
            #applying reverse complementation when the frame is negative
            if not frame:
               contig_cut.seq = contig_cut.seq.reverse_complement()
               
               
            #alignment
            for j in range(len(msa_bioseq[i])):
                if (msa_bioseq[i][j] == '-'):
                    contig_cut = contig_cut[:j*3]+'---'+contig_cut[j*3:]
                    
            if contig_cut.index() != -1:
                msa_nt.append(contig_cut)
            else:
                mas_nt_empty.append(contig_cut)  

    return msa_nt


def prot_to_nt_ma_w_output(path_aln,
                           path_subject_nt,
                           path_contig_nt,
                           msa_bioseq):
    """
    :param aln: path to aln file
    :param subject_bioseq_fasta: path to fasta file containg
                                the  subject sequences in nt
    :param contig_bioseq_fasta: path to fasta file containg
                                the  contig sequences in nt
    :param msa_bioseq: multi aligned sequences in protein, 
                    output of addind_indels function
    :type aln: String
    :type subject_bioseq_fasta: String
    :type contig_bioseq_fasta: String
    :type msa_bioseq: list of BioPython sequences

    :return: nothing
    :output: ..._multi_alignment.msa fasta file
        in same directory containing the multi-alignment 
        in nucleotids

    .. note:: meant to be used in command lines
    """
  
    """
    opening the files
    """
    #aln file containing the output of blastx
    aln = open_aln(path_aln)
    
    #fasta file containing the contigs
    contig_bioseq_fasta, contig_id_fasta = open_fasta(path_contig_nt, "dna")

    #fasta file containing the subjects
    subject_bioseq_fasta, subject_id_fasta = open_fasta(path_subject_nt, "dna")

    #fasta file containing the multiple alignment in proteins
    msa_bioseq = open_fasta(adding_indels_output, "protein")[0]
    
    """
    running the alignment function
    """
    msa_nt = prot_to_nt_ma(aln,
                           subject_bioseq_fasta,
                           contig_bioseq_fasta,
                           msa_bioseq)
        
    """
    outputing the result
    """
    output_file  = path_aln[:-4] + "_multi_alignment.msa"
    with open(output_file, "w+") as f_out:
        for seq in msa_nt:
            SeqIO.write(seq, f_out, "fasta")

if __name__ == "__main__":
    prot_to_nt_ma_w_output(sys.argv[1],
                           sys.argv[2],
                           sys.argv[3],
                           sys.argv[4])
            
