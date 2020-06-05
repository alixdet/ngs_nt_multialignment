#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 09:59:50 2020

@author: Alix de Thoisy
@Institut Pasteur de la Guyane

"""
import sys
import os

def open_aln(path_file):
    """
  read the alignments in aln file with following pattern:
  aln : query=k141_10006 [ length =279] [frame = 3 ] [9 140] : subject=WP_105045112 [44 1 True] Bacteria
  LRDIDVLYINMDKSVDRRKHMEKQLSRLPVRTTRISGVDGKKLE
  MKTYKVYYINLDKSTERRDFMEKQFQKLDIPITRFSAVYGKELE
  ---
  aln : query=k141_10090_sect:1006_1285 [ length =279] [frame = 2 ] [5 277] : subject=AHN52697 [89 11 True] Viruses
  LPLSGRVPVVYEPTSSAQKTTLVHTTGKIPFAAGNVTATGNAGDQTTGPLLFPQESPSAPTQGFIDLNGTHYADLSQATSTTVNEFRRAIR
  LPLGDRANIVYEHDNVAQYIRL-RSSGGLP------STATDLRNSSTG--IFEQVDNDDV---HLDLNGTHYADLSTATAATINSLRNSFQ
  ----

  :param path_file: path of aln file
  :type path_file: String

  :return: 
    nb_seq : int, number of sequences found

    contig_query : list of contig queries
    contig_frame : list of contig frames
    contig_beg : list of contig start positions (not corrected)
    contig_end : list of contig end positions
    contig_prot : list of contig raw protein sequences

    subject_acc : list of subject accessions
    subject_sense : list of subject reading senses
    subject_beg : list of subject start positions (not corrected)
    subject_end : list of subject end positions
    subject_prot : list of subject raw protein sequences

  :exceptions: exits if file not found
               exits if file empty
               exits if no sequence found
    """
    try :
        with open(path_file,"r") as file_aln :
            aln = file_aln.read()
            
    except FileNotFoundError as file_err:
        raise FileNotFoundError("Open_aln function unable to open file. \n \nError : {}".format(file_err))
    
    if os.path.getsize(path_file) == 0:
        raise IOError("Error while processing open_aln function.\nAln file is empty")

    # splitting the alignments
    # see case if protein sequence ends by '---'
    # nb of alignments in this example : 454
    aln_split = aln.split("aln : query=")[1:]
    
    
    
    """
    extracting info from the list of alignments
    """
    # initialization of lists
    nb_seq = 0
    
    contig_query = []   # str
    contig_frame = []   # str
    contig_beg = []     # int
    contig_end = []     # int
    contig_prot = []    # should be BioPython protein sequence, actually str
    
    subject_acc = []    # str
    subject_sense = []  # bool, str if error
    subject_beg = []    # int
    subject_end = []    # int
    subject_prot = []    # should be BioPython protein sequence, actually str
    
    
    for i in range(len(aln_split)) :
        nb_seq += 1
        # ================================
        #     
        # on contig
        #     
        # ================================
        
        # query
        # 'k141_10283_sect:205_1165' : cut section ?
        q = aln_split[i].split('[')
        q = q[0].split('=')
        contig_query.append(q[0].strip())
        
        # coordinates
        c= aln_split[i].split('[')
        c= c[3].split(']')
        contig_beg0, contig_end0 = c[0].split(' ')
        contig_beg.append(int(contig_beg0))
        contig_end.append(int(contig_end0))
        
        # frame
        f = aln_split[i].split('[')
        f = f[2].split(' ')
        contig_frame.append(int(f[2]))
        
        # ================================
        #     
        # on subject
        #     
        # ================================
        
        # accession number
        a = aln_split[i].split('[')
        a = a[3].split('=')
        subject_acc.append(a[1][:-1])
    
        # sense
        s = aln_split[i].split('[')
        s = s[-1].split(']')
        s = s[0].split(' ')
        if s[-1] == "True" :
            subject_sense.append(bool(True))
        elif s[-1] == "False" :
            subject_sense.append(bool(False))
        else : subject_sense.append(' ERROR ')
    
        # coordinates
        c = aln_split[i].split('[')
        c = c[-1].split(' ')
        c1, c2 = int(c[0]),int(c[1])
        if c2 < c1 : 
            subject_beg.append(c2)
            subject_end.append(c1)
        else :
            subject_beg.append(c1)
            subject_end.append(c2)
            
            
        # ================================
        #     
        # protein sequence
        #     
        # ================================
        p=aln_split[i].split('\n')[1]
        p2=aln_split[i].split('\n')[2]
        #pr = Seq(p,ProteinAlphabet)
        subject_prot.append(p2)
        contig_prot.append(p)
        
    if contig_query == [] or subject_acc == []:
        raise Exception("Error while processing open_aln function.\nNo sequence found in aln file")

    return  (nb_seq,
            contig_query,
            contig_frame,
            contig_beg,
            contig_end,
            contig_prot,
            subject_acc,
            subject_sense,
            subject_beg, 
            subject_end,
            subject_prot)
