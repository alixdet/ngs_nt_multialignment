#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:15:58 2020

@author: alix de thoisy
@Institut Pasteur de la Guyane
"""
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA, ExtendedIUPACProtein

from open_aln import open_aln
from open_fasta import open_fasta
from realigning_functions import indels_front

def realign(aln_file, path_contig_ma_prot):
    """
    :param aln_file: aln file opened by open_aln function
    :param path_contig_ma_prot: path to fasta file containing
                            the contigs multi-aligned in proteins
    :type aln: list of lists
    :type path_contig_ma_prot: String

    :return: msa_prot
    :rtype: list of BioPython sequences multi aligned in nt

    .. note:: requires BioPython library
            meant to be used in main function
            to have an output, use realign_w_output
    """
    #initialization
    msa_prot= []

    #opening fasta file
    prot_bioseq, prot_id0 = open_fasta(path_contig_ma_prot, ab=ExtendedIUPACProtein())

    prot_id = []
    for i in range(len(prot_id0)):
        prot_id.append(prot_id0[i].split('_ORGANISM')[0])

    """	
    opening the aln file and creating dictionaries	
    """
    nb_seq_aln = aln_file[0]
    contig_query = aln_file[1]
    contig_frame = aln_file[2]
    contig_beg, contig_end = aln_file[3:5]
    contig_prot = aln_file[5]
    subject_acc = aln_file[6]
    subject_beg = aln_file[8]
    subject_prot = aln_file[10]

    #some contig queries are the same so we change the second name
    used_queries = []
    for i in range(len(contig_query)):
        if contig_query[i] in used_queries:
            n = used_queries.count(contig_query[i])
            used_queries.append(contig_query[i])
            #defining the frame
            frame = True
            if contig_frame[i] < 0: 
                frame = False
            contig_query[i] = contig_query[i]+"__"+str(n+1)+"__coord_"+\
                str(contig_beg[i])+"_"+str(contig_end[i])+"_"+str(frame)
        else :
            used_queries.append(contig_query[i])
            
    """
    we build nested dictionaries, one dictionary per protein acc
    each one contains one dictionary per subject_query having
    start, contig_prot, contig_end as keys
    aln name of the main dict
    
    aln = {
            acc1 = {
                "entire",
                query1 = {"start",
                          "contig",
                          "subject"},
                query2:
                ...
                    },
            acc2 = { 
                    ...
                    },
            ...
                
        }
    """
    aln = {}

    #find all unique subject_acc
    unique_acc = []
    for acc in subject_acc:
        if acc not in unique_acc:
            unique_acc.append(acc) 

    #creating the dictionaries
    for acc in unique_acc:
        if acc in prot_id:
            acc_in_fasta = prot_id.index(acc)
            aln.update ({ 
                acc : {"entire" : prot_bioseq[acc_in_fasta].seq}
                })
            indexes = [i for i in range(nb_seq_aln) if subject_acc[i] == acc]
            for index in indexes:
                aln[acc].update({
                    contig_query[index]: {
                    "start": subject_beg[index] - 1,
                    "contig": contig_prot[index],
                    "subject": subject_prot[index]
                        }
                    })
    """
    counting the number of indels that need to be added in front of the partial 
    sequence
    using indels_front function from realigning_functions module
    """
    used_acc = []
    for acc in aln:
        used_queries = []
        for query in aln[acc]:
            if query != "entire":
                aln[acc][query]["indels_front"] = indels_front(aln[acc][query]["start"],
                                                               aln[acc][query]["subject"],
                                                               aln[acc]["entire"])
                #adding the indels of this query on both the entire and partial sequence
                j = 0
                no_indels_front = aln[acc][query]["indels_front"]
                while j <= (len(aln[acc][query]["subject"]))-1:
                    if aln[acc]["entire"][no_indels_front + j] == "-" and aln[acc][query]["subject"][j] != "-":
                        aln[acc][query]["subject"] = aln[acc][query]["subject"][:j]+"-"+aln[acc][query]["subject"][j:]
                        aln[acc][query]["contig"] = aln[acc][query]["contig"][:j]+"-"+aln[acc][query]["contig"][j:]
                    if aln[acc][query]["subject"][j] == "-" and aln[acc]["entire"][no_indels_front + j] != "-":
                        aln[acc]["entire"] = aln[acc]["entire"][:(no_indels_front + j)]\
                            +"-" \
                            +aln[acc]["entire"][(no_indels_front + j):]                                        
                        #we add the same indel on the previous queries from this acc
                        for prev_query in used_queries:
                            aln[acc][prev_query]["subject"] = \
                                aln[acc][prev_query]["subject"][:(no_indels_front + j)] \
                                    + "-" \
                                    + aln[acc][prev_query]["subject"][(no_indels_front + j):]
                            aln[acc][prev_query]["contig"] = \
                                aln[acc][prev_query]["contig"][:(no_indels_front + j)] \
                                    + "-" \
                                    + aln[acc][prev_query]["contig"][(no_indels_front + j):]
                                
                        #we add the same indel on the entires of acc not processed
                        for acc0 in aln:
                            if acc0 != acc and acc0 not in used_acc:
                                        aln[acc0]["entire"] = aln[acc0]["entire"][:(no_indels_front + j)] \
                                            + "-" \
                                            + aln[acc0]["entire"][(no_indels_front + j):]
      
                        #we add the same indel everywhere on the sequences we already processed
                        for acc0 in aln:
                            if acc0 != acc and acc0 in used_acc:
                                for query0 in aln[acc0]:
                                    if query0 == "entire":
                                        aln[acc0][query0] = aln[acc0][query0][:(no_indels_front + j)] \
                                            + "-" \
                                                + aln[acc0][query0][(no_indels_front + j):]
                                    else :   
                                            aln[acc0][query0]["subject"] = \
                                            aln[acc0][query0]["subject"][:(no_indels_front + j)] \
                                                + "-" \
                                                 + aln[acc0][query0]["subject"][(no_indels_front + j):]
                                            aln[acc0][query0]["contig"] = \
                                            aln[acc0][query0]["contig"][:(no_indels_front + j)] \
                                                 + "-" \
                                                 + aln[acc0][query0]["contig"][(no_indels_front + j):]
                    j += 1
                #counting the indels at the back of subject
                dico = aln[acc][query]
                dico["indels_back"] = len(aln[acc]["entire"]) - \
                    (dico["indels_front"] + len(dico["subject"]))
                aln[acc][query] = dico
                    
                #adding the indels on the subject
                dico = aln[acc][query]
                dico["subject"] = "-" * dico["indels_front"] + dico["subject"] \
                                + "-" * dico["indels_back"]
                dico["contig"] = "-" * dico["indels_front"] \
                        + dico["contig"] \
                        + "-" * dico["indels_back"]
                used_queries.append(query)
                aln[acc][query] = dico
        used_acc.append(acc)


    """ 
    converting to BioPython protein sequences 
    """
    for acc in aln:
        aln[acc]["entire_bioseq"] = SeqRecord(Seq(str(aln[acc]["entire"]),
                                              Alphabet), 
                                              id=acc, 
                                              name="No name",
                                              description="subject")
        msa_prot.append(aln[acc]["entire_bioseq"])
        for query in aln[acc]:
            if query != "entire" and query != "entire_bioseq" and "coord" not in query:
                dico = aln[acc][query]
                dico["contig_bioseq"] = SeqRecord(Seq(dico["contig"],Alphabet),
                                                  id=query,
                                                  name="No name",
                                                  description="contig")
                aln[acc][query] = dico
                msa_prot.append(aln[acc][query]["contig_bioseq"])
                            
            elif query != "entire" and query != "entire_bioseq":
                dico = aln[acc][query]
                dico["contig_bioseq"] = SeqRecord(Seq(dico["contig"],Alphabet),
                                  id=(query.split("__coord")[0]),
                                  name="No name",
                                  description="_contig_" + query.split("__")[-1])
                aln[acc][query] = dico
                msa_prot.append(aln[acc][query]["contig_bioseq"])
                
    return msa_prot


def realign_w_output(path_aln,
                     path_contig_ma_prot,
                     output_directory):
    """
    :param path_aln: path to aln file
    :param path_contig_ma_prot: path to fasta file containing
                            the contigs multi-aligned in proteins
    :param output_directory: output directory
    :type aln: String
    :type path_contig_ma_prot: String
    :type output_directory: String

    :return: nothing
    :output: fasta file called "adding_indels_output.msa"
            multi aligned sequences in protein

    .. note:: requires BioPython library
            meant to be used in command lines
    """
    #opening aln file
    aln_file = open_aln(path_aln)
    
    #processing to alignment
    seq_list = realign(aln_file, path_contig_ma_prot)
    
    #writing
    output_file  = output_directory + "/adding_indels_output.msa"
    with open(output_file, "w+") as f_out:
        for seq in seq_list:
            SeqIO.write(seq, f_out, "fasta")


if __name__ == "__main__":
    realign_w_output(sys.argv[1], sys.argv[2], sys.argv[3])
