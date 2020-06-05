#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:20:20 2020

@author: alix de thoisy
@ Institut Pasteur de Guyane
"""
def indels_front(subject_start, subject, entire):
    """
    inputs
    :subject_start: integer, start position on the entire
    :subject: str, subject sequence
    :entire: str, entire sequence
    
    returns the nomber of indels that need to be added in front of the 
    SUBJECT sequence in order to have an alignement with the entire
    sequence
    
    could also take previous_position as parameter to not browse through
    the same parts
    """
    #number of characters met so far
    no_char = 0
    index = 0
    while no_char < subject_start:
        if entire[index] !="-":
            no_char += 1
        index += 1
        
    return index

def align_1subject_1entire (no_indels_front, subject, entire):
    """    
    Parameters
    ----------
    no_indels_front : integer
        DESCRIPTION.
    subject : str
        DESCRIPTION.
    entire : str
        DESCRIPTION.

    Returns
    -------
    subject : str
        DESCRIPTION.
    entire : str
        DESCRIPTION.

    """
    #index of indels that are added on the ENTIRE sequence
    added_indels = []

    j = 0    
    while j < (len(subject))-1:
        if entire[no_indels_front + j] == "-" and subject[j] != "-":
            subject = subject[:j]+"-"+subject[j:]
        if subject[j] == "-" and entire[no_indels_front + j] != "-":
            entire = entire[:(no_indels_front + j)]\
                +"-" \
                +entire[(no_indels_front + j):]
            added_indels.append(no_indels_front + j)
            
        j += 1
    
    return entire, subject, added_indels
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    