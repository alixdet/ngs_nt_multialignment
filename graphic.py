#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:07:56 2020

@author: alix de thoisy
@Institut Pasteur de la Guyane
"""
#unable to be imported by Pyinstaller
import pkg_resources.py2_warn

import time

import tkinter as tk
from tkinter.filedialog import askopenfilename, askdirectory
from tkinter import messagebox

from main import main_w_output

"""
defining the several functions to get the paths
"""

def open_aln_file():
    aln_file = askopenfilename(title="Select aln file",
                               filetypes=[('aln file', '.aln'),
                                          ('text file','.txt'),
                                          ('all files','.*')])
    aln_file_var.set(aln_file)
    aln_file_print.set(aln_file.split("/")[-1])
    
def open_contig_ma_prot_file():
    path_contig_ma_prot = askopenfilename(title="Select fasta file",
                               filetypes=[('fasta file', '.fas'),
                                          ('text file','.txt'),
                                          ('all files','.*')])
    path_contig_ma_prot_var.set(path_contig_ma_prot)
    path_contig_ma_prot_print.set(path_contig_ma_prot.split("/")[-1])
    
def open_subject_nt_file():
    path_subject_nt = askopenfilename(title="Select fasta file",
                               filetypes=[('fasta file', '.fas'),
                                          ('text file','.txt'),
                                          ('all files','.*')])
    path_subject_nt_var.set(path_subject_nt)
    path_subject_nt_print.set(path_subject_nt.split("/")[-1])
    
def open_contig_nt_file():
    path_contig_nt = askopenfilename(title="Select fasta file",
                               filetypes=[('fasta file', '.fas'),
                                          ('text file','.txt'),
                                          ('all files','.*')])
    path_contig_nt_var.set(path_contig_nt)
    path_contig_nt_print.set(path_contig_nt.split("/")[-1])

    
def chooseDirectory():
    directory = askdirectory()
    directory_var.set(directory)
    print(directory)

    
def align(path_aln,
          path_contig_ma_prot,
          path_subject_nt,
          path_contig_nt,
          output_directory):
    if aln_file_print.get() == "No aln file selected":
        messagebox.showerror("Error", "Aln file missing")
    elif path_contig_ma_prot_print.get() == "No file selected":
        messagebox.showerror("Error", "Multi alignment file missing")
    elif path_subject_nt_print.get() == "No file selected":
        messagebox.showerror("Error", "Subject file missing")
    elif path_contig_nt_print.get() == "No file selected":
        messagebox.showerror("Error", "Contig file missing")   
    else:
        try :
            t0 = time.time()
            main_w_output(path_aln,
                  path_contig_ma_prot,
                  path_subject_nt,
                  path_contig_nt,
                  output_directory)
            t1 = time.time()
            messagebox.showinfo("Success !","Succesfully done \nRunning time : {} seconds".format(t1-t0))
        except Exception as e: 
            messagebox.showerror("Error", str(e))
      
#main window        
root = tk.Tk()
root.title("Realign")
root.geometry("700x175")
root.minsize(700,175)
root.resizable(True, False)


"""
initializing all the StringVar
"""
aln_file_var = tk.StringVar()
aln_file_print = tk.StringVar()
aln_file_print.set("No aln file selected")

path_contig_ma_prot_var = tk.StringVar()
path_contig_ma_prot_print = tk.StringVar()
path_contig_ma_prot_print.set("No file selected")

path_subject_nt_var = tk.StringVar()
path_subject_nt_print = tk.StringVar()
path_subject_nt_print.set("No file selected")

path_contig_nt_var = tk.StringVar()
path_contig_nt_print = tk.StringVar()
path_contig_nt_print.set("No file selected")

directory_var = tk.StringVar()
#by default set as "." so that it can be run in this directory if nothing
#specified
directory_var.set(".")

"""
buttons and labels
"""
#aln
tk.Button(root, text ="Select aln file", command = open_aln_file)\
    .grid(row=0, column=0, padx = 0)
tk.Label(root, textvariable = aln_file_print).grid(row=0, column=1)

#multi alignment fasta file
tk.Button(root,
          text ="Select multi alignment file",
          command = open_contig_ma_prot_file).grid(row=1, column=0)
tk.Label(root, textvariable = path_contig_ma_prot_print).grid(row=1, column=1)

#subjects in nt fasta file
tk.Button(root,
          text ="Select subject file",
          command = open_subject_nt_file).grid(row=2, column=0)
tk.Label(root, textvariable = path_subject_nt_print).grid(row=2, column=1)

#contigs in nt fasta file
tk.Button(root,
          text ="Select contig file",
          command = open_contig_nt_file).grid(row=3, column=0)
tk.Label(root, textvariable = path_contig_nt_print).grid(row=3, column=1)


#output directory
tk.Button(root, text="Select output directory", command = chooseDirectory)\
    .grid(row=4, column=0)
tk.Label(root, textvariable = directory_var).grid(row=4, column=1)    

#whitespace
tk.Label(root, text ="                                                     ")\
    .grid(row=5, column=1)

#realign button    
tk.Button(root, text="Realign",
          command= lambda: align(aln_file_var.get(),
                                 path_contig_ma_prot_var.get(),
                                 path_subject_nt_var.get(),
                                 path_contig_nt_var.get(),
                                 directory_var.get())
          ).grid(row=6, column=2)

         
root.mainloop()
