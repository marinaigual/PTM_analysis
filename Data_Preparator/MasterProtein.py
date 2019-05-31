__author__ = "Marina Igual"

import os
from pathlib import Path

"""Program to crate a column "Master Protein" with the name 
   of the protein most frequent in each peptide 
   (counting number of scans)"""


############################ FUNCTIONS ####################################

def accessPath(MasterPath, fileTollok, subfolders):
    list_final = []
    file_cod = "**/" + fileTollok
    results = list(Path(MasterPath).glob(file_cod))
    
    for result in results:
        result = str(result)
        elements = result.split("\\")
        n = -1 - len(subfolders)
        accept = 0
        for i in range(len(subfolders)):
            subfol = elements[n]
            n += 1
            if subfol == subfolders[i]:
                accept += 1
            if accept == len(subfolders):
                list_final.append(result)
    print("Total number of files: ", len(list_final))
    return list_final

def MasProt(masterFiles, names, subf):
    
    FileList = []
    for master in masterFiles:
        files = accessPath(master, names, subf)
        
        for file in files:
            FileList.append(file)

    SeqProt = {}
    for file in FileList:
        file_in = open(file, 'r')
        next(file_in)
        for line in file_in:
            if line != "\n":
                splits2 = line.split("\t")
                seq = splits2[7].strip()
                prot_id = splits2[24].strip()
                
                try:
                    SeqProt[seq][prot_id]
                    SeqProt[seq][prot_id] += 1
                except:
                    try:
                        SeqProt[seq]
                        SeqProt[seq][prot_id] = 1
                    except:
                        SeqProt[seq]={}
                        SeqProt[seq][prot_id] = 1
        file_in.close()
    
    SeqProt_def = {}      
    for seq in SeqProt:
        proteins = sorted(SeqProt[seq].items(), key=lambda x: x[1], reverse=True)
        SeqProt_def[seq] = proteins[0][0]
    
    #SeqProt_def is a dictionary with the good protein id
    final_list=[]
    name = names.split(".")
    for file in FileList:
        print(file)
        file_str = str(file)
        parts_file = file_str.split("\\")
        parts_file2 = parts_file[:-1]
        title = str()
        for part in parts_file2:
            title = title + os.sep + part
        complete_name = name[0] + "_MasterProtein.txt"
        title_file = title[1:] + os.sep + complete_name
    
        file_out = open(title_file, 'w')
        file_in = open(file, 'r')
        
        header = next(file_in)
        header = header.strip()
        header = header +  "\t" + "MasterProtein"
        file_out.write(header)
        file_out.write("\n")
        
        for line2 in file_in:
            if line2 != "\n":
                splits2 = line2.split("\t")
                seq = splits2[7].strip()
                prot = SeqProt_def[seq]
                
                for i in range(len(splits2)):
                    file_out.write(splits2[i].strip())
                    file_out.write("\t")
            file_out.write(prot)
            file_out.write("\n")
                
        file_in.close()
        file_out.close()
        
        final_list.append(title_file)
    
    return(final_list)
############################ MAIN ##################################
""" Program to match all the peptides with same sequence and the most abundant protein_id.
     The program works with SHIFTS files, they could be generated with or without isotopic mass correction.
     The program output returns a file with the same structure of the input plus an extra column with the most abundand protein ID """






