__author__ = "Marina Igual"

import os

"""Program to crate a column "Master Protein" with the name 
   of the protein most frequent in each peptide 
   (counting number of scans)"""


############################ FUNCTIONS ####################################

def accessPath(MasterPath, fileTollok):
    Listofpath = []
    intialList = ["C2", "C3", "C4", "C5"]

    for root, dirs, files in os.walk(MasterPath):
        for dirNames in dirs:

            if dirNames in intialList:

                Level1 = root + os.sep + dirNames

                for root1, dirs1, files1 in os.walk(Level1):

                    for dirNames1 in dirs1:
                        if dirNames1.startswith("TMT"):

                            Level2 = root1 + os.sep + dirNames1 + os.sep + "cXcorr_Len_Rank_Results" + os.sep + "Vertex" + os.sep + str(fileTollok)

                            print (Level2)

                            if Level2 not in Listofpath:
                                Listofpath.append(Level2)

    print (len(Listofpath))

    return Listofpath

def MasterProtein(masterFiles, names, IO):
    
    FileList = []
    for master in masterFiles:
        files = accessPath(master, names)
        
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
    name = names.split(".")
    for file in FileList:
        print(file)
        parts_file = file.split("\\")
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

############################ MAIN ##################################
""" Program to match all the peptides with same sequence and the most abundant protein_id.
     The program works with SHIFTS files, they could be generated with or without isotopic mass correction.
     The program output returns a file with the same structure of the input plus an extra column with the most abundand protein ID """
     
     #PROGRAM INPUTS

#MasterFile => Write one or more path to the folder where are located all the sub-folders of SHIFTS 
masterFile = [r"S:\U_Proteomica\PROYECTOS\PESA_omicas\Comet-PTM-2a-5ta_Cohortes_V1\DigPar_cometonly\SHIFTS_C2-3-4-5_V1_FDRper-experiment",
              r"S:\U_Proteomica\PROYECTOS\PESA_omicas\Comet-PTM-2a-5ta_Cohortes_V2\DigPar_cometonly\SHIFTS_C2-3-4-5_V2_FDRper-experiment"]

#fileName => Name of the SHIFT output file. All files must have the same name
fileName = "IsotopCorrection_TargetData_withSequence-massTag.txt"

MasterProtein(masterFile, fileName)




