__author__ = "Marina Igual"

""" Program to create a relation table with 
    all the peptides and the protein description.
    The program works with the MasterProtein output files"""
    
### IMPORTS
from pathlib import Path

################ FUNCTIONS #######################


def Relation_File(listfoFiles, out, isotopic):

    file_out = open(out, 'w')
    
    if isotopic == "Yes":
        dp = 28 #Column for deltaPeptide
        lab = 18 #Column for Target or Decoy
        tg = 23 #Column for tags
        fasta = 31 #Column with masterProtein
    else:
        dp = 27 #Column for deltaPeptide
        lab = 18 #Column for Target or Decoy
        tg = 23 #Column for tags
        fasta = 28 #Column with masterProtein
        
    file_out = open(out, "w")  
    
    for filename in listfoFiles:
        with open(filename) as f2:
            next(f2)
            for line2 in f2:
                if line2 != "\n":
                    splits2 = line2.split("\t")
                    deltaPEP = splits2[dp].strip()  #For isotopic correction file 28
                    label = splits2[lab].strip()
                    tags = splits2[tg].strip()
                    fasta_description = splits2[fasta].strip() #For isotopic correction file 31

                    if label == "Target" and tags == "NA":
                        file_out.write(fasta_description)
                        file_out.write('\t')
                        file_out.write(deltaPEP)
                        file_out.write('\n')
    file_out.close()


################ MAIN ######################
#
##masterFiles => Path to one or more master folders. 
#
#masterFiles = ["S:\\U_Proteomica\\PROYECTOS\\PESA_omicas\\Comet-PTM-2a-5ta_Cohortes_V1\\DigPar_cometonly\\SHIFTS_C2-3-4-5_V1_FDRper-experiment",
#               "S:\\U_Proteomica\\PROYECTOS\\PESA_omicas\\Comet-PTM-2a-5ta_Cohortes_V2\\DigPar_cometonly\\SHIFTS_C2-3-4-5_V2_FDRper-experiment"]
#
##fileName => Name of the files to create the table deltap2q. It must be enter the file and sub-folder
#filename = "Vertex/IsotopCorrection_TargetData_withSequence-massTag_MasterProtein.txt"
#
##fileOut => Path and name of the output file
#fileOut = r"S:\U_Proteomica\UNIDAD\DatosCrudos\laboratorio\Marina_Igual\HIPERMOD\p2site\PESA_V1V2_deltap2q_CorrMass.txt"
#
##IsotopicCorrection => Specify if the files come from the last SHIFTS module: Yes or No
#IsotopicCorrection = "Yes"
#
#Relation_File_DeltaP2Q(masterFiles, filename, fileOut, IsotopicCorrection)


                            
            