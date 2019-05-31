__author__ = "Marina Igual"

"""Program to prepare MS data in order to perform different proteomics analysis"""

############### IMPORT #############
import sys
sys.path.append(r"C:\Probes_CO\Code\Data_Preparator")
import MasterProtein
import Table_DeltaP2Q

############## MAIN ################ 

if __name__ == "__main__":
    
    ############ MASTERPROTEIN###############
    
    #MasterFile => Write one or more path to the folder where are located all the sub-folders of SHIFTS 
    masterFile = [r"S:\U_Proteomica\PROYECTOS\PESA_omicas\Comet-PTM-2a-5ta_Cohortes_V1",
                  r"S:\U_Proteomica\PROYECTOS\PESA_omicas\Comet-PTM-2a-5ta_Cohortes_V2"]
    #fileName => Name of the SHIFT output file. All files must have the same name
    fileName = "IsotopCorrection_TargetData_withSequence-massTag.txt"
    #Name of the folders where the SHIFTS files are located
    subfolders = ["cXcorr_Len_Rank_Results", "Vertex"]
    list_master = MasterProtein.MasProt(masterFile, fileName, subfolders)
    
    ########### DELTAP2Q ############
    
    fileOut = r"S:\U_Proteomica\UNIDAD\DatosCrudos\laboratorio\Marina_Igual\HIPERMOD\p2site\PESA_V1V2_deltap2q_CorrMass.txt"
    IsotopicCorrection = "Yes"
    Table_DeltaP2Q.Relation_File(list_master, fileOut, IsotopicCorrection)
    
    ### NEXT STEP => DATA MUST UNDERGO P2SITE PROGRAM