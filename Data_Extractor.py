__author__ = "Marina Igual"

import re
import os
import sys
import math
import csv

sys.path.append(r"C:\Probes_CO\Code")
import RT_1_PosDM_PSM

def truncate(f, n):
    return math.floor(f * 10 ** n) / 10 ** n    
    
def check_mass(mass_i):
    mass_i = float(mass_i)
    mass = truncate(mass_i, 1)
    value = 1
    if -0.8 <= mass <= 0.8:                                 #Considering only modifications
        value = 0
    elif abs(mass) == 1.0 or abs(mass) == 2.0:
        value = 0
    elif -1.3 <= mass <= -1.0:
        value = 0
    else:
        value = 1
    return(value)

#Function to remove artefacts
def remove_Nter_art(mass, position, length):
    crit_value = ((length/4)+1)
    mass_num = float(mass)
    mass_i = truncate(mass_num, 2)
    values_del = [-229.16, -229.19, -229.10, -229.07, -230.18, -230.14, -203.12, -203.18, -186.15, -186.1, -186.06, -172.14, -172.08, -172.05, -172.04, -157.17, -157.14, -157.18, -147.15, -147.16, -147.06, -147.07, -129.14, -129.15, -129.04, -114.17, -114.13, -114.13, -114.14, -114.04,-72.12, -72.02, -72.13, 100.01, 100.02] 
    if mass_i in values_del and position < crit_value:  #Search if the mass is in the 1st Quartil of the sequence
        return(0)
    else:
        return(1)

#Function to remove partial digestions
def check_pardig(list_seqs):
    list_noParGid = []
    list_s = sorted(list_seqs)
    for i in range(len(list_s)):
        result = [x for x in list_s[i+1:] if x.startswith(list_s[i])]
        if result == []:
            list_noParGid.append(list_s[i])
    return(list_noParGid)

def change_dec(num):
    decimals = re.search("\.+\d+", num)
    if decimals:
        if len(decimals[0]) < 7:
            n_0 = '0'*(7 - len(decimals[0]))
            num = num + n_0
    return (num)

def short_lines(list_lines):
    new_list=[]
    i=len(list_lines[0])
    for line in list_lines:
        new_list.append(line[0:i+1])
        i = i-1
    return(new_list)

def change_num(list_str):
    list_new = []
    list_def = []
    for string in list_str:
        num = float(string)
        list_new.append(num)
    list_sort = sorted(list_new)
    for float_num in list_sort:
        num_str = change_dec(str(float_num))
        list_def.append(num_str)
    return(list_def)

def DataExtr(file_p2q, folderOut):
    ##### PART 1. EXTRACT THE SPECIFIC MASSES IN POSITION 0
    
    for RESIDUE in ["W", "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "Y", "V", "U"]:
        print(RESIDUE)          
        ### Function to obtain the frequencies of Am in position -3,-2, -1, 0, 1, 2, 3
        
        ResMassPos = {}
        list_AA = ["W", "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "Y", "V", "U"]
        for AA in list_AA:
            ResMassPos[AA]={}
        
        for filename in [file_p2q]:
            with open(filename) as f2:
                next(f2)
                for line2 in f2:
                    if line2 != "\n":
                        
                        splits2 = line2.split("\t")
                        prot_id = splits2[0].strip()
                        deltaPEP = splits2[1].strip()
                        tags = splits2[2].strip()
                        seq = re.sub('\[\S+\]', '', deltaPEP)
            
                        start = deltaPEP.index("[")
                        end = deltaPEP.index("]")
                        fisrtHalf = deltaPEP[0:start]
                        secondHalf = deltaPEP[end+1:]
                        massMod = deltaPEP[(start + 1): (end)]
                        
                        if tags == "0" and check_mass(massMod) == 1:                    
                            i = start - 1
                            for n in range(len(fisrtHalf)):
                                res = fisrtHalf[n]
                                try:
                                    ResMassPos[res][massMod][i]
                                    ResMassPos[res][massMod][i] += 1
                                except:
                                    try:
                                        ResMassPos[res][massMod]
                                        ResMassPos[res][massMod][i] = 1
                                    except:
                                        ResMassPos[res]
                                        ResMassPos[res][massMod] = {}
                                        ResMassPos[res][massMod][i] = 1 
                                i = i - 1
        
                            e = - 1
                            for n in range(0, len(secondHalf)):      
                                res = secondHalf[n]
                                try:
                                    ResMassPos[res][massMod][e]
                                    ResMassPos[res][massMod][e] += 1
                                except:
                                    try:
                                        ResMassPos[res][massMod]
                                        ResMassPos[res][massMod][e] = 1
                                    except:
                                        ResMassPos[res]
                                        ResMassPos[res][massMod] = {}
                                        ResMassPos[res][massMod][e] = 1
                                e = e - 1
        
        #Create a list with all the masses specific for position 0 in W:
        mass_0 = []
        
        for mass in ResMassPos[RESIDUE]:
            total_inc = 0
 
            for pos1 in [-3, -2, -1, 0, 1, 2, 3]:
                try:
                    ResMassPos[RESIDUE][mass][pos1]
                    total_inc = total_inc + ResMassPos[RESIDUE][mass][pos1]
                except:
                    total_inc = total_inc + 0
            total = total_inc/7
            try:
                ResMassPos[RESIDUE][mass][0]
                if total < ResMassPos[RESIDUE][mass][0]:
                    mass_0.append(mass)
                else:
                    continue
            except:
                continue
        
        print("Number of masses selected: ", len(mass_0)) 
        #479 with position -3, -2, -1, 1, 2, 3
        
        ##### PART 2. EXTRACT THE SPECIFIC CO-OCURRENCES 
        
        
        mod_pept = {}       #Dictionary with all proteins, its modified residues and the sequences. 
        
        with open(file_p2q, 'r') as file_in:
            next(file_in)
            for line2 in file_in:
                if line2 != "\n":
                    splits2 = line2.split("\t")
                    prot_id = splits2[0].strip()
                    deltaPEP = splits2[1].strip()
                    tags = splits2[2].strip()
                    seq = re.sub('\[\S+\]', '', deltaPEP)
                    
                    position1 = deltaPEP.index("[")              #Modification position 0
                    position2 = deltaPEP.index("]")
                    massMod = deltaPEP[(position1 + 1): (position2)]
                    
                    if remove_Nter_art(massMod, position1 - 1, len(seq)) == 1:
                        try:
                            b = mass_0.index(massMod)
                            if tags == '0' and check_mass(massMod) == 1 and seq[position1 - 1] == RESIDUE:
                                try:
                                    mod_pept[seq][position1 - 1]
                                    mod_pept[seq][position1 - 1].append(massMod)
                                except:
                                    try:
                                        mod_pept[seq]
                                        mod_pept[seq][position1 - 1]=[massMod]
                                    except:
                                        mod_pept[seq]={}
                                        mod_pept[seq][position1 - 1]=[massMod]
                        except:
                            continue
                                        
        print("Number of peptides selected: ",len(mod_pept))    #There are 832 different peptides
        
        ############## SUM FREQ WITH A DIFFERENCE OF C13 #############
        mod_pept_uniq = {}
        
        for seq in mod_pept:
            mod_pept_uniq[seq] = {}
            for pos in mod_pept[seq]:
                masses = list(set(mod_pept[seq][pos]))
                mod_pept_uniq[seq][pos] = masses
                                
        for seq in mod_pept_uniq:
            for pos in mod_pept_uniq[seq]:
                
                masses = change_num(list(set(mod_pept_uniq[seq][pos])))
                masses.append("0")
                mod_pept_uniq[seq][pos] = []
                list_diff = []
                mod_pept_uniq[seq][pos].append(masses[0])
                if len(masses) == 2:
                    continue
                else:
                    for i in range(len(masses)-1):
                        diff = abs(float(masses[i]) - float(masses[i+1]))
                        if 1.00 <= diff <= 1.01:
                            continue
                        else:
                            if masses[i+1] != "0":
                                mod_pept_uniq[seq][pos].append(masses[i+1])    
        
        
        list_masses=[]
        for seq in mod_pept_uniq:
            for pos in mod_pept_uniq[seq]:
                for mass in list(set(mod_pept_uniq[seq][pos])):
                    list_masses.append(mass)
        
        list_masses = list(set(list_masses))
        
        dict_masses = {}    #Masses frequency
        dict_pept_freq = {} #Peptide Frequency
        for mass in list_masses:
            dict_masses[mass] = 0
        
        for seq in mod_pept_uniq:
            for pos in mod_pept_uniq[seq]:
                try:
                    dict_pept_freq[seq]
                    dict_pept_freq[seq][pos]=len(mod_pept_uniq[seq][pos])
                except:
                    dict_pept_freq[seq]={}
                    dict_pept_freq[seq][pos]=len(mod_pept_uniq[seq][pos])
                
                for mass in mod_pept_uniq[seq][pos]:
                    dict_masses[mass] += 1
        ####################### TABLE PEPT_Am ##################
        
        folderOut = folderOut + "\\"
        name_file = folderOut + RESIDUE + "_PosDM_Observed.csv"
        file_out = open(name_file, 'w', encoding="ISO-8859-1", newline='')
        wr = csv.writer(file_out)
        
        val_row1 = ["","Sum"]
        val_row2 = ["","Peptide_numberResidue"]
            
        for mass in list_masses:
            val_row1.append(dict_masses[mass])
            val_row2.append(mass)
            
        wr.writerow(val_row1)
        wr.writerow(val_row2)
                
        for seq in mod_pept_uniq:
            for pos in mod_pept_uniq[seq]:
                seq_id = seq + "_" + str(pos)
                val_row = [dict_pept_freq[seq][pos],seq_id]
                for i in range(2, len(val_row2)):
                    mass = val_row2[i]
                    try:
                        b = mod_pept_uniq[seq][pos].index(mass)
                        val_row.append(1)
                    except:
                        val_row.append(0)
                wr.writerow(val_row)
        file_out.close()
        
############ MAIN ##############

if __name__ == "__main__":

    ###### PosDM OBSERVED ##########
    #file_in => Enter input file. It must be a deltap2q file, generated previously
    file_in = r"S:\U_Proteomica\LABS\LAB_JMR\Marfan\Human\plasma\HFviejo\Marfan_plasma_human_Alvaro\Plasma_HUman-Marfan\OS-PTM\SHIFTS\Cooccurrences\data_input\p2q_rpos_v2.tsv" 
    folder_out = r"S:\U_Proteomica\LABS\LAB_JMR\Marfan\Human\plasma\HFviejo\Marfan_plasma_human_Alvaro\Plasma_HUman-Marfan\OS-PTM\SHIFTS\Cooccurrences\Observed"
    DataExtr(file_in, folder_out)
    
    ###### PosDM PSM ##########
    RT_1_PosDM_PSM.PosDM(file_in, folder_out)    
    