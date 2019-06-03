__author__ = "Marina Igual"

"""Program to read all tables with significant Co-ocurrences masses 
    and create the relation table for SanXoT"""
    
    
############# IMPORTS ###############
import re
import csv
import math

############ FUNCTIONS ##############


def extract_id(fasta_id):
    prot = re.search("\|{1}\w{6}\-*\w*\|{1}", fasta_id)
    prot_id = prot[0].replace("|","").strip()
    gen = re.search("GN=\w*\s", fasta_id)
    if gen:
        gen_id = gen[0].replace("GN=","").strip()
        return(prot_id, gen_id)
    else:
        prot_name = re.search("\|{1}\w+\s{1}", fasta_id)
        name = prot_name[0].replace("|","").strip()
        if prot_name:
            return(prot_id, name)

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


def PosDM(file_p2q, folderOut):   
### Function to obtain the frequencies of Am in position -3,-2, -1, 0, 1, 2, 3
    ############################ PART 1. EXTRACT THE SPECIFIC MODIFICATION FOR W ####################################
    
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
    mass_0 = {}
    
    for RESIDUE in ResMassPos:
        mass_0[RESIDUE]=["0"]
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
                    mass_0[RESIDUE].append(mass)
                else:
                    continue
            except:
                continue
    ######################### PART 2 #####################################
    
    pos_DM = {}    #Dictionary with the peptides of each protein and its deltamases
    masses = {}
    res_list = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "U"]
    for AA in res_list:
        masses[AA]=[]
    
    with open(file_p2q, 'r') as file_in:
        next(file_in)
        for line2 in file_in:
            if line2 != "\n":
                splits2 = line2.split("\t")
                prot_fasta = splits2[0].strip()
                deltaPEP = splits2[1].strip()
                tags = splits2[2].strip()
                
                if tags == '0':
                    seq = re.sub('\[\S+\]', '', deltaPEP)       #Sequence
                    prot, name = extract_id(prot_fasta)         #Prot_ID and Gen_ID
                    id_prot = prot + "|" + name
                    id_info = prot_fasta.split("_")
                    id_residue = id_info[1] + '_' + id_info[0]  #Position (in the protein) and Residue
                    
                    position1 = deltaPEP.index("[")              #Modification position 0
                    position2 = deltaPEP.index("]")
                    massMod = deltaPEP[(position1 + 1): (position2)]
                    
                    if check_mass(massMod) == 1:
                        massMod = massMod
                    else:
                        massMod = "0"
                    
                    try:
                        b = mass_0[id_info[1]].index(massMod)
                        masses[id_info[1]].append(massMod)
                        try:
                            pos_DM[id_prot][id_residue][massMod]
                            pos_DM[id_prot][id_residue][massMod] += 1
                            
                        except:
                            try:
                                pos_DM[id_prot][id_residue]
                                pos_DM[id_prot][id_residue][massMod]=1
                            except:
                                try:
                                    pos_DM[id_prot]
                                    pos_DM[id_prot][id_residue]={}
                                    pos_DM[id_prot][id_residue][massMod]=1
                                except:
                                    pos_DM[id_prot]={}
                                    pos_DM[id_prot][id_residue]={}
                                    pos_DM[id_prot][id_residue][massMod]=1
                    except:
                        continue
    
    for res in masses:
        masses[res] = list(set(masses[res]))
    folderOut = folderOut + "\\"
    for RESIDUE in masses:
        file_name = folderOut + RESIDUE + "_PosDM.csv"
        file_out = open(file_name, 'w', encoding="ISO-8859-1", newline='')
        wr = csv.writer(file_out)
        
        line1 = ['']
        for mass in masses[RESIDUE]:
            line1.append(mass)
        wr.writerow(line1)
        
        for prot in pos_DM:
            for id_residue in pos_DM[prot]:
                res, pos = id_residue.split("_")
                
                if res == RESIDUE:
                    identifier = pos + "_" + prot
                    line = [identifier]
                    for mass in line1:
                        try:
                            pos_DM[prot][id_residue][mass]
                            line.append(pos_DM[prot][id_residue][mass])
                        except:
                            line.append(0)
                    wr.writerow(line)               