__author__ = "Marina Igual"

""" Program to select the delta_masses pairs that are co-occurring not by chance.
    Expected co-occurring values are calculated from a random dataset which is created also in this program"""

############### TRYING TO CREATE THE RANDOM DATA ###############

import csv
import scipy.stats as stats
from itertools import combinations

############## FUNCTIONS ####################

def return_smallest(list_1, list_2):
    if len(list_1) > len(list_2):
        return(len(list_2))
    else:
        return(len(list_1))


def comb_count(mod_pept_dict):
    mass_list = []

    for seq in mod_pept_dict:
        for res in mod_pept_dict[seq]:
            for mass in list(set(mod_pept_dict[seq][res])):
                mass_list.append(mass)
            
    comb_dict = {}
    
    comb_calc = list(combinations(list(set(mass_list)), 2))

    for comb in comb_calc:
        comb_dict[comb] = 0
    
    for seq in mod_pept_dict:
        for pos in mod_pept_dict[seq]:
            masses_comb = list(combinations(list(set(mod_pept_dict[seq][pos])), 2))
            for comb in masses_comb:
                comb2 = (comb[1], comb[0])
                try:
                    a=comb_calc.index(comb) 
                    comb_dict[comb] += 1
                except:
                    try:
                        b=comb_calc.index(comb2)
                        comb_dict[comb2] += 1
                    except:
                        continue
    return(comb_dict)

def count_freq(dict_i, mass1mass2):
    m1m2 = 0
    nc_m2 = 0   #Frequency co-ocurrence without mass2
    nc_m1 = 0   #Frequency co-ocurrence without mass1
    nc_m1m2 = 0 #Frequency co-ocurrence without mass2 neither mass1
    mass1, mass2 = mass1mass2
    
    for seq in dict_i:
        for pos in dict_i[seq]:
            if mass1 in dict_i[seq][pos] and mass2 in dict_i[seq][pos]:
                m1m2 += 1
            elif mass1 in dict_i[seq][pos] and mass2 not in dict_i[seq][pos]:
                nc_m2 += 1
            elif mass2 in dict_i[seq][pos] and mass1 not in dict_i[seq][pos]:
                nc_m1 += 1
            elif mass1 not in dict_i[seq][pos] and mass1 not in dict_i[seq][pos]:
                nc_m1m2 += 1
    return(m1m2, nc_m2, nc_m1, nc_m1m2)

def list_masses(s_peptfreq, pept_mod, pept1, i):
    for n in range(i, len(s_peptfreq)):
        pept2 = s_peptfreq[n][0]
            
        masses_interes1 = [] ## Masses in pept2 but not in pept1
        masses_interes2 = [] ## Masses in pept1 but not in pept2
            
        for mass in pept_mod[pept2]:
            if mass not in pept_mod[pept1]:
                masses_interes1.append(mass)
            
        for mass in pept_mod[pept1]:
            if mass not in pept_mod[pept2]:
                masses_interes2.append(mass)
    return(pept_mod)

def random_data(pept_freq, pept_mod):        
    s_peptfreq = sorted(pept_freq.items(), key=lambda kv: kv[1], reverse = True)    ### Dictionary sorted by frequency
    
    for z in range(10):
        for i in range(len(s_peptfreq)-1):
            pept1 = s_peptfreq[i][0]
            for n in range(0, len(s_peptfreq)):
                pept2 = s_peptfreq[n][0]
                
                masses_interes1 = [] ## Masses in pept2 but not in pept1
                masses_interes2 = [] ## Masses in pept1 but not in pept2
                
                for mass in pept_mod[pept2]:
                    if mass not in pept_mod[pept1]:
                        masses_interes1.append(mass)
                
                for mass in pept_mod[pept1]:
                    if mass not in pept_mod[pept2]:
                        masses_interes2.append(mass)
                
                max_value = return_smallest(masses_interes1, masses_interes2)
                
                for x in range(max_value):
                    mass1 = masses_interes1[x]
                    mass2 = masses_interes2[x]
                    pept_mod[pept1].remove(mass2)
                    pept_mod[pept2].append(mass2)
                    pept_mod[pept2].remove(mass1)
                    pept_mod[pept1].append(mass1)
    
    return(pept_mod)

def peptmod_dict(pept_mod):
    pept_mod_ran={}
    for peptide in pept_mod:
        elements = peptide.split("_")
        try:
            pept_mod_ran[elements[0]]
            pept_mod_ran[elements[0]][elements[1]]=pept_mod[peptide]
        except:
            pept_mod_ran[elements[0]]={}
            pept_mod_ran[elements[0]][elements[1]]=pept_mod[peptide]
    return(pept_mod_ran)


def compare_values(combination):
    mass1, mass2 = combination
    mass_i = float(mass1)
    mass_n = float(mass2)
    diff = abs(mass_i - mass_n)
    if diff > 1.005:
        return(1)
    else:
        return(0)
        
def clean(comb):
    comb1, comb2 = comb
    chars = "(')"
    for char in chars:
        comb1 = comb1.replace(char, "")
        comb2 = comb2.replace(char, "")
    return(comb1, comb2)

def stat_test(o1,o2, o3, o4, e1, e2, e3, e4):
    if e1 >10 and e2 > 10 and e3> 10 and e4 > 10:
        chi_squared_stat = (((o1-e1)**2)/e1)+(((o2-e2)**2)/e2)+(((o3-e3)**2)/e3)+(((o4-e4)**2)/e4)
        p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,  df=1)
        return("Xi", p_value)
    else:
        oddsratio, p_value = stats.fisher_exact([[o1, o2], [o3, o4]])
        return("FET", p_value)
    
#################### MAIN ##########################

def CO_Selection_all(observed_forlder, out_folder):
    observed_forlder = observed_forlder + "\\"
    for RESIDUE in ["W", "Y", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "V", "U"]:
        #"A", "R", "N", "D",
        print(RESIDUE)
        
        observed_file = observed_forlder + RESIDUE + "_PosDM_Observed.csv"
        
        file_in = open(observed_file, "r")
        line1 = next(file_in)   #Am frequencies
        line2 = next(file_in)   #Am list
        
        line2 = line2.strip()
        masses = line2.split(",")
        masses = masses[2:]
        
        pept_freq = {}
        pept_mod = {}
        
        for line in file_in:
            if line != "\n":
                line = line.strip()
                splits2 = line.split(",")
                freq = splits2[0].strip()
                pept_res = splits2[1].strip()
                pept_mod[pept_res]=[]
                tf = splits2[2:]        #Values with 0 and 1 for DeltaMasses in the peptide
                pept_freq[pept_res]=int(freq)
                for i in range(len(tf)):
                    if tf[i] == "1":
                        pept_mod[pept_res].append(masses[i])                
        file_in.close()    
        #### PART 1. CREATING RANDOM DATA  
        
        
        
        ## DictionarY with observed data
        
        pept_mod_obs = peptmod_dict(pept_mod)
        observed_comb = comb_count(pept_mod_obs)
        
        #Create a list with the comb of the final file
        
        comb_final = []
        for comb in observed_comb:
            if comb[0] != comb[1] and observed_comb[comb] != 0 and compare_values(comb) == 1:
                comb_final.append(comb)
        
        observed_mean = {}       
        for comb in comb_final:
                fm1m2, fm1, fm2, fnm = count_freq(pept_mod_obs, comb)
                observed_mean[comb]=[fm1m2, fm1, fm2, fnm]
             
        ### Random Data
                
        expected_all = {}
        for comb in comb_final:
            expected_all[comb]=[]
              
        for i in range(100):
            peptmod_mod = random_data(pept_freq, pept_mod)
            peptmod_ran = peptmod_dict(peptmod_mod)
            for comb in comb_final:
                fm1m2_r, fm1_r, fm2_r, fnm_r = count_freq(peptmod_ran, comb)
                expected_all[comb].append([fm1m2_r, fm1_r, fm2_r, fnm_r])
        
        expected_mean = {}
        for comb in expected_all:
            val_m1m2 = 0
            val_m1Nm2 = 0
            val_m2Nm1 = 0
            val_Nm1Nm2 = 0
            num = 0
            for list_i in expected_all[comb]:
                val_m1m2 += list_i[0]
                val_m1Nm2 += list_i[1]
                val_m2Nm1 += list_i[2]
                val_Nm1Nm2 += list_i[3]
                num += 1
            expected_mean[comb] = [val_m1m2/num, val_m1Nm2/num, val_m2Nm1/num, val_Nm1Nm2/num]
            
        #### PRINTING RESULTS
        out_folder = out_folder + "\\"
        name_file = out_folder + RESIDUE + "_CO.csv"
        file_out = open(name_file, 'w', encoding="ISO-8859-1", newline='')
        wr = csv.writer(file_out)
        line = ["m1","m2", "m1m2_Obs","m1Nm2_Obs","m2Nm1_Obs", "Nm1Nm2_Obs", "m1m2_Exp(10c)","m1Nm2_Exp(10c)","m2Nm1_Exp(10c)", "Nm1Nm2_Exp(10c)", "Test(Xi_square/FET)", "p_value"]
        wr.writerow(line)
        for comb in comb_final:
            c1, c2 = clean(comb)
            c3 = observed_mean[comb][0]
            c4 = observed_mean[comb][1]
            c5 = observed_mean[comb][2]
            c6 = observed_mean[comb][3]
            c7 = expected_mean[comb][0]
            c8 = expected_mean[comb][1]
            c9 = expected_mean[comb][2]
            c10 = expected_mean[comb][3]
            c11, c12 = stat_test(c3,c4,c5,c6,c7,c8,c9,c10)
            
            line = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12]
            wr.writerow(line)      
                
        file_out.close()
        
        ### PRINTING ONE OF THE RANDOM MODELS
        list_masses=[]
        for seq in peptmod_ran:
            for pos in peptmod_ran[seq]:
                for mass in list(set(peptmod_ran[seq][pos])):
                    list_masses.append(mass)
        
        list_masses = list(set(list_masses))
        
        name_file = out_folder + RESIDUE + "_PeptDM_Random.csv"
        file_out = open(name_file, 'w', encoding="ISO-8859-1", newline='')
        wr = csv.writer(file_out)
            
        val_row1 = ["Peptide_numberResidue"]
            
        for mass in list_masses:
            val_row1.append(mass)
            
        wr.writerow(val_row1)
                
        for seq in peptmod_ran:
            for pos in peptmod_ran[seq]:
                seq_id = seq + "_" + str(pos)
                val_row = [seq_id]
                for i in range(1, len(val_row1)):
                    mass = val_row1[i]
                    try:
                        b = peptmod_ran[seq][pos].index(mass)
                        val_row.append(1)
                    except:
                        val_row.append(0)
                wr.writerow(val_row)
        file_out.close()

################# MAIN ####################

observed_master = r"S:\U_Proteomica\LABS\LAB_JMR\Marfan\Human\plasma\HFviejo\Marfan_plasma_human_Alvaro\Plasma_HUman-Marfan\OS-PTM\SHIFTS\Cooccurrences\Observed"   ##Master folder with observed data

out = r"S:\U_Proteomica\LABS\LAB_JMR\Marfan\Human\plasma\HFviejo\Marfan_plasma_human_Alvaro\Plasma_HUman-Marfan\OS-PTM\SHIFTS\Cooccurrences\Selected"     ##Master folder for the output files
    
CO_Selection_all(observed_master, out)
