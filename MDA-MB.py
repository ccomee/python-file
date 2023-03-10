from cobra import Model, Reaction, Metabolite
from cobra.sampling import sample
from PIL import Image, ImageOps
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

metabolite_excel = pd.read_excel("C:/Users/courauc/Desktop/python file/data/mammalian_cell_model_final.xlsx", sheet_name=0)
reaction_excel = pd.read_excel("C:/Users/courauc/Desktop/python file/data/mammalian_cell_model_final.xlsx", sheet_name=1)
YSI_excel = pd.read_excel("C:/Users/courauc/Desktop/python file/data/YSI.xlsx")
seahorse_excel = pd.read_excel("C:/Users/courauc/Desktop/python file/data/Seahorse_atp.xlsx")

#%% Create model

model = Model("MDA-MB")

model.add_metabolites([Metabolite(metabolite_excel.at[i,'ID'], name=metabolite_excel.at[i,'Name']) for i in range(metabolite_excel.shape[0])])


for j in range(reaction_excel.shape[0]):
    consumed_str, produce_str = reaction_excel.at[j, 'Formula'].replace(' ', '').split('-->')
    if  consumed_str[-1] == '<' :
        consumed_str = consumed_str.replace('<', '')
        reversible = True
    else :
        reversible = False 
    
    
    consumed_str = consumed_str.split('+')
    produce_str = produce_str.split('+')
    produce_str = list(filter(None, produce_str))
    
    reaction = Reaction(reaction_excel.at[j, 'ID'], name = reaction_excel.at[j, 'Name'])
    reaction.lower_bound = -1000*reversible
    
    for metabolite_str in consumed_str :
        nb_str, id_str = metabolite_str.split(')')
        reaction.add_metabolites({model.metabolites.get_by_id(id_str):-float(nb_str[-3:])})

    for metabolite_str in produce_str :
        nb_str, id_str = metabolite_str.split(')')
        reaction.add_metabolites({model.metabolites.get_by_id(id_str):float(nb_str[-3:])})
    
    model.add_reactions([reaction])


#%% Image tools

def get_concat_h(im1, im2):
    if im1 == None :
        return im2
    dst = Image.new('RGB', (im1.width + im2.width, im1.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst

def get_concat_v(im1, im2):
    if im1 == None :
        return im2
    dst = Image.new('RGB', (im1.width, im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

def rogner(im):
    width, height = im.size
    x = width/9
    y = width*10/11
    return im.crop((x, 0, y, height))

#%% YSI & Seahorse data used for Flux sampling

name = 'Flux distribution'
border_size = 10
colors = ['lightcoral', 'lightskyblue', 'lightgreen']

solutions = [] 
for i in range(3):
    mean_lac_c, std_lac_c = np.mean(YSI_excel['Lactate'][0+4*i:3+4*i]),np.std(YSI_excel['Lactate'][0+4*i:3+4*i])
    model.reactions.ex_lac.bounds = mean_lac_c - std_lac_c, mean_lac_c + std_lac_c

    mean_glc_c, std_glc_c = np.mean(YSI_excel['Glucose'][0+4*i:3+4*i]),np.std(YSI_excel['Glucose'][0+4*i:3+4*i])
    model.reactions.ex_glc.bounds = mean_glc_c - std_glc_c, mean_glc_c + std_glc_c

    mean_gln_c, std_gln_c = np.mean(YSI_excel['Glutamine'][0+4*i:3+4*i]),np.std(YSI_excel['Glutamine'][0+4*i:3+4*i])
    model.reactions.ex_gln.bounds = mean_gln_c - std_gln_c, mean_gln_c + std_gln_c

    mean_glu_c, std_glu_c = np.mean(YSI_excel['Glutamate'][0+4*i:3+4*i]),np.std(YSI_excel['Glutamate'][0+4*i:3+4*i])
    model.reactions.ex_glu.bounds = mean_glu_c - std_glu_c, mean_glu_c + std_glu_c
    
    mean_o2_c, std_o2_c = np.mean(seahorse_excel["ATP OCR (pmol/cell/min)"][0+8*i:3+8*i]),np.std(seahorse_excel["ATP OCR (pmol/cell/min)"][0+8*i:3+8*i])
    model.reactions.ex_o2.bounds = -mean_o2_c - std_o2_c, -mean_o2_c + std_o2_c
    
    solutions.append(sample(model, 1000))

titles = ['Parental', 'BrM2', 'LM2']
reaction_IDs = [reaction_excel.at[i,'ID'] for i in range(reaction_excel.shape[0])]

range_min = np.min([[np.min(solutions[k][reaction_IDs[l]]) for k in range(3)] for l in range(24)])
range_max = np.max([[np.max(solutions[k][reaction_IDs[l]]) for k in range(3)] for l in range(24)])

for j in range(1,25):
    plt.figure(figsize=(15,5), dpi = 150, edgecolor='black')
    #range_min = np.min([np.min(solutions[k][reaction_IDs[j-1]]) for k in range(3)])
    #range_max = np.max([np.max(solutions[k][reaction_IDs[j-1]]) for k in range(3)])
    for i in range(1,4):
        plt.subplot(1,3,i)
        plt.hist(solutions[i-1][reaction_IDs[j-1]],range=(range_min,range_max), rwidth=0.7, color=colors[i-1], ec = 'k', lw = 0.8)
        plt.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
        plt.title(titles[i-1], wrap = True)
        plt.xlabel('Flux [pmol/cell/hr]', style = 'italic')
    plt.suptitle(model.reactions.get_by_id(reaction_IDs[j-1]).name.capitalize(), fontsize = 15)
    plt.savefig('C:/Users/courauc/Desktop/python file/images/misc/{}.png'.format(j))
    plt.close()

im_total = None
for j in range(6):
    im = None
    for i in range(4):
        im = get_concat_h(im, ImageOps.expand(rogner(Image.open('C:/Users/courauc/Desktop/python file/images/misc/{}.png'.format(i*6+j+1))), border=border_size//2, fill='black'))
    im_total = get_concat_v(im_total, im)

im_total.save('C:/Users/courauc/Desktop/python file/images/{}.png'.format(name))