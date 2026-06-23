import matplotlib.pyplot as plt
import numpy as np

compartment_names = ['adipose', 'bone', 'brain' , 'gut', 'heart', 'kidney', 'liver', 'lung', 'muscle', 'pancreas', 'skin', 'spleen', 'stomach', 'testes', 'arterial_blood', 'venous_blood']

def identify_compartment(arg):
    if type(arg) == int:
        return compartment_names[arg]
    else:
        return compartment_names.index(arg)

def matrix_index():
    for i in range(4):
        for j in range(4):
            yield i, j

def graph_whole_helper(time, concentrations, save_as):
    plt.rcParams['font.size'] = 8
    fig, axs = plt.subplots(4, 4, layout='constrained', sharex=True, figsize=(6.27, 3.5), dpi=300)
    comp = 0
    for i in matrix_index():
        axs[i].plot(time, concentrations[comp, :])
        axs[i].set_title(identify_compartment(comp))
        axs[i].set_facecolor('#F2F2F2')
        comp += 1
    plt.savefig(save_as, dpi=300, bbox_inches="tight")