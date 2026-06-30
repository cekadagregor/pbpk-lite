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
    fig.supxlabel('Time (min)')
    fig.supylabel('Concentration (ng/mL)')
    plt.savefig(save_as, dpi=300, bbox_inches="tight")

def graph_venous_helper(time, concentrations, save_as, limit_of_detection, log):
    plt.rcParams['font.size'] = 8
    fig, axs = plt.subplots(1, 1, layout='constrained', sharex=True, figsize=(6.27, 3.5), dpi=300)
    axs.plot(time, concentrations[15, :])
    if limit_of_detection is not None:
        axs.axhline(y=limit_of_detection, color='r', linestyle='--')
    axs.set_facecolor('#F2F2F2')
    axs.set_xlabel('Time (min)')
    axs.set_ylabel('Concentration (ng/mL)')
    if log == True:
        axs.set_yscale('log')
    plt.savefig(save_as, dpi=300, bbox_inches="tight")

def graph_compartments_helper(time, concentrations, compartments, save_as):
    plt.rcParams['font.size'] = 8
    fig, axs = plt.subplots(1, 1, layout='constrained', sharex=True, figsize=(6.27, 3.5), dpi=300)
    for compartment in compartments:
        axs.plot(time, concentrations[identify_compartment(compartment), :], label=compartment)
    axs.set_facecolor('#F2F2F2')
    axs.set_xlabel('Time (min)')
    axs.set_ylabel('Concentration (ng/mL)')
    axs.legend()
    plt.savefig(save_as, dpi=300, bbox_inches="tight")