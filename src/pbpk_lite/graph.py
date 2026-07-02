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


def _resolve_time_axis(time, time_unit='min'):
    time = np.asarray(time)
    if time_unit == 'min':
        return time, 'Time (min)'
    if time_unit == 'hours':
        return time / 60, 'Time (hours)'
    if time_unit == 'days':
        return time / (60 * 24), 'Time (days)'
    raise ValueError("time_unit must be one of 'min', 'hours', or 'days'")


def graph_whole_helper(time, concentrations, save_as, time_unit='min'):
    plt.rcParams['font.size'] = 8
    fig, axs = plt.subplots(4, 4, layout='constrained', sharex=True, figsize=(6.27, 3.5), dpi=300)
    time_values, time_label = _resolve_time_axis(time, time_unit)
    comp = 0
    for i in matrix_index():
        axs[i].plot(time_values, concentrations[comp, :])
        axs[i].set_title(identify_compartment(comp))
        axs[i].set_facecolor('#F2F2F2')
        comp += 1
    fig.supxlabel(time_label)
    fig.supylabel('Concentration (ng/mL)')
    plt.savefig(save_as, dpi=300, bbox_inches="tight")
    plt.close(fig)

def graph_venous_helper(time, concentrations, save_as, limit_of_detection, log, time_unit='min'):
    plt.rcParams['font.size'] = 8
    fig, axs = plt.subplots(1, 1, layout='constrained', sharex=True, figsize=(6.27, 3.5), dpi=300)
    time_values, time_label = _resolve_time_axis(time, time_unit)
    axs.plot(time_values, concentrations[15, :])
    if limit_of_detection is not None:
        axs.axhline(y=limit_of_detection, color='r', linestyle='--')
    axs.set_facecolor('#F2F2F2')
    axs.set_xlabel(time_label)
    axs.set_ylabel('Concentration (ng/mL)')
    if log == True:
        axs.set_yscale('log')
    plt.savefig(save_as, dpi=300, bbox_inches="tight")
    plt.close(fig)


def graph_compartments_helper(time, concentrations, compartments, save_as, time_unit='min'):
    plt.rcParams['font.size'] = 8
    fig, axs = plt.subplots(1, 1, layout='constrained', sharex=True, figsize=(6.27, 3.5), dpi=300)
    time_values, time_label = _resolve_time_axis(time, time_unit)
    for compartment in compartments:
        axs.plot(time_values, concentrations[identify_compartment(compartment), :], label=compartment)
    axs.set_facecolor('#F2F2F2')
    axs.set_xlabel(time_label)
    axs.set_ylabel('Concentration (ng/mL)')
    axs.legend()
    plt.savefig(save_as, dpi=300, bbox_inches="tight")
    plt.close(fig)