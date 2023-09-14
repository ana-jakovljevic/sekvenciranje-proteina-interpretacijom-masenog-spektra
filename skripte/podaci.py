import gc, copy, os, time, random, re
import random, math, numpy as np
import matplotlib.pyplot as plt
import pandas as pd

mase_amino_kiselina = {
    "celobrojne":{
        "G": 57,
        "A": 71,
        "S": 87,
        "P": 97,
        "V": 99,
        "T": 101,
        "C": 103,
        "I": 113,
        "L": 113,
        "N": 114, 
        "D": 115,
        "K": 128,
        "Q": 128,
        "E": 129,
        "M": 131,
        "H": 137,
        "F": 147,
        "R": 156,
        "Y": 163,
        "W": 186
    },
    
    "monoizotopske": {
        'G': 57.02146,
        'A': 71.03711,
        'S': 87.03203,
        'P': 97.05276,
        'V': 99.06841,
        'T': 101.04768,
        'C': 103.00919,
        'I': 113.08406,
        'L': 113.08406,
        'N': 114.04293,
        'D': 115.02694,
        'Q': 128.05858,
        'K': 128.09496,
        'E': 129.04259,
        'M': 131.04049,
        'H': 137.05891,
        'F': 147.06841,
        'R': 156.10111,
        'Y': 163.06333,
        'W': 186.07931
    }, 
    
    "srednje": {
        'G': 57.0519,
        'A': 71.0788,
        'S': 87.0782,
        'P': 97.1167,
        'V': 99.1326,
        'T': 101.1051,
        'C': 103.1388,
        'I': 113.1594,
        'L': 113.1594,
        'N': 114.1038,
        'D': 115.0886,
        'Q': 128.1307,
        'K': 128.1741,
        'E': 129.1155,
        'M': 131.1926,
        'H': 137.1411,
        'F': 147.1766,
        'R': 156.1875,
        'Y': 163.176,
        'W': 186.2132
    }
}

celobrojna_masa_u_amino_kiselinu = {
    57: "G",
    71: "A",
    87: "S",
    97: "P",
    99: "V",
    101: "T",
    103: "C",
    113: "I",
    113: "L",
    114: "N", 
    115: "D",
    128: "K",
    128: "Q",
    129: "E",
    131: "M",
    137: "H",
    147: "F",
    156: "R",
    163: "Y",
    186: "W"
}

def prikazi_recnik(recnik):
    display(pd.DataFrame.from_dict(recnik, orient='index', columns=["masa"]))
    
def prikaz_primera_masenog_spektra():
    
    with open("../podaci/t-rex", "r") as f:
        data = [x.replace("\t", " ") for x in f.read().split("\n")]

    spektar = []
    redni_broj = 0

    for line in data:        
        if line.startswith("BEGIN IONS"):
            redni_broj += 1

        if np.any([line.startswith(x) for x in ["TITLE", "PEPMASS", "CHARGE", "END", "BEGIN"]]):
            continue

        if redni_broj == 6439:
            spektar.append([float(x) for x in line.split(" ")])

        if redni_broj > 6439:
            break
    
    maseni_spektar = pd.DataFrame(spektar, columns=["m_z", "intenzitet"])

    fig, ax = plt.subplots(figsize=(50, 20))
    ax.bar(maseni_spektar["m_z"], height=maseni_spektar["intenzitet"], width=0.8,  linewidth=1)
    
    ax.set_xlabel("masa/naelektrisanje")
    ax.set_ylabel("intenzitet")
    
    # oznake_pikova = np.where(maseni_spektar["intenzitet"] >= 50, maseni_spektar["m_z"], 0)
    #for indeks, oznaka in enumerate(list(oznake_pikova)):
    #    if oznaka > 100:
    #        ax.annotate(oznaka, (maseni_spektar.loc[indeks, "m_z"], maseni_spektar.loc[indeks, "intenzitet"]))

    plt.show()