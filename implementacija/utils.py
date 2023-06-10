from IPython.display import display, HTML
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

display(HTML("<style>.container { width:90% !important; }</style>"))

amino_kiseline = {
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
}


amino_kiseline_monoizotopske = {
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
}

amino_kiseline_srednje = {
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

celobrojne_mase_amino_kiselina = {
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

def masa_amino_kiseline(amino_kiselina):
    return amino_kiseline[amino_kiselina]

def amino_kiselina_mase(masa):
    return mase_amino_kiselina[masa]