import numpy as np

def Kauffmann(x):
    return 0.61/(x-0.05)+1.3

def Kewley(x):
    return 0.61/(x-0.47)+1.19

def Stasinska(x):
    return (-30.787+1.1358*x+0.27297*x**2)*np.tanh(5.7409*x)-31.093

def Schawinski(x):
    return 1.05*x+0.45

def main_AGN_BPT2(x):
    return 0.73/(x+0.59)+1.33

def LINER_SY2_BPT2(x):
    return 1.18*x+1.30

def main_AGN_BPT3(x):
    return 0.72/(x-0.32)+1.30

def LINER_SY2_BPT3(x):
    return 1.89*x+0.76