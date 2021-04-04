import numpy as np
from os import path
import os
import matplotlib.pyplot as plt
import Resonator
import generate
import logger_module
logger = logger_module.getlogger('Q_Factor')

def single_measurement():
    # path to csv file from vna
    csv_path = r"C:\Users\physicsuser\Dropbox (Weizmann Institute)\Quantum Circuits Lab\Ofir\Cavity Qubit\cooldowns\cooldown3\Read out\RR_S21_-20DdBm_Delay_corrected_55.4ns_2042021.csv"
    # creating measurement object
    # "config" argument can take the value 'T' for T connectoe configuration, and 'circulator' for circulator configuration
    # "s_mat_element" is for reading the the desried data from the csv file
    Readout = Resonator.Measurement(csv_path, config='circulator', s_mat_element='21')
    # running the calculationd, use "plot_data=True" to see the algorithm converged properly
    Ql_calc, Qc_calc, Qi_calc, fr = Readout.measure(plot_data=True)


if  __name__=="__main__":
    # multi()
    single_measurement()
plt.show()
