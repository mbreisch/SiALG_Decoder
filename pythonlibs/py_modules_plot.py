import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("./")
import py_modules_fits as py_fit

def PlotCalibrationCurve(volt,inputlists,idlist):
    [spe_list,spe_err_list,fitlist,valuelist] = inputlists

    fig = plt.figure(0,figsize=[10,10], facecolor='white')
    plt.title("SiPM Pulseheight vs Bias Voltage")

    colorlist = ["red","peru","orange","olivedrab","darkgreen","aqua"]

    for enum,id in enumerate(idlist):    
        plt.errorbar(volt,spe_list[enum],xerr=0.094,yerr=np.array(spe_err_list[enum])/2,color=colorlist[enum],marker='x',linestyle='',label="SiPM-ID:{}".format(id))
        a,b = fitlist[enum]
        fl = py_fit.linear(volt,a,b)
        plt.plot(volt,fl,c=colorlist[enum],linestyle='-',label="SiPM-ID:{},U={:.3f}V".format(id,valuelist[enum]))

    plt.xlabel("Bias Voltage [V]")
    plt.ylabel("Detected Pulse Height [ADC]")
    plt.legend()