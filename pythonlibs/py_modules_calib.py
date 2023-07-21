import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import scipy.ndimage as filters

sys.path.append("./")
import py_modules_fits as py_fit
import py_modules_plot as py_plot


"""
Calculates the sum of multiple Gaussian functions at the given x-values.

Parameters:
    x (array-like): The x-values at which the Gaussian functions are evaluated.
    *params: Variable-length argument list containing the parameters for the Gaussian functions.
                The parameters are organized in groups of three: amplitude, mean, and standard deviation.
                The number of parameters should be a multiple of 3.

Returns:
    float: The sum of the Gaussian functions evaluated at the given x-values.
    
"""
def GetCalibrationCurve(AllDictonary,print_single_volt,id):
    pedestal = []
    spe = []
    dpe = []
    err_pedestal = []
    err_spe = []
    err_dpe = []
    volt = []

    if print_single_volt==True:
        fig, axes = plt.subplots(5, 4, figsize=(20, 20))
        plt.subplots_adjust(hspace=0.3)
    run_index = 0
    for run,data in AllDictonary.items():
        rows = run_index // 4
        cols = run_index % 4
        ldata = data["TTree_Analysis_PHD;1"]["PulseHeight_Ch2"]
        bin_width = 1/4.096
        if print_single_volt==True:
            axes[rows, cols].set_title("SiPM {} at {} V".format(id,run))
            axes[rows, cols].hist(ldata,bins=np.arange(min(ldata),max(ldata)+bin_width,bin_width),label="Data",alpha=0.5,log=False)
            axes[rows, cols].set_xlim([-50,0])
            axes[rows, cols].set_xlabel("Pulse Height [ADC]")
            axes[rows, cols].set_ylabel("Counts")

        hist, bins = np.histogram(ldata,bins=np.arange(min(ldata),max(ldata)+bin_width,bin_width))
        x_data = np.array((bins[:-1] + bins[1:]) / 2)
        y_data = np.array(hist)

        # Step 2: Smooth the histogram using a Gaussian kernel
        smoothed_histogram = filters.gaussian_filter1d(hist, sigma=1)

        # Step 3: Identify local maxima
        p,_ = find_peaks(smoothed_histogram,distance=20,height=50)

        printflag = "N"
        if len(p)==1:
            initial_guess = []
            bounds_low = []
            bounds_up = []

            initial_guess.extend([-5,1,1000])
            initial_guess.extend([-6,1,500])

            bounds_low.extend([-20,0,0,-20,0,0])
            bounds_up.extend([0,np.inf,np.inf,0,2.5,np.inf])

            params, _ = curve_fit(py_fit.gaussian_2, x_data, y_data, p0=initial_guess,bounds=[bounds_low,bounds_up],maxfev=2000000)
            printflag = "3"
        elif len(p)==2:
            initial_guess = []
            bounds_low = []
            bounds_up = []
            for enum_p,p_val in enumerate(p):
                initial_guess.extend([bins[p_val],1,hist[p_val]])
                bounds_low.extend([2*bins[p_val],0,0])
                bounds_up.extend([0,1.75,np.inf])     
            params, _ = curve_fit(py_fit.gaussian_2, x_data, y_data, p0=initial_guess,bounds=[bounds_low,bounds_up],maxfev=2000000)
            printflag = "3"      
        else:
            initial_guess = []
            bounds_low = []
            bounds_up = []
            for enum_p,p_val in enumerate(p):
                initial_guess.extend([bins[p_val],2,hist[p_val]])
                bounds_low.extend([2*bins[p_val],0,0])
                bounds_up.extend([0,np.inf,np.inf])
            params, _ = curve_fit(py_fit.gaussian_X, x_data, y_data, p0=initial_guess,bounds=[bounds_low,bounds_up],maxfev=2000000)
            printflag = "X"
    
        retvalues = []
        for i_p in range(0,len(params),3):
            retvalues.append(params[i_p])
        retvalues = sorted(retvalues)
        retvalues = retvalues[::-1]

        pedestal.append(retvalues[0])
        spe.append(retvalues[1])
        dpe.append(0)
        volt.append(run) 

        err_index_ped = np.where(params==retvalues[0])[0][0]+1
        err_pedestal.append(params[err_index_ped])
        err_index_spe = np.where(params==retvalues[1])[0][0]+1
        err_spe.append(params[err_index_spe])
        err_dpe.append(0)

        if print_single_volt==True:
            axes[rows, cols].plot(bins[np.array(p)],hist[np.array(p)],'xr')
            if printflag == "X":
                axes[rows, cols].plot(x_data,py_fit.gaussian_X(x_data,*params),'-r')
            elif printflag == "3":
                axes[rows, cols].plot(x_data,py_fit.gaussian_2(x_data,*params),'-r')
            for i in range(0,len(params),3):
                axes[rows, cols].plot(x_data,params[2+i]*py_fit.gaussian(x_data,*params[i:2+i]),'--',label=params[i:2+i])
            axes[rows, cols].legend()
        run_index += 1
    

    return volt,pedestal,spe,dpe,err_pedestal,err_spe,err_dpe


def GetVoltageValues(threshold,data,id_list):
    [volt,spe_list,spe_err_list] = data

    fitlist = []
    valuelist = []
    for enum_id,id in enumerate(id_list):
        tmp_x = volt[5:]
        tmp_y = spe_list[enum_id][5:]

        params, _ = curve_fit(py_fit.linear, tmp_x, tmp_y, p0=[-4,0],method='lm',maxfev=2000000)
        fitlist.append(params)

        value = (threshold-params[1])/params[0]
        valuelist.append(value)

    py_plot.PlotCalibrationCurve(np.array(volt),[spe_list,spe_err_list,fitlist,valuelist],id_list)