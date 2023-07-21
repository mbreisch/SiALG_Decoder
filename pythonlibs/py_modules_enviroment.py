import matplotlib.pyplot as plt
import numpy as np

def GetVoltageStability(file_path):
    with open(file_path, 'r') as file:
        path_list = file.read().splitlines()

    dict = {}

    for entry in path_list:
        date,channel,value = entry.split(";")

        if channel not in dict:
            dict.update({channel:[]})
        else:
            dict[channel].append(float(value))

    set_to_channel = [26,26.5,27,27.5,28,28.5,29,29.5]    
    bin_width = 0.005
    for channel,list_of_values in dict.items():
        fig = plt.figure(channel,figsize=[7,7], facecolor='white')
        list_of_values = np.array(list_of_values)
        plt.hist(list_of_values,bins=np.arange(min(list_of_values),max(list_of_values)+bin_width,bin_width),label="Dist={}".format(set_to_channel[int(channel)]-np.mean(list_of_values)))
        plt.axvline(x=set_to_channel[int(channel)])
        plt.legend()
