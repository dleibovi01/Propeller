import numpy as np
import matplotlib.pyplot as plt
import os
import csv



csvfile = open('drag_coeffs.csv', 'w', newline='')
fieldnames = ['Re (based on length)', 'case', 'cd_avg', 'cd_std']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()



for i in range(39):
    data_folder=os.getcwd()+"/detailed_car_"+str(i+1)+"/simpleFoam_steady"
    data = np.loadtxt(data_folder+'/postProcessing/forceCoeffs1/0/coefficient.dat', skiprows=13)

    #print(data[:,0].shape)
    plt.plot(data[:,0], data[:,1], label="cd")
    # plt.plot(data[:,0], data[:,4], label="cl")
    plt.legend()
    plt.ylim((0.2, 0.8))
    plt.xlim((1500, 3000))
    plt.grid()
    plt.savefig('detailed_car_'+str(i+1)+'/simpleFoam_steady/cd_plot')

    print("Average drag coefficient over last 1000 steps", np.mean(data[-1000:,1]))
    line_data = [{'Re (based on length)': 1, 'case': i+1, 'cd_avg': np.mean(data[-1000:,1]), 'cd_std': np.std(data[-1000:,1])}]
    # with open('drag_coeffs.csv', 'w', newline='') as csvfile:
    writer.writerows(line_data)


csvfile.close()

