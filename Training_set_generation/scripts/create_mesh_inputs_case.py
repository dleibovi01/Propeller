import os
import numpy as np
import math



def read_file(fp):
    fd = open( fp, "rt")
    data = fd.read()
    fd.close()
    return data

def write_file(fp,data):
    fd = open( fp, "wt")
    fd.write(data)
    fd.close()



if __name__ == "__main__":


    inlet_velocity = "-4.96"

    print("Creating mesh input data")
    data_folder=os.getcwd()
    data=np.genfromtxt(data_folder+'/case_info.txt',delimiter=": ",usecols=(1)).tolist()
    print("\n")
    stl_file = data_folder+"/simpleFoam_steady/constant/triSurface/propeller.stl"
    print(stl_file)
    if os.path.exists(stl_file):

        print("Creating the necessary mesh and solver inputs")

        # REFINEMENT_INPUTS=create_mesh_input(data)
        model_dir = data_folder + "/simpleFoam_steady/0_org/include/"

        print("Replacing in initialConditions")
        txt = read_file(model_dir+"initialConditions")
        txt = txt.replace(inlet_velocity, str(data[0]))

        write_file(model_dir+"initialConditions",txt)
        print("Finished wiriting in initialConditions")

    else:
        print("No inputs needed as geometry is not available\n\n")
        # os.system("rm -r  f "+data_folder+"/simpleFoam_steady")

