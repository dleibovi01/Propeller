from blade import blade_generation
from hub import hub_generation
import gmsh
import numpy as np
import os
import csv
import pandas as pd
import argparse

# Masctsr script for generating the propeller

# Input some parameters for NACA, chord length, pitch, skew, Diameter

def propeller_generation(case, Diam = 0.2286):

    # Read data from CSV file
    Nz = 60
    sizeNACA = 1e-4


    # df0 = pd.read_csv('case_info.csv', usecols=[1])
    df = pd.read_csv('case_info.csv')
    radiusq = np.array(df['radiusq']) 
    # Diam = radiusq[-1]

    # print(Diam)
    Rmax = Diam / 2.0
    # radiusq = radiusq[:-1] * Rmax
    radiusq = radiusq[:-1] * Rmax

    r0 = radiusq[0]
    # print(r0)
    # input()
    betas = np.array(df['beta']) * np.pi / 180
    betas = betas[:-1]
    Pitches_Dq = np.tan(betas) * 2 * np.pi * radiusq / Diam
    Cq = np.array(df['Cq']) * Diam * 0.5
    Cq = Cq[:-1]
    t0_D = 60 * Cq / Diam

    NACA1s = np.array(df['NACA1s'])
    NACA2s = np.array(df['NACA2s'])
    NACA1s = NACA1s[:-1]
    NACA2s = NACA2s[:-1]
    NACA3s = np.floor(100 * t0_D * Diam / Cq / 1000 / 10)
    NACA4s = 100 * t0_D * Diam / Cq / 1000 - 10 * NACA3s  

    NACAsq = np.zeros((len(radiusq), 4))  

    NACAsq[:, 0] = NACA1s
    NACAsq[:, 1] = NACA2s
    NACAsq[:, 2] = NACA3s
    NACAsq[:, 3] = NACA4s

    Skewsq = np.array(df['Skewsq'])
    Skewsq = Skewsq[:-1]
 
    # Generate blade

    blade_generation(radiusq, betas, Pitches_Dq, Cq, t0_D, NACAsq, Skewsq, diameter=Diam, sizeNACA=sizeNACA)

    # Call OpenFoam script to duplicate / rotate it

    os.system("surfaceTransformPoints -rotate-y 180 blade1.stl blade2.stl")

    # Create cylinder

    hub_generation(radiusq[0], sizeNACA)

    # Merge the three pieces into one stl

    gmsh.initialize()
    gmsh.model.add("propeller")

    gmsh.merge('blade1.stl')
    gmsh.merge('blade2.stl')
    gmsh.merge('hub.stl')

    gmsh.model.occ.synchronize()

    # ... and save it to disk
    gmsh.write("propeller.stl")

    # # Launch the GUI to see the model:
    # if '-nopopup' not in sys.argv:
    #     gmsh.fltk.run()


    gmsh.finalize()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that generates the blade for the case given in CMD"
    )

    parser.add_argument("--case", default=0, type=int)
    args = parser.parse_args()
    case = args.case

    propeller_generation(case)