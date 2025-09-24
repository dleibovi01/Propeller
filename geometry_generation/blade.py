import gmsh
import sys
import numpy as np
from scipy.interpolate import PchipInterpolator, pchip_interpolate


def generate_blade(Diam, r, Pitch_D, skew, NACA, c, nline, npoint, nradii, npoint_circle=10000, Ni=100, Nsplit=20, rounding=0.65, sizeNACA=5e-5):
        
        Pitch = Diam * Pitch_D # Input pitch data
        alpha = np.arctan(Pitch / (2 * np.pi * r)) # pitch angle
        theta_s = np.pi / 180 * skew # 10.76 degree skew
        S = r * theta_s / np.cos(alpha) # Skew

        x_center = 2
        y_center = 0
        R_circle = 4

        beta=np.linspace(0,rounding * np.pi,Ni)
        x = c*(0.5*(1-np.cos(beta)))
        z_c = np.zeros_like(x)
        z_t = np.zeros_like(x)
        theta = np.zeros_like(x)

    
        m = NACA[0]/100  # Values of m, p and t
        p = NACA[1]/10
        t = NACA[2]*10 + NACA[3]/100

        z_t = (t*c/0.2) * (0.2969*(x/c)**0.5 - 0.1260*(x/c) - 0.3516*(x/c)**2 + 0.2843*(x/c)**3 - 0.1036*(x/c)**4)

        
        if (p > 0): 
            z_c = z_c + (m*x/p**2) * (2*p - x/c) * (x < p*c) # Calculate camber
            z_c = z_c + (m*(c-x)/(1-p)**2) * (1 + x/c - 2*p) * (x >= p*c)

            theta = theta + np.arctan( (m/p**2) * (2*p - 2*x/c) ) * (x < p*c) #  % Calculate theta-value
            theta = theta + np.arctan( (m/(1-p)**2) * (-2*x/c + 2*p) ) * (x >= p*c)
        
        
        Xu = x - z_t*np.sin(theta) # % Calculate coordinates of upper surface
        Zu = z_c + z_t*np.cos(theta)

        Xl = x + z_t*np.sin(theta) # Calculate coordinates of lower surface
        Zl = z_c - z_t*np.cos(theta)

        upper = np.concatenate(([Xu], [Zu]), axis=0)
        lower = np.concatenate(([Xl], [Zl]), axis=0)

        dXu = Xu[-1] - Xu[-2]
        dZu = Zu[-1] - Zu[-2]

        dXl = Xl[-1] - Xl[-2]
        dZl = Zl[-1] - Zl[-2]

        dU = np.concatenate(([dXu], [dZu]), axis=0)
        dL = np.concatenate(([dXl], [dZl]), axis=0)

        M = np.concatenate(([np.transpose(dU)], [np.transpose(dL)]), axis=0)
        b = np.concatenate(([dU[0] * Xu[-1] + dU[1] * Zu[-1]], [dL[0] * Xl[-1] + dL[1] * Zl[-1]]), axis=0)

        intersect = np.linalg.solve(M, b)   

        # print(intersect)


        mid = 0.5 * np.array([Xu[-1] + Xl[-1], Zu[-1] + Zl[-1]])
        R = np.linalg.norm((intersect - np.array([Xu[-1], Zu[-1]])))
        edge = intersect + R * (mid - intersect) / np.linalg.norm(mid - intersect)

        # print(edge)

        X = np.concatenate((upper[0, :], lower[0, Ni-1:0:-1], np.array([intersect[0]]), np.array([edge[0]])), axis=0)
        Z = np.concatenate((upper[1, :], lower[1, Ni-1:0:-1], np.array([intersect[1]]), np.array([edge[1]])), axis=0)

        # print(X.shape)

        theta_u = np.arctan((Zu[-1] - intersect[1]) / (Xu[-1] - intersect[0]))
        theta_mid = np.arctan((edge[1] - intersect[1]) / (edge[0] - intersect[0]))
        theta_l = np.arctan((Zl[-1] - intersect[1]) / (Xl[-1] - intersect[0]))

        N_circle = 100
        theta_up = np.flip(np.linspace(theta_mid, theta_u, N_circle + 2))
        theta_down = np.linspace(theta_l, theta_mid, N_circle + 2)
        circle_up = intersect + np.stack((R * np.cos(theta_up[1 : -1]), R * np.sin(theta_up[1 : -1])), axis=1)
        circle_down = intersect + np.stack((R * np.cos(theta_down[1 : -1]), R * np.sin(theta_down[1 : -1])), axis=1)


        Rot = np.array([[np.cos(alpha), np.sin(alpha)], [- np.sin(alpha), np.cos(alpha)]])
        V = np.array([[S], [0]]) + np.matmul(Rot, np.stack((X - 0.5 * c, Z), axis=0))

        X = V[0, :]
        Z = V[1, :]

        circle_up = np.array([[S], [0]]) + np.matmul(Rot, np.stack((circle_up[:, 0] - 0.5 * c, circle_up[:, 1]), axis=0)) 
        circle_down = np.array([[S], [0]]) + np.matmul(Rot, np.stack((circle_down[:, 0] - 0.5 * c, circle_down[:, 1]), axis=0)) 

        N = len(X)


        X_cyl = r / np.sqrt(X**2 + Z**2 + r**2) * X
        # X_cyl = r * X / np.sqrt(X**2 + r**2)
        Y_cyl = r / np.sqrt(X**2 + Z**2 + r**2) * Z
        # Y_cyl = Z
        Z_cyl = r / np.sqrt(X**2 + Z**2 + r**2) * r
        # Z_cyl = r**2 / np.sqrt(X**2 + r**2)

        for i in range(N):
            gmsh.model.occ.addPoint(X_cyl[i], Y_cyl[i], Z_cyl[i], sizeNACA, npoint + i + 1)

        X_cup = r / np.sqrt(circle_up[0, :]**2 + circle_up[1, :]**2 + r**2) * circle_up[0, :]
        # X_cup = r * circle_up[0, :] / np.sqrt(circle_up[0, :]**2 + r**2)
        Y_cup = r / np.sqrt(circle_up[0, :]**2 + circle_up[1, :]**2 + r**2) * circle_up[1, :]
        # Y_cup = circle_up[1, :]
        Z_cup = r ** 2 / np.sqrt((circle_up[0, :]**2 + circle_up[1, :]**2 + r**2)) 
        # Z_cup = r ** 2 / np.sqrt((circle_up[0, :]**2 + r**2)) 
        
        for i in range(N_circle):
            gmsh.model.occ.addPoint(X_cup[i], Y_cup[i], Z_cup[i], sizeNACA, npoint_circle + i + 1)

        X_cdown = r / np.sqrt(circle_down[0, :]**2 + circle_down[1, :]**2 + r**2) * circle_down[0, :]
        # X_cdown = r * circle_down[0, :] / np.sqrt(circle_down[0, :]**2 + r**2)
        Y_cdown = r / np.sqrt(circle_down[0, :]**2 + circle_down[1, :]**2 + r**2) * circle_down[1, :]
        # Y_cdown = circle_down[1, :]
        Z_cdown = r ** 2 / np.sqrt((circle_down[0, :]**2 + circle_down[1, :]**2 + r**2))
        # Z_cdown = r ** 2 / np.sqrt((circle_down[0, :]**2 + r**2))
        
        for i in range(N_circle):
            gmsh.model.occ.addPoint(X_cdown[i], Y_cdown[i], Z_cdown[i], sizeNACA, npoint_circle + N_circle + i + 1)        

        gmsh.model.occ.synchronize()

        gmsh.model.occ.addSpline([*range(npoint + 1, int(npoint + Nsplit + 1), 1)], nline + 1)
        gmsh.model.occ.addSpline([*range(npoint + Nsplit, int(npoint + Ni + 1), 1)], nline + 2)
        gmsh.model.occ.addSpline([*range(npoint + Ni + 1, int(npoint + 2 * Ni - Nsplit + 1), 1)], nline + 5)
        elems = [*range(int(npoint + 2 * Ni - Nsplit), int(npoint + 2 * Ni), 1)]
        elems.append(npoint + 1)
        gmsh.model.occ.addSpline(elems, nline + 6)

        elems = [*range(int(npoint_circle + 1), int(npoint_circle + N_circle + 1), 1)]
        elems.insert(0, npoint + Ni)
        elems.append(npoint + N)
        gmsh.model.occ.addSpline(elems, nline + 3)

        elems = [*range(int(npoint_circle + N_circle + 1), int(npoint_circle + 2 * N_circle + 1), 1)]
        elems.insert(0, npoint + Ni + 1)
        elems.append(npoint + N)
        gmsh.model.occ.addSpline(elems, nline + 4)


        npoint += N
        nline += 6
        npoint_circle += 2 * N_circle

        return npoint, nline, npoint_circle    


# def body_generation(diameter=0.2286, sizeNACA=1e-4):
def blade_generation(radiusq, betas, Pitches_Dq, Cq, t0_D, NACAsq, Skewsq, diameter=0.2286, sizeNACA=1e-4):

    gmsh.initialize()

    gmsh.model.add("propeller")

    # Inputs: For now this is just hardcoded, but will have to be passed to the function
    # Input for GWS EP-9050

    # Diam = diameter
    # Diam = 0.2286
    # Rmax = Diam / 2.0
    # Nz = 60

    # radiusq = np.array([0.15, 0.2, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1]) * Rmax
    # r0 = radiusq[0]
    # betas = np.array([25.65, 26.35, 26.75, 26.49, 25.57, 23.95, 21.72, 19.38, 17.53, 15.80, 14.40, 13.19, 12.09, 11.38, 10.13, 8.70, 6.5, 4.26]) * np.pi / 180
    # Pitches_Dq = np.tan(betas) * 2 * np.pi * radiusq / Diam
    # Cq = np.array([0.169, 0.181, 0.183, 0.183 , 0.182, 0.180, 0.176, 0.173, 0.170, 0.167, 0.162, 0.156, 0.150, 0.142, 0.132, 0.117, 0.095, 0.073]) * Diam * 0.5

    # t0_D = 60 * Cq / Diam


    # NACAsq = np.array([[0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2], 
    #         [0, 0, 1, 2]])
    # NACAsq[:, 2] = NACA3s
    # NACAsq[:, 3] = NACA4s
    # Skewsq = np.zeros_like(radiusq)


    # NACA3s = np.floor(100 * t0_D * Diam / Cq / 1000 / 10)
    # NACA4s = 100 * t0_D * Diam / Cq / 1000 - 10 * NACA3s
    # NACAsq[:, 2] = NACA3s
    # NACAsq[:, 3] = NACA4s

    Diam = diameter
    Rmax = Diam / 2.0
    Nz = 60    
    correction_ratio = Pitches_Dq / Cq * Diam
    rounding = 0.65



    # Maybe some pre-processing to do from the inputs

    rmin = radiusq[0]
    rmax = radiusq[-1]
    radius1 = np.linspace(rmin, rmax, Nz)
    radius = radius1
    Pitches_D = pchip_interpolate(radiusq, Pitches_Dq, radius)
    C = pchip_interpolate(radiusq, Cq, radius)   
    NACAs = np.zeros((Nz, 4))
    for i in range(4):
        NACAs[:, i] = pchip_interpolate(radiusq, NACAsq[:, i], radius)
    Skews = pchip_interpolate(radiusq,Skewsq,radius)
    correction_ratios = C / (NACAs[:, 2]*10 + NACAs[:, 3]) / Diam
    correction_ratios = correction_ratios / correction_ratios[0]




    # Building the blade by superposing the airfoils

    nline = 0
    npoint = 0
    nradii = len(radius)
    npoint_circle = 100000
    Ni = 100; 
    Nsplit = int(np.round(0.2* Ni)) 
    rounding = 0.65   

    for i in range(nradii):
        Pitch_D = Pitches_D[i]
        r = radius[i]
        skew = Skews[i]
        NACA = NACAs[i, :]
        c = C[i]  
        npoint, nline, npoint_circle = generate_blade(Diam, r, Pitch_D, skew, NACA, c, nline, npoint, nradii, npoint_circle, Ni=100, Nsplit=20, rounding=0.65, sizeNACA=sizeNACA)
            
    
    gmsh.model.occ.synchronize()

    # ## 2 new lines
    # gmsh.model.occ.addLine(Nsplit, 2 * Ni - Nsplit, nline + 1)
    # gmsh.model.occ.addLine(Ni, Ni + 1, nline + 2)

    # gmsh.model.occ.synchronize()

    # # Launch the GUI to see the model:
    # if '-nopopup' not in sys.argv:
    #    gmsh.fltk.run()
    
    nlines_airfoil = nline

    for i in range(nradii - 1):
        gmsh.model.occ.addLine(i * (2 * Ni + 1) + 1, (i + 1) * (2 * Ni + 1) + 1, nline + 1)
        gmsh.model.occ.addLine(i * (2 * Ni + 1) + Nsplit, (i + 1) * (2 * Ni + 1) + Nsplit, nline + 2)
        gmsh.model.occ.addLine(i * (2 * Ni + 1) + Ni, (i + 1) * (2 * Ni + 1) + Ni, nline + 3)
        gmsh.model.occ.addLine(i * (2 * Ni + 1) + 2 * Ni + 1, (i + 1) * (2 * Ni + 1) + 2 * Ni + 1, nline + 4)
        gmsh.model.occ.addLine(i * (2 * Ni + 1) + Ni + 1, (i + 1) * (2 * Ni + 1) + Ni + 1, nline + 5)
        gmsh.model.occ.addLine(i * (2 * Ni + 1) + 2 * Ni - Nsplit, (i + 1) * (2 * Ni + 1) + 2 * Ni - Nsplit, nline + 6)
        nline += 6

    gmsh.model.occ.synchronize()



    nline_per_level = 6
    nloop_per_level = 24
    nline_per_blade = 6
    buffer = 1000

    for i in range(nradii - 1):
        print(i)
        # gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 1, -(nlines_airfoil + nline_per_level * i + 2), -(nline_per_blade * i + 1), nlines_airfoil + nline_per_level * i + 1], nloop_per_level * i + 1)        
        # gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 2, -(nlines_airfoil + nline_per_level * i + 3), -(nline_per_blade * i + 2), nlines_airfoil + nline_per_level * i + 2], nloop_per_level * i + 3)
        # gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 3, -(nlines_airfoil + nline_per_level * i + 4), -(nline_per_blade * i + 3), nlines_airfoil + nline_per_level * i + 3], nloop_per_level * i + 4)
        # gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 4, -(nlines_airfoil + nline_per_level * i + 5), -(nline_per_blade * i + 4), nlines_airfoil + nline_per_level * i + 4], nloop_per_level * i + 5)
        # gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 5, -(nlines_airfoil + nline_per_level * i + 6), -(nline_per_blade * i + 5), nlines_airfoil + nline_per_level * i + 5], nloop_per_level * i + 6)
        # gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 6, -(nlines_airfoil + nline_per_level * i + 1), -(nline_per_blade * i + 6), nlines_airfoil + nline_per_level * i + 6], nloop_per_level * i + 8)


        # gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 1, nloop_per_level * i + 1) 
        # gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 3, nloop_per_level * i + 3)  
        # gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 4, nloop_per_level * i + 4)   
        # gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 5, nloop_per_level * i + 5)      
        # gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 6, nloop_per_level * i + 6) 
        # gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 8, nloop_per_level * i + 8)          


        gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 1, -(nlines_airfoil + nline_per_level * i + 2), -(nline_per_blade * i + 1), nlines_airfoil + nline_per_level * i + 1], nloop_per_level * i + 1)        
        gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 2, -(nlines_airfoil + nline_per_level * i + 3), -(nline_per_blade * i + 2), nlines_airfoil + nline_per_level * i + 2], nloop_per_level * i + 2)
        gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 3, -(nlines_airfoil + nline_per_level * i + 4), -(nline_per_blade * i + 3), nlines_airfoil + nline_per_level * i + 3], nloop_per_level * i + 3)
        gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 4, -(nlines_airfoil + nline_per_level * i + 5), -(nline_per_blade * i + 4), nlines_airfoil + nline_per_level * i + 4], nloop_per_level * i + 4)
        gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 5, -(nlines_airfoil + nline_per_level * i + 6), -(nline_per_blade * i + 5), nlines_airfoil + nline_per_level * i + 5], nloop_per_level * i + 5)
        gmsh.model.occ.addCurveLoop([nline_per_blade * (i + 1) + 6, -(nlines_airfoil + nline_per_level * i + 1), -(nline_per_blade * i + 6), nlines_airfoil + nline_per_level * i + 6], nloop_per_level * i + 6)

        gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 1, nloop_per_level * i + 1) 
        gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 2, nloop_per_level * i + 2)  
        gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 3, nloop_per_level * i + 3)   
        gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 4, nloop_per_level * i + 4)      
        gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 5, nloop_per_level * i + 5) 
        gmsh.model.occ.addSurfaceFilling(nloop_per_level * i + 6, nloop_per_level * i + 6)   
         
    gmsh.model.occ.synchronize()

    # # Launch the GUI to see the model:
    # if '-nopopup' not in sys.argv:
    #    gmsh.fltk.run()

    # Add the surface at the bottom and at the top
    gmsh.model.occ.addCurveLoop([1, 2, 3, 4, 5, 6], (nradii - 1) * nloop_per_level + 1) 
    gmsh.model.occ.addSurfaceFilling((nradii - 1) * nloop_per_level + 1, (nradii - 1) * nloop_per_level + 1) 
    gmsh.model.occ.addCurveLoop([nlines_airfoil - 5, nlines_airfoil - 4, nlines_airfoil - 3, nlines_airfoil - 2, nlines_airfoil - 1, nlines_airfoil], (nradii - 1) * nloop_per_level + 3) 
    gmsh.model.occ.addSurfaceFilling((nradii - 1) * nloop_per_level + 3, (nradii - 1) * nloop_per_level + 3)    


    gmsh.model.occ.synchronize()

    # gmsh.model.occ.synchronize()

    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)

    # ... and save it to disk
    gmsh.write("blade1.stl")

    # # Launch the GUI to see the model:
    # if '-nopopup' not in sys.argv:
    #    gmsh.fltk.run()


    gmsh.finalize()


if __name__ == "__main__":
    body_generation()