import gmsh
import sys
import numpy as np

def hub_generation(r, sizeNACA):

    gmsh.initialize()

    gmsh.model.add("hub")

    fileName = "hub.stl"
    extension = 0.99
    R = r * extension
    zMax = 0.005
    H = 0.012

    gmsh.model.occ.addPoint(0, -zMax, 0, sizeNACA, 1)
    gmsh.model.occ.addPoint(R, -zMax, 0, sizeNACA, 2)
    gmsh.model.occ.addPoint(0, -zMax, R, sizeNACA, 3)
    gmsh.model.occ.addPoint(-R, -zMax, 0, sizeNACA, 4)
    gmsh.model.occ.addPoint(0, -zMax, -R, sizeNACA, 5)

    gmsh.model.occ.addCircleArc(2, 1, 3, 1)
    gmsh.model.occ.addCircleArc(3, 1, 4, 2)
    gmsh.model.occ.addCircleArc(4, 1, 5, 3)
    gmsh.model.occ.addCircleArc(5, 1, 2, 4)

    gmsh.model.occ.addCurveLoop([1, 2, 3, 4], 1) 
    gmsh.model.occ.addSurfaceFilling(1, 1) 


    extrusion_vector = [0, H, 0]
    hub = gmsh.model.occ.extrude([(2, 1)], dx=extrusion_vector[0], dy=extrusion_vector[1], dz=extrusion_vector[2])


    gmsh.model.occ.synchronize()

    # gmsh.model.occ.synchronize()

    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)

    # ... and save it to disk
    gmsh.write(fileName)

    # # Launch the GUI to see the model:
    # if '-nopopup' not in sys.argv:
    #    gmsh.fltk.run()

    gmsh.finalize()
