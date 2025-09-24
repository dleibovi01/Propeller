# Propeller dataset generation

A framework for generating geometries and running RANS simulations for rotating propellers

## Geometry generation

First, generate the cases by running:

```
bash generate_cases [case_start] [case_end]
```

This will generate test cases in the `geometry_generation/cases` directory. 
This includes a csv file detailing the geometric parameters of the propeller to produce, and the simulation parameters (inlet and rotation speed)

Then, generate the propellers by running 
```
bash timeout_launch_range_stl_generations.sh [case_start] [case_end]
```

This will attempt to generate the geometries, in the `geometry_generation/propeller_geometries` directory. Some generations might fail because of surface self-intersection; the generation script
will atempt to generate them, and timeout after 120s if the generation fails. It is then necessary to run the following script:
```
bash cleanup_geometries_dir.sh [case_start] [case_end]
```
This will automatically remove the subdirectories of `geometry_generation/propeller_geometries` for which no .stl file was generated
