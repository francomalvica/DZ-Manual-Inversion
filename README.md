# DZ-Manual-Inversion
Julia-based framework for 1D geoelectrical forward modeling and inversion, featuring apparent resistivity computation via O’Neill filter and Sunde resistivity transform, CRV/CRA generation, Dar Zarrouk parameterization, iterative layer merging, manual model editing, and full tools for synthetic and observational VES analysis.


# Modeling and Fitting Apparent Resistivity Curves

This project provides a set of Julia tools for:

  Generating synthetic   apparent resistivity curves (ARC)   from a user-defined model.
  Applying the   Dar–Zarrouk   method to reduce the number of layers.
  Performing   manual adjustment   of both resistivities and thicknesses.
  Comparing the observed ARC with the modeled one to evaluate the inversion.

---

## 1. Synthetic Data Generation (run_filter.jl)

The file   run_filter.jl   allows you to generate apparent resistivity curves from any invented true-resistivity model.
This is useful for producing test datasets to validate the inversion routines in the Dar–Zarrouk and manual adjustment module.

The script:

  Accepts a layered model of true resistivities.
  Applies the   O’Neil filter   to compute the corresponding apparent resistivity curve.
  Produces an output file containing the ARC and its AB/2 values.
  Plots the resulting synthetic curve.

This tool is ideal for experimenting with different models before applying the inversion to real data.

---

## 2. Dar–Zarrouk and Manual Adjustment (main.jl)

Inside the  Dar–Zarrouk and manual  folder you will find   main.jl  , which performs the inversion workflow.

This program allows you to:

### • Apply the Dar–Zarrouk method

Automatically reduces the number of layers while maintaining electrical equivalence.

### • Perform interactive manual fitting

You can adjust:

  Layer resistivities
  Layer thicknesses

until the modeled ARC matches the observed ARC.
The program continuously displays the   RMS error   to guide the inversion process.

---

## 3. Recommended Workflow

1. Generate a synthetic ARC by running:

   ```
   julia run_filter.jl
   ```
2. This creates   RA.dat   (or a similar file), which serves as input for the inversion.
3. Navigate to the  Dar–Zarrouk and manual  directory and run:

   ```
   julia main.jl
   ```
4. Use the Dar–Zarrouk automatic reduction or manually adjust the model until you obtain an acceptable RMS error.

---

## 4. Main Files

| File                   | Purpose                                                            |
| ---------------------- | ------------------------------------------------------------------ |
|   run_filter.jl        | Generates synthetic apparent resistivity curves.                   |
|   algoritmo_zohdy.jl   | Computes an initial layered model based on Zohdy’s approach.       |
|   filtro_oneil.jl      | Implements the O’Neil filter for ARC calculation.                  |
|   main.jl              | Performs the Dar–Zarrouk method and interactive manual adjustment. |
|   funcion_dar.jl       | Functions used for the Dar–Zarrouk reduction.                      |
