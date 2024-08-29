
# pH-dependentCG

Welcome to the **pH-dependentCG** repository! This project provides tools to run molecular dynamics (MD) simulations with pH-dependent coarse-grained models.

## How to Run the Simulation

To compile and run the simulation:

1. Use the `MDmake` file by executing one of the following commands in your terminal:
   ```bash
   make -f MDmake
   ```
   or
   ```bash
   ./MDmake
   ```
2. This will generate an executable file named `MD.exe`.

## Important Notes

1. **Settings and Parameters:**
   - The `ss.f` file reads the `settings.dat` file, which contains all the parameters for the simulation, including salt concentration, friction constant, and various flags.
   - You can find an example `settings.dat` file [here](./settings_54237.dat).

2. **Input Files:**
   - The `init.f` file reads an input file named `Hst5_up.dat` (as specified in `settings.dat`). This file defines the beads' mass, electrostatic charge, and short-range potential coefficients.
   - If you want to use initial conditions from a previous run, set the option in the `settings.dat` file as follows:
     ```text
     3 : whether to use initial conditions
     ```
   - You can also use a final step file, such as `Final_Hst5_up_0.45_54237.dat`, which includes velocity data.

3. **Reduced Units:**
   - This simulation uses reduced units for parameters like temperature, mass, friction coefficient, and ionic strength. Note that the ionic strength is scaled to be approximately three times larger than the input value. This approach is consistent with the methodologies established in the following references:
     - [Reference 1](https://doi.org/10.1093/nar/gkz1202)
     - [Reference 2](https://doi.org/10.1021/acs.jpcb.1c00757)
     - [Reference 3](https://doi.org/10.1016/j.jmb.2009.08.010)

4. **pH-Dependent Histidine Protonation:**
   - An updated version of the simulation, which allows for switching the charge state of Histidine residues (and their short-range interaction strength) based on pH and their chemical environment, is currently in the testing phase. Stay tuned for updates!


    