# Action Phase Jump analysis (APJ) for particle accelerators

This repository holds a ported version of the original python 2.7 __APJ__ project for python 3.

Project originally designed by Fernando Cardona, associated to Universidad Nacional de Colombia, and currently supported by universitys' Particle Accelerator Research Group.


# What's in the repository?

Inside the `apj/` directory the main python scripts can be found:

- __`ActionPhaseJump.py`:__ version 7. Script that implements the APJ method.
- __`get_nom_files.py`:__ version 4. Script that prepares the nominal files to be used.
- __`get_IR_corrs.py`:__ version 2. Script that analyses the error and corrections to apply.
- __`tbt_simulator.py`:__ version 5. Simulates the beams.

Followed by utilery scripts:

- `outliers.py`
- `plot_actionPhase.py`: version 1.
- `plot_expVssim.py`" version 2.
- `utils_ActPhase10.py`



# How to use it?

1. Download the project (or pull it from GitHub) into your machine.
2. Download the `madx` software from CERN official webpage into the projects' folder outside of the `apj/` folder and make it executable.
3. Inside `apj/lhc/` you can find 3 different examples
  - `./commands_IR1exp_set`
  - `./commands_IR1simu_est`
  - `./commands_IR5simu_est`
4. Open the example you want to run and modify the variable `madx_program` to where you saved the `madx` executable. Take into account `home_dir` is automatically set to the current project folder (i.e. `apj/` parent folder).
5. Run the example by executing the file.
