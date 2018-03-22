# Configure a Simulation
HYDRAD requires setting many different configuration files, making it difficult and tedious to configure a simulation by hand. Fortunately, the `hydrad_tools` package makes it easy to quickly configure a new simulation using Python.

## Setting a Default Configuration

## Running the Simulation

## Configuration Parameters
The tables below give an exhaustive list of all of the different HYDRAD configuration options. If the units are listed, the input must have units that can be converted to the listed unit with the [Astropy units module](http://docs.astropy.org/en/stable/units/), e.g. `loop_length` can be input in Mm.

### General
| Name | Description | Type | Units |
|:----:|:-----------|:----:|:-----:|
| total_time | Total duration of the simulation | `int` | s |
|output_interval | How often results are printed to file | `int` | s |
| loop_length | Footpoint-to-footpoint distance of the coronal loop | `float` | cm |
| loop_inclination | Angle between loop and surface normal | `float` | degree |
| footpoint_height | Height of loop footpoint above the solar surface | `float` | cm |

### Initial Conditions
| Name | Description | Type | Units |
|:----:|:-----------|:----:|:-----:|
| footpoint_temperature | Temperature at the loop footpoint | `float` | K
| footpoint_density | Density at the loop footpoint | `float` | cm<sup>-3</sup>|

### Heating

### Radiation

### Solver 

### Grid

