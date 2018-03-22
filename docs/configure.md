# Configure a Simulation
HYDRAD requires setting many different configuration files, making it difficult and tedious to configure a simulation by hand. Fortunately, the `hydrad_tools` package makes it easy to quickly configure a new simulation using Python.

## Setting a Default Configuration

## Running the Simulation

## Configuration Parameters
The table below gives an exhaustive list of all of the different HYDRAD configuration options

<table>
<!-- -->
<tr><th>General</th></tr>
<tr><th>Name</th><th>Description</th><th>Type</th></tr>
<tr>
    <td> total_time </td>
    <td> Total duration of the simulation </td>
    <td> int </td>
</tr>
<tr>
    <td> output_interval </td>
    <td> How often results are printed to file </td>
    <td> int </td>
</tr>
<tr>
    <td> loop_length </td>
    <td> Footpoint-to-footpoint distance of the coronal loop </td>
    <td> float </td>
</tr>
<tr>
    <td> loop_inclination </td>
    <td> Angle between loop and surface normal </td>
    <td> float </td>
</tr>
<tr>
    <td> footpoint_height </td>
    <td> Height of loop footpoint above the solar surface </td>
    <td> float </td>
</tr>
<!-- -->
<tr><th>Initial Conditions</th></tr>
<tr>
    <td> footpoint_temperature </td>
    <td> Temperature at the loop footpoint </td>
    <td> float </td>
</tr>
<tr>
    <td> footpoint_density </td>
    <td> Density at the loop footpoint </td>
    <td> float </td>
</tr>
<!-- -->
<tr><th>Heating</th></tr>
<!-- -->
<tr><th>Radiation</th></tr>
<!-- -->
<tr><th>Solver</th></tr>
<!-- -->
<tr><th>Grid</th></tr>
</table>
