# pbpk-lite

`pbpk-lite` is a lightweight Python package for basic physiologically based pharmacokinetic (PBPK) modeling.

It provides a simple programmatic interface for defining substance and metabolite properties, patient physiology, administration regimens, and elimination kinetics, and for solving the resulting ODE system.

The implementation uses nanograms (ng) for doses and amounts, milliliters (mL) for volumes, milliliters per minute (mL/min) for blood flows and clearances, nanograms per milliliter (ng/mL) for concentrations, and minutes for time. The simulation uses minutes internally, while the plotting helpers can display the x-axis in minutes, hours, or days via the optional `time_unit` argument.

## Features

- Substance and metabolite partition coefficient calculation from logP and fraction unbound
- Patient blood flows and tissue volumes derived from body weight
- Linear metabolism and excretion pathways
- Oral absorption support with tunable Weibull parameters
- ODE solution via `scipy.integrate.solve_ivp`
- Plotting helpers for whole-model, venous-blood, and selected-compartment concentration profiles
- Support for intravenous, inhalation, and oral administration
- Simple optimization helpers for fitting clearances and oral parameters

The model uses 16 compartments for each substance. Simulation results contain
one block of 16 concentration rows per substance; administered doses enter the
parent substance block.

## Installation

Install from PyPI:

```bash
pip install pbpk-lite
```

## Quick Start

```python
from pbpk_lite import model

m = model()
m.set_substance_and_metabolites(
    (6.97, 5.33, 5.24),
    (0.0022448, 0.01209, 0.101),
    (314.469, 330.468, 344.451),
)
m.set_patient(bw=70)
dose = 25e6  # dose in ng
m.set_administration([('iv', dose, 0, 2)])
m.set_metabolism_and_excretion(
    ((0, 1), (1, 2), (2, -1), (2, -1)),
    ((6, 247), (6, 998), (6, 8.94), (5, 1.48)),
)

t, c = m.simulate(additional_time=60 * 24)

m.graph_whole('concentrations.png')
m.graph_venous('venous.png', limit_of_detection=0.15, time_unit='hours')
m.graph_compartments([(0, 6), (0, 5)], 'selected.png', time_unit='days')

# The simulation uses minutes internally; only the displayed axis changes.
```

## Route of Administration

Pass an administration regimen to `set_administration()` to control where and
when doses enter the model. Each entry has the form
`(route, dose, start_time, administration_time)`.

Doses can only be administered to the parent substance. Metabolites are
produced by the model and cannot be dosed directly.

Supported values are:

- `iv`: intravenous infusion into the venous blood compartment
- `inh`: inhalation into the lung compartment
- `oral`: oral administration through the gut compartment

For oral dosing, the optional `oral_params` argument accepts a tuple of the form
`(t_lag, scaling_factor, shape_parameter)`. The default values are `(20, 450, 2)`.

Example:

```python
m.set_administration([('inh', dose, 0, 10)])
t, c = m.simulate(additional_time=24 * 60)

m.set_administration([('oral', dose, 0, 0)], oral_params=(20, 450, 2))
t, c = m.simulate(additional_time=24 * 60)
```

## Dosing Schedule

The `set_administration()` method expects:

- `regimen`: a sequence of `(route, dose, start_time, administration_time)` entries
- `oral_params`: an optional tuple `(t_lag, scaling_factor, shape_parameter)` for oral dosing only

Times, durations, and the `additional_time` simulation argument are interpreted
in minutes. Doses are expressed in nanograms.

Example with one dose:

```python
# one intravenous infusion at time 0, followed by 24 hours of simulation
regimen = [('iv', dose, 0, 2)]
m.set_administration(regimen)
t, c = m.simulate(additional_time=60 * 24)
```

Example with two doses:

```python
# two identical doses one hour apart, followed by a final observation at 24 hours
regimen = [('iv', dose, 0, 2), ('iv', dose, 60, 2)]
m.set_administration(regimen)
t, c = m.simulate(additional_time=60 * 24)
```

## API Summary

### `pbpk_lite.model`

#### `set_substance_and_metabolites(log_ps, fus, mms)`

Set the physicochemical properties and molecular masses for the substance and
its metabolites. All three sequences should contain one value for each
substance.

- `log_ps`: log octanol-water partition coefficients
- `fus`: fractions unbound in blood
- `mms`: molecular masses; the first entry is the parent substance

#### `set_patient(bw)`

Set patient physiological parameters using body weight in kilograms. Blood
flows and tissue volumes are recalculated for the new body weight.

#### `set_administration(regimen, oral_params=(20, 450, 2))`

Set the administration regimen. Each entry is
`(route, dose, start_time, administration_time)`. Supported routes are `iv`,
`inh`, and `oral`; times and doses are in minutes and nanograms.

For oral dosing, `oral_params` is a tuple of `(t_lag, scaling_factor,
shape_parameter)` used to shape the absorption profile.

#### `set_metabolism_and_excretion(connections, clearances)`

Set the metabolism and excretion pathways. Each connection is a pair of
substance indices, and each clearance is a `(compartment, clearance)` pair.

- `connections`: parent-to-metabolite connections; use `-1` for excretion
- `clearances`: compartment indices and linear clearance values

#### `simulate(additional_time=60)`

Simulate the configured PBPK model and return time points `t` and compartment
concentrations `c`. The simulation continues for `additional_time` minutes
after the final administration time. Administered doses enter the parent
substance block; returned concentrations have 16 rows per substance.

#### `fit_clearances(t, c, report=True)`

Optimize the linear clearance values against experimental concentration data.
This helper is most useful when fitting an IV dataset before refining oral
absorption parameters.

#### `fit_oral(t, c, report=True)`

Optimize the oral absorption parameters against concentration-time data.
The method updates the default oral parameters used by `set_administration()`.

#### `graph_whole(name, time_unit='min')`

Save a multi-panel plot of concentrations across all compartments. The optional `time_unit` argument can be `'min'`, `'hours'`, or `'days'` to change the x-axis label and values. When multiple substances are modeled, one figure is saved per substance with an index appended to the file name.

#### `graph_venous(name, limit_of_detection=None, log=True, time_unit='min', substance_index=0)`

Save a plot of venous blood concentrations, optionally marking a detection limit. The optional `log` argument controls whether the y-axis uses a logarithmic scale, and `time_unit` can be `'min'`, `'hours'`, or `'days'` to change the x-axis label and values.

#### `graph_compartments(substance_and_compartment_indexes, name, time_unit='min')`

Save a plot of selected compartments using pairs of the form `(substance_index,
compartment_index)`. The optional `time_unit` argument can be `'min'`,
`'hours'`, or `'days'` to change the x-axis label and values.

## License

`pbpk-lite` is licensed under the MIT License.
