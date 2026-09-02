# pbpk-lite

`pbpk-lite` is a lightweight Python package for basic physiologically based pharmacokinetic (PBPK) modeling.

It provides a simple programmatic interface for defining substance and metabolite properties, patient physiology, elimination kinetics, and solving the resulting ODE system.

The implementation uses nanograms (ng) for doses and amounts, milliliters (mL) for volumes, milliliters per minute (mL/min) for blood flows and clearances, nanograms per milliliter (ng/mL) for concentrations, and minutes for time. The simulation uses minutes internally, while the plotting helpers can display the x-axis in minutes, hours, or days via the optional `time_unit` argument.

## Features

- Substance and metabolite partition coefficient calculation using logP and fraction unbound
- Patient blood flows and tissue volumes derived from body weight
- Linear liver and kidney elimination pathways
- ODE solution via `scipy.integrate.solve_ivp`
- Plotting helpers for whole-model, venous-blood, and selected-compartment concentration profiles
- Support for different administration routes, including intravenous, intra-arterial, and inhalation dosing

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
	log_ps=(6.97, 5.33, 5.24),
	fus=(0.0022448, 0.01209, 0.101),
	mms=(314.469, 330.468, 344.451),
)
m.set_patient(bw=70)
m.set_elimination(cl_ls=(247, 998, 8.94), cl_ks=(0, 0, 1.48))

dose = 25e6  # ng

doses = [dose]
times = [0, 60*24]

t, c = m.simulate(doses, times, route_of_administration='iv')

m.graph_whole('concentrations.png')
m.graph_venous('venous.png', limit_of_detection=0.15, time_unit='hours')
m.graph_compartments(['liver', 'kidney'], 'selected.png', time_unit='days')

# The underlying simulation still uses minutes internally; only the displayed axis changes.
```

## Route of Administration

The `simulate()` method accepts a `route_of_administration` argument to control where each dose is introduced into the model.

Doses can only be administered to the parent substance. Metabolites are
produced by the model and cannot be dosed directly.

Supported values are:

- `iv`: intravenous dosing into the venous blood compartment (default)
- `ia`: intra-arterial dosing into the arterial blood compartment
- `inh`: inhalation dosing into the lung compartment

Example:

```python
m.simulate([dose], [0, 24*60], route_of_administration='inh')
```

## Dosing Schedule

The `simulate()` method expects:

- `doses`: array-like of administered doses
- `times`: array-like of dosing times plus a final endpoint

Important: `times` must have exactly one more element than `doses` and must be strictly increasing.
Each dose at index `i` is administered at `times[i]`, and the final value in `times` is the last observation or endpoint. Times are interpreted in minutes, so dosing schedules and simulation endpoints should be provided in minutes.

Example with one dose:

```python
# one dose at time 0, observation at 24 hours
doses = [dose]
times = [0, 60*24]
```

Example with two doses:

```python
# two identical doses spaced one hour apart, with a final observation at 24 hours
doses = [dose, dose]
times = [0, 60, 60*24]
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

#### `set_elimination(cl_ls, cl_ks)`

Set linear clearance from the liver and kidney compartments for each
substance.

- `cl_ls`: liver clearance values
- `cl_ks`: kidney clearance values

#### `simulate(doses, times, route_of_administration='iv')`

Simulate the PBPK model and return time points `t` and compartment concentrations `c`.

- `route_of_administration`: administration route used for each dose (`'iv'`, `'ia'`, or `'inh'`)
- administered doses enter the parent substance block; returned concentrations have 16 rows per substance

#### `graph_whole(name, time_unit='min')`

Save a multi-panel plot of concentrations across all compartments. The optional `time_unit` argument can be `'min'`, `'hours'`, or `'days'` to change the x-axis label and values.

#### `graph_venous(name, limit_of_detection=None, log=True, time_unit='min')`

Save a plot of venous blood concentrations, optionally marking a detection limit. The optional `log` argument controls whether the y-axis uses a logarithmic scale, and `time_unit` can be `'min'`, `'hours'`, or `'days'` to change the x-axis label and values.

#### `graph_compartments(compartments, name, time_unit='min')`

Save a plot of selected compartments by name or index. The optional `time_unit` argument can be `'min'`, `'hours'`, or `'days'` to change the x-axis label and values.

## License

`pbpk-lite` is licensed under the MIT License.
