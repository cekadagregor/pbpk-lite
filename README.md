# pbpk-lite

`pbpk-lite` is a lightweight Python package for basic physiologically based pharmacokinetic (PBPK) modeling.

It provides a simple programmatic interface for defining substance properties, patient physiology, elimination kinetics, and solving the resulting ODE system.

The implementation is intended to be used with a consistent unit system based on mass units for doses and amounts, milliliters (mL) for volumes, and minutes for time. When using the model, doses and clearances should be expressed in compatible units within that framework. The simulation itself uses minutes internally, while the plotting helpers can display the x-axis in minutes, hours, or days via the optional `time_unit` argument.

## Features

- Substance partition coefficient calculation using logP and fraction unbound
- Patient blood flows and tissue volumes derived from body weight
- Linear liver and kidney elimination pathways
- ODE solution via `scipy.integrate.solve_ivp`
- Plotting helpers for whole-model, venous-blood, and selected-compartment concentration profiles
- Support for different administration routes, including intravenous, intra-arterial, and inhalation dosing

## Installation

Install from PyPI:

```bash
pip install pbpk-lite
```

## Quick Start

```python
from pbpk_lite import model

m = model()
m.set_substance(log_p=6.97, fu=0.0022448)
m.set_patient(bw=70)
m.set_elimination(cl_l=748.6643482986731, cl_k=0)

dose = 25e6

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

#### `set_substance(log_p, fu)`

Set the substance physicochemical properties.

- `log_p`: log octanol-water partition coefficient
- `fu`: fraction unbound in blood

#### `set_patient(bw)`

Set patient physiological parameters using body weight in kilograms.

#### `set_elimination(cl_l=0, cl_k=0)`

Set linear clearance from liver and kidney compartments.

#### `simulate(doses, times, route_of_administration='iv')`

Simulate the PBPK model and return time points `t` and compartment concentrations `c`.

- `route_of_administration`: administration route used for each dose (`'iv'`, `'ia'`, or `'inh'`)

#### `graph_whole(name, time_unit='min')`

Save a multi-panel plot of concentrations across all compartments. The optional `time_unit` argument can be `'min'`, `'hours'`, or `'days'` to change the x-axis label and values.

#### `graph_venous(name, limit_of_detection=None, log=True, time_unit='min')`

Save a plot of venous blood concentrations, optionally marking a detection limit. The optional `log` argument controls whether the y-axis uses a logarithmic scale, and `time_unit` can be `'min'`, `'hours'`, or `'days'` to change the x-axis label and values.

#### `graph_compartments(compartments, name, time_unit='min')`

Save a plot of selected compartments by name or index. The optional `time_unit` argument can be `'min'`, `'hours'`, or `'days'` to change the x-axis label and values.

## License

`pbpk-lite` is licensed under the MIT License.
