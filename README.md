# pbpk-lite

`pbpk-lite` is a lightweight Python package for basic physiologically based pharmacokinetic (PBPK) modeling.

It provides a simple programmatic interface for defining substance properties, patient physiology, elimination kinetics, and solving the resulting ODE system.

The implementation is intended to be used with a consistent unit system based on mass units for doses and amounts, milliliters (mL) for volumes, and hours for time. When using the model, doses, volumes, flows, and clearances should all be expressed in compatible units within that framework.

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
m.set_elimination(cl_l=10, cl_k=0)

doses = [10]
times = [0, 24]

t, c = m.simulate(doses, times, route_of_administration='iv')
print(t)
print(c.shape)

m.graph_whole('concentrations.png')
m.graph_venous('venous.png', limit_of_detection=1.0)
m.graph_compartments(['liver', 'kidney'], 'selected.png')
```

## Route of Administration

The `simulate()` method accepts a `route_of_administration` argument to control where each dose is introduced into the model.

Supported values are:

- `iv`: intravenous dosing into the venous blood compartment (default)
- `ia`: intra-arterial dosing into the arterial blood compartment
- `inh`: inhalation dosing into the lung compartment

Example:

```python
m.simulate([10], [0, 24], route_of_administration='inh')
```

## Dosing Schedule

The `simulate()` method expects:

- `doses`: array-like of administered doses
- `times`: array-like of dosing times plus a final endpoint

Important: `times` must have exactly one more element than `doses`.
Each dose at index `i` is administered at `times[i]`, and the final value in `times` is the last observation or endpoint. Times are interpreted in hours, so dosing schedules and simulation endpoints should be provided in hours.

Example with one dose:

```python
# one dose at time 0, observation at 24 hours
doses = [10]
times = [0, 24]
```

Example with two doses:

```python
# doses at 0 and 12 hours, observation at 24 hours
doses = [10, 10]
times = [0, 12, 24]
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

#### `graph_whole(name)`

Save a multi-panel plot of concentrations across all compartments.

#### `graph_venous(name, limit_of_detection=None)`

Save a plot of venous blood concentrations, optionally marking a detection limit.

#### `graph_compartments(compartments, name)`

Save a plot of selected compartments by name or index.

## License

`pbpk-lite` is licensed under the MIT License.
