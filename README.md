# pbpk-lite

`pbpk-lite` is a lightweight Python package for basic physiologically based pharmacokinetic (PBPK) modeling.

It provides a simple programmatic interface for defining substance properties, patient physiology, elimination kinetics, and solving the resulting ODE system.

## Features

- Substance partition coefficient calculation using logP and fraction unbound
- Patient blood flows and tissue volumes derived from body weight
- Linear liver and kidney elimination pathways
- ODE solution via `scipy.integrate.solve_ivp`
- Simple plotting helper for concentration visualization

## Installation

Install from PyPI:

```bash
pip install pbpk-lite
```

Install from source:

```bash
pip install .
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

t, c = m.simulate(doses, times)
print(t)
print(c.shape)
```

## Dosing Schedule

The `simulate()` method expects:

- `doses`: array-like of administered doses
- `times`: array-like of dosing times plus a final endpoint

Important: `times` must have exactly one more element than `doses`.
Each dose at index `i` is administered at `times[i]`, and the final value in `times` is the last observation or endpoint.

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

#### `simulate(doses, times)`

Simulate the PBPK model and return time points `t` and compartment concentrations `c`.

## License

`pbpk-lite` is licensed under the MIT License.
