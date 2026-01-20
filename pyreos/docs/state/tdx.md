## A Thermodynamic State at (T,ρ,x)

Units in SI.

Create a State object, given `temperature`, `density` and `composition`.
If the conditions passed are physically wrong, an error will be raised 
when trying to compute pressure.

### Parameters
----------

- `eos` : EquationOfState
- `temperature` : float
- `pressure` : float
- `x` : numpy.ndarray[float]

### Returns
-------
State at (T,ρ,x)