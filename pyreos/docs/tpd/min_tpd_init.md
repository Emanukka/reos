## Min TPD initialization

Create a MinTPD object from the mother phase state at '(T,P,z)' condition.
After the creation, two methods are available:

- `dg()`: returns the minimum TPD value (ΔG formation of the incipient phase from mother phase)
- `state`: returns the incipient phase State at (T,P,x)

See [Michelsen, 1982](https://doi.org/10.1016/0378-3812(82)85001-2)

### Parameters
----------

- `xphase` : 'liquid' or 'vapor'
- `xguess` : numpy.ndarray[float]
- `tol` : float, optional
- `it_max` : int, optional

### Returns
-------
MinTPD 
