# CPA Parameters

This directory contains files with parameters for CPA (including gc-PC-SAFT).
The files named according to the pattern `NameYear.json` correspond to published parameters. 
The corresponding publication is provided in the [`literature.bib`](literature.bib) file.

```python
from reos.cpa import CPAParameters

params = CPAParameters.from_json(
    ["water", "carbon dioxide"], 
    "tsivintzelis2011.json", 
    "tsivintzelis2011_binary.json"
)
```

## Individual Parameters

|file|file binary|description|publication(s)|
|-|-|-|:-:|
[`kontogeorgis2010.json`](kontogeorgis2010.json) |   | water 4C, acids 1A, alcohols 2B, glycols 4C, alkanolamines 4C and non self-associating substances | [&#128279;](https://doi.org/10.1002/9780470747537.ch9)
[`kontogeorgis2006.json`](kontogeorgis2006.json) |   | water 4C and alcohols 3B | [&#128279;](https://doi.org/10.1002/9780470747537.ch9)
[`tsivintzelis2011.json`](tsivintzelis2011.json) | [`tsivintzelis2011_binary.json`](tsivintzelis2011_binary.json) | water, hydrocarbons and $\mathrm{CO_2}$ 1ea | [&#128279;](https://doi.org/10.1002/9780470747537.ch9)

