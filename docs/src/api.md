# API reference

## Index

```@index
Pages = ["api.md"]
```

## Types

```@docs
PsPFile
UpfFile
Psp8File
HghFile
```

## Functions

### Loading and saving
```@docs
load_psp_file
save_psp_file
```

### Conversion
```@docs
UpfFile(::Psp8File)
```

### Metadata
Note that these functions are not exported, but part of the public API.

```@docs
PseudoPotentialIO.identifier
PseudoPotentialIO.format
PseudoPotentialIO.element
PseudoPotentialIO.functional
PseudoPotentialIO.valence_charge
PseudoPotentialIO.is_norm_conserving
PseudoPotentialIO.is_ultrasoft
PseudoPotentialIO.is_paw
PseudoPotentialIO.has_spin_orbit
PseudoPotentialIO.has_model_core_charge_density
PseudoPotentialIO.formalism
```
