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
```@docs
identifier
format
element
functional
valence_charge
is_norm_conserving
is_ultrasoft
is_paw
has_spin_orbit
has_model_core_charge_density
formalism
```
