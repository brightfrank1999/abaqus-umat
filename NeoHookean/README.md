## Abaqus UMAT for Neo-Hookean material

### Single Element Verification
```
abaqus job=SingleElemExt user=NeoHookean_f90.for -interactive
```

#### Verification results
<div align=center>
<img src="https://github.com/brightfrank1999/abaqus-umat/blob/main/NeoHookean/imgs/Verification.jpg">
</div>

### Inflation of a Hollow Cylinder
```
abaqus job=HollowCyliner user=NeoHookean_f90.for -interactive
```

#### Verification results
<div align=center>
<img src="https://github.com/brightfrank1999/abaqus-umat/blob/main/NeoHookean/imgs/HollowCylinder.png">
</div>
