## Abaqus UMAT for Arruda-Boyce material

### Single Element Verificaiton
```
abaqus job=SingleElemExt user=ArrudaBoyce_f77.for -interactive
abaqus job=HollowCylinder user=ArrudaBoyce_f77.for cpus=2 -interactive
```

#### Verification results
![image](https://github.com/brightfrank1999/abaqus-umat/blob/main/ArrudaBoyce/imgs/Verification.jpg)

