## Abaqus UMAT for Light Activated Shape Memory Polymer

### Single Element Verification
```
abaqus job=SingleExtension user=umat_LASMP_f77.for -interactive
```

#### Verification results
<div align=center>
<img src="https://github.com/brightfrank1999/abaqus-umat/blob/main/LASMP/imgs/StrainStressCurve.jpg" width="360" height="300"/><img src="https://github.com/brightfrank1999/abaqus-umat/blob/main/LASMP/imgs/ExtentofReaction.jpg" width="360" height="300">
</div>

### References
1. Modeling the mechanics of light activated shape memory polymers, JS Sodhi, IJ Rao, International Journal of Engineering Science 48 (11), 1576-1589
