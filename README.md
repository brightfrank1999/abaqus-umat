# abaqus-umat

## Requirements
Abaqus 2021  
Visual Studio Community 2019  
Intel oneAPI 2023  

## System Information
<div align=center>
<img src="https://github.com/brightfrank1999/abaqus-umat/blob/main/imgs/SystemInfo.jpg">
</div>
  
## Fortran 90 free format
Some UMAT subroutines are written in Fortran 90 free format. For these subroutines, we have to modify the environmental variable of the compiler. In "C:\SIMULIA\EstProducts\2021\win_b64\SMA\site\win86_64.env", add the term '/free' in "compile_fortran".
<div align=center>
<img src="https://github.com/brightfrank1999/abaqus-umat/blob/main/imgs/win_env.jpg">
</div>

## STOP Function
