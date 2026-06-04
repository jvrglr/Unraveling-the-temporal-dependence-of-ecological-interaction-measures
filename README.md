# Unraveling the temporal dependence of ecological interaction measures
Codes and datasets used in Ref. [1].

## Codes

* **Grasshopper.py**: Compute interactions from grasshopper abundance data (Figs.2 and A.1 in Ref. [1]).
* **synthetic_data**: Fortran codes to generate data through integration of consumer-resource model and compute measures of interaction from it.
  * **M_declarations.f90**:  Define public variables to be used in the rest of modules. Parameters and variables are defined in Eqs.  (9) and (10) of Ref.[1].
  * **M_functions.f90**: Mathematical functions used in the main codes. Here resource renewal and consumer consumption functions are defined.
  * **M_subroutines.f90**: Main script containing codes for numerical integration of consumer-resource model and computation of interaction measures from trajectories.
  * **GH_main_inference.f90**: Example of code to generate synthetic data and compute interactions from it used to generate Fig. 7 in Ref.[1].
  
## Data
Data rights belong to the parties responsible for data acquisition. 

Any use or distribution of the data should be accompanied by citation of the original data sources (see below) and discussed with the parties responsible for data acquisition.

* **Grasshopper data**: Grasshopper abundances over time extracted from Fig.1 in ref. [2].

## License
This project is shared for **academic and research purposes**. 

The codes are free to use, redistribute, modify, and share for research purposes, provided that proper credit is given to the authors through citation of [1].

## Citation

```
@article{8fb1-jmqh,
  title = {Unraveling the temporal dependence of ecological interaction measures},
  author = {Aguilar, Javier and Suweis, Samir and Maritan, Amos and Azaele, Sandro},
  journal = {PRX Life},
  pages = {},
  year = {2026},
  month = {Jun},
  publisher = {American Physical Society},
  doi = {10.1103/8fb1-jmqh},
  url = {https://link.aps.org/doi/10.1103/8fb1-jmqh}
}
```

## References

[1] Aguilar J., Suweis S., Maritan A., Azaele S. Phys. Rev. X Life (2026). [10.1103/8fb1-jmqh](https://journals.aps.org/prxlife/accepted/10.1103/8fb1-jmqh)

[2] Ritchie, M. E., & Tilman, D. (1993). Oecologia, 94(4), 516-527. [10.1007/BF00566967](https://doi.org/10.1007/BF00566967)
