# DAbI: Digital defocus Aberration Interference for automated optical microscopy

 We observed that the digitally summed Fourier spectrum of two images acquired from two-angle illumination exhibits interference-like fringe modulation when the sample is defocused. These digital fringes correlate directly with defocus through a physics-based relation. Based on this principle, we developed an automatic, efficient, and generalizable defocus detection method termed digital defocus aberration interference (DAbI).
 
📑 **[Paper](https://doi.org/10.1038/s41467-026-72287-x)** 
🌐 **[Project Page](https://hwzhou2020.github.io/DAbI-Web/)** 
📄 **[arXiv Paper](https://arxiv.org/abs/2507.10867)**
🧪 **[Data on OSF](https://osf.io/dvztc/)**  
💻 **[Code on GitHub](https://github.com/hwzhou2020/DAbI)**


## 📌 Publication Info

**Digital defocus aberration interference for automated optical microscopy**  
*Haowen Zhou\*, Shi Zhao\*, Yujie Fan, Zhenyu Dong, Oumeng Zhang, Viviana Gradinaru, Changhuei Yang*  
(\* Equal contribution)

- **Paper**: https://doi.org/10.1038/s41467-026-72287-x
- **arXiv**: https://arxiv.org/abs/2507.10867  
- **Project Page**: https://hwzhou2020.github.io/DAbI-Web/  
- **Code**: https://github.com/hwzhou2020/DAbI  
- **Data**: https://osf.io/dvztc/

---

## Prerequisites

- MATLAB
- Image Processing Toolbox
- Parallel Computing Toolbox (optional, only needed for GPU acceleration)

The MATLAB scripts use relative paths. Run each script from its own folder rather than from the repository root.

## 🚀 Usage Instructions

### 🔬 Simulations

The simulation examples are self-contained. This repository already includes the bundled sample assets:

- `NSCLC.mat` and `USAF-pc200nm.png` for 2D simulation
- `cell.mat` for 3D simulation

`NSCLC.mat` and `cell.mat` are in the Subfucntions_Simulation folder in "https://osf.io/dvztc/". 

Please move them to DAbI_simulation/Subfunctions_Simulation folder in this repository.

Open MATLAB, change into the simulation folder, then run either script:

```
cd DAbI_simulation
DAbI_Simulation_2D
% or
DAbI_Simulation_3D
```

Customize system parameters directly in:

```
DAbI_simulation/DAbI_Simulation_2D.m
DAbI_simulation/DAbI_Simulation_3D.m
```

### 🧫 Experiments
Download experimental data from:
```
https://osf.io/dvztc/
```
Create the following folders if they do not already exist, then place the downloaded data there:
```
./DAbI_experiments_2D/Data/ 
./DAbI_experiments_3D/Data/ 
```

Then open MATLAB, change into the corresponding experiment folder, and run:
```
cd DAbI_experiments_2D
DAbI_main_2D_experiments

cd ../DAbI_experiments_3D
DAbI_main_3D_experiments
```

---

## 📁 File Structure and Explanation
```
├── DAbI_experiments_2D
│   ├── subFunctions
│   │   ├── add_aberration_zernike.m        # Add additional aberrations
│   │   ├── Data_generator.m
│   │   ├── findDefocus_DAbI_Direction.m    # Get defocus direction
│   │   ├── findDefocus_DAbI_FFT.m          # FFT-based method
│   │   ├── findDefocus_DAbI.m              # Main DAbI function
│   │   ├── subPixelFit.m
│   │   └── system_parameters.mat
│   ├── DAbI_main_2D_experiments.m          # Main script for 2D thin sample experiments 
│   └── Data                                # Put 2D experimental raw data here
│
├── DAbI_experiments_3D
│   ├── subFunctions
│   │   ├── add_aberration_zernike.m
│   │   ├── findDefocus_DAbI_3D.m
│   │   ├── findDefocus_DAbI_Direction_3D.m
│   │   ├── findDefocus_DAbI_FFT_3D.m
│   │   ├── rowRangeCrop.m
│   │   ├── subPixelFit.m
│   │   └── system_parameters_3D.mat
│   ├── DAbI_main_3D_experiments.m          # Main script for 3D thick sample experiments
│   └── Data                                # Put 3D experimental raw data here
│
├── DAbI_simulation
│   ├── Subfunctions_DAbI
│   │   ├── DAbI_2D
│   │   │   ├── findDefocus_DAbI_Direction.m
│   │   │   ├── findDefocus_DAbI_FFT.m
│   │   │   ├── findDefocus_DAbI.m
│   │   │   └── subPixelFit.m
│   │   └── DAbI_3D
│   │       ├── findDefocus_DAbI_3D.m
│   │       ├── findDefocus_DAbI_Direction_3D.m
│   │       ├── findDefocus_DAbI_FFT_3D.m
│   │       ├── rowRangeCrop.m
│   │       └── subPixelFit.m
│   ├── Subfunctions_Simulation
│   │   ├── add_aberration_zernike.m
│   │   ├── calBoundary.m
│   │   ├── cNeoAlbedo.m                    # Colormap
│   │   ├── imagingMultiSlice.m             # 3D simulation forward model
│   │   ├── cell.mat                        # 3D simulation sample
│   │   ├── NSCLC.mat                       # 2D simulation sample
│   │   └── USAF-pc200nm.png                # USAF 2D simulation data
│   ├── DAbI_Simulation_2D.m                # Main script for 2D simulation
│   └── DAbI_Simulation_3D.m                # Main script for 3D simulation
│
├── DAbI_simulation/Simulated_Data          # Created automatically by 2D simulation when saving is enabled
├── DAbI_simulation/Simulated_Data_3D       # Created automatically by 3D simulation when saving is enabled
├── LICENSE
└── README.md
```

---
## 🤝 Contact & Collaboration
For questions, feedback, or collaboration inquiries, feel free to:

Open a GitHub issue: https://github.com/hwzhou2020/DAbI/issues

Contact via both emails:

hzhou7@caltech.edu

szhao5@caltech.edu

## 📖 Citation
If you find this work useful, please cite:
```
@article{Zhou2026DigitalDefocus,
  author  = {Zhou, Haowen and Zhao, Shi and Fan, Yujie and Dong, Zhenyu and Zhang, Oumeng and Gradinaru, Viviana and Yang, Changhuei},
  title   = {Digital defocus aberration interference for automated optical microscopy},
  journal = {Nature Communications},
  year    = {2026},
  date    = {2026-04-24}
}
```

```
@misc{zhou2025DAbI,
  title={Digital defocus aberration interference for automated optical microscopy}, 
  author={Haowen Zhou and Shi Zhao and Yujie Fan and Zhenyu Dong and Oumeng Zhang and Viviana Gradinaru and Changhuei Yang},
  year={2025},
  eprint={2507.10867},
  archivePrefix={arXiv},
  primaryClass={physics.optics},
  url={https://arxiv.org/abs/2507.10867},
}
```
---
## ⚖️ License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.  
See [`LICENSE`](./LICENSE) for details.

© 2025 Biophotonics Laboratory, Caltech
