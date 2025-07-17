# DAbI: Digital defocus Aberration Interference

 We discovered an optical phenomenon that the digitally summed Fourier spectrum of two images acquired from two-angle illumination exhibits interference-like fringe modulation when the sample is out-of-focus. These digital fringes correlate directly with defocus through a physics-based relation. Based on this principle, we developed an automatic, efficient, and generalizable defocus detection method termed digital defocus aberration interference (DAbI).

🌐 **[Project Page](https://hwzhou2020.github.io/DAbI-Web/)** 
📄 **[arXiv Paper](https://arxiv.org/abs/2507.10867)**
🧪 **[Data on OSF](https://osf.io/dvztc/)**  
💻 **[Code on GitHub](https://github.com/hwzhou2020/DAbI)**


Note: `.p` files are included for protected functions. They do not affect reproducibility. Full code will be released after publication.


## 📌 Publication Info

**Digital defocus aberration interference for automated optical microscopy**  
*Haowen Zhou\*, Shi Zhao\*, Yujie Fan, Zhenyu Dong, Oumeng Zhang, Viviana Gradinaru, Changhuei Yang*  
(\* Equal contribution)

- **arXiv**: https://arxiv.org/abs/2507.10867  
- **Project Page**: https://hwzhou2020.github.io/DAbI-Web/  
- **Code**: https://github.com/hwzhou2020/DAbI  
- **Data**: https://osf.io/dvztc/

---

## 🚀 Usage Instructions

### 🔬 Simulations

Ready to go — download and run the code directly.

Customize your system parameters in the main simulation scripts:

```
DAbI_simulation/DAbI_Simulation_2D.m
DAbI_simulation/DAbI_Simulation_3D.m
```

### 🧫 Experiments
Download experimental data from:
```
https://osf.io/dvztc/
```
Place the data in the following folders:
```
./DAbI_experiments_2D/Data/ 
./DAbI_experiments_3D/Data/ 
```
Then run:
```
DAbI_experiments_2D/DAbI_main_2D_experiments.m
DAbI_experiments_3D/DAbI_main_3D_experiments.m
```
Note: `.p` files are included for protected functions. They do not affect reproducibility. Full code will be released after publication.

---

## 📁 File Structure and explanation
```
├── DAbI_experiments_2D
│   ├── subFunctions
│   │   ├── add_aberration_zernike.m        # Ddd additional aberrations
│   │   ├── findDefocus_DAbI_Direction.p    # Get defocus direction
│   │   ├── findDefocus_DAbI_FFT.p          # FFT method
│   │   ├── findDefocus_DAbI.p              # DAbI function
│   │   ├── subPixelFit.m                   
│   │   └── system_parameters.mat
│   ├── DAbI_main_2D_experiments.m          # Main script for 2D thin sample experiments 
│   └── Data                                # Put 2D experimental raw data in this folder 
│
├── DAbI_experiments_3D
│   ├── subFunctions
│   │   ├── add_aberration_zernike.m
│   │   ├── findDefocus_DAbI_3D.p
│   │   ├── findDefocus_DAbI_Direction_3D.p
│   │   ├── findDefocus_DAbI_FFT_3D.p
│   │   ├── rowRangeCrop.m
│   │   ├── subPixelFit.m
│   │   └── system_parameters_3D.mat
│   ├── DAbI_main_3D_experiments.m          # Main script for 3D thick sample experiments
│   └── Data                                # Put 3D experimental raw data in this folder 
│
│
├── DAbI_simulation
│   ├── Simulated_Data                      # Folder to save intensity images
│   ├── Subfunctions_DAbI
│   │   ├── DAbI_2D
│   │   │   ├── findDefocus_DAbI_Direction.p
│   │   │   ├── findDefocus_DAbI_FFT.p
│   │   │   ├── findDefocus_DAbI.p
│   │   │   └── subPixelFit.m
│   │   └── DAbI_3D
│   │       ├── findDefocus_DAbI_3D.p
│   │       ├── findDefocus_DAbI_Direction_3D.p
│   │       ├── findDefocus_DAbI_FFT_3D.p
│   │       ├── rowRangeCrop.m
│   │       └── subPixelFit.m
│   ├── Subfunctions_Simulation
│   │   ├── add_aberration_zernike.m
│   │   ├── calBoundary.m
│   │   ├── cNeoAlbedo.m                    # Colormap
│   │   ├── imagingMultiSlice.m             # 3D simulation forward model
│   │   ├── NSCLC.mat                       # 2D simulation data non-small-cell lung cancer
│   │   ├── sam_1.mat                       # 3D simulation data part 1-5
│   │   ├── sam_2.mat
│   │   ├── sam_3.mat
│   │   ├── sam_4.mat
│   │   ├── sam_5.mat
│   │   └── USAF-pc200nm.png                # USAF 2D simulation data
│   ├── DAbI_Simulation_2D.m                # Main script for 2D simulation
│   └── DAbI_Simulation_3D.m                # Main script for 3D simulation
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
