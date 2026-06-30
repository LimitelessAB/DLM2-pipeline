<img width="125" height="109" alt="image" src="https://github.com/user-attachments/assets/a08e1fa2-5453-4e94-ae57-cd2936397581" />
# DLM2-pipeline (DustPy Linked to MCFOST and Dark Lane Measurements pipeline)
### Attention! Please note this code requires the user to have already installed DustPy python module, PyMCFOST python module, and MCFOST!
- DustPy: <https://github.com/stammler/dustpy.git>
- MCFOST: <https://github.com/cpinte/mcfost.git>
- PyMCFOST: <https://github.com/cpinte/pymcfost.git>
## About:
DLM2-pipeline is a tool that assembles into one function call: 
1. Launching DustPy
2. Converting results for MCFOST
3. Launching MCFOST
4. Convolving resulting images with instrument PSF
5. Measuring dark lanes and radius following a method described in Duchêne et al. 2024: <https://iopscience.iop.org/article/10.3847/1538-3881/acf9a7>.

## Installing:
Run into your terminal : $pip install git+https://github.com/LimitelessAB/DLM2-pipeline.git
