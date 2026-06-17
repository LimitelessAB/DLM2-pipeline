# DLM2-pipeline (DustPy Linked to MCFOST and Dark Lane Measurements pipeline)
### Attention!!! Please note this code requires the user to have already installed DustPy python module, PyMCFOST python module, and MCFOST!
- DustPy: <https://github.com/stammler/dustpy.git>
- MCFOST: <https://github.com/cpinte/mcfost.git>
- PyMCFOST: <https://github.com/cpinte/pymcfost.git>
## About:
DLM2-pipeline is a tool that assembles into one function call: launching DustPy, converting results for MCFOST, launching MCFOST, convolving resulting images with instrument PSF, and measuring dark lanes and radius following a method described in Duchêne et al. 2024: <https://iopscience.iop.org/article/10.3847/1538-3881/acf9a7>.
