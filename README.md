ILEE_CSK is a Python library that provides a platform for unguided and automated quantitative analysis of cytoskeletal images at high accuracy and free from bias, supporting both 2D and native 3D data structure. ILEE means Implicit Laplacian of Enhanced Edge, which is a local thresholding algorithm that detects and utilize the edge of cytoskeleton filaments as the local baseline for super-fine segmentation. Our robust algorithm can compute tens of indices reflecting diverse quantitative features of eukaryotic cytoskeleton for biological interpretation, which include (de/)polymerization, bundling, severing, branching, and local directionality. 

<br/>

For fast and convenient access to ILEE_CSK for everyone (without requirement for coding experience), please use our Google Colab pipeline:

(2D mode)
<br/>
[https://colab.research.google.com/github/phylars/ILEE_CSK/blob/ipynb/Copy_of_ILEE_2D_mode_(V_1_10).ipynb](https://colab.research.google.com/github/phylars/ILEE_CSK/blob/ipynb/ILEE_2D_mode_(V_1_10).ipynb)
<br/>
(3D CPU mode)
<br/>
https://colab.research.google.com/github/phylars/ILEE_CSK/blob/ipynb/ILEE_3D_mode_(V_1_0).ipynb
<br/>
(3D GPU mode) Google Colab currently not supported. Please use the library locally for Matlab empowered 3D GPU mode.


<br/>
ATTENTION for users who are handling cell samples with true blank (exactly nothing inside, not including tissue, such as leaf epidermis, or extremely maximized single cell with full field occupied) area in the image:
Your effective area (area of the cell) is very important to the accuracy of the result. While we do have the automatic cell segment algorithm, but we cannot promise the area is absolutely accurate. To obtain accurate result, please first run the pipeline 

<br/>
Github:

https://github.com/phylars/ILEE_CSK

Documentation:

https://phylars.github.io/ILEE_CSK/

