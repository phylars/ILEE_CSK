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
<br/>

<br/>
ATTENTION for users who are handling cell samples with true blank (exactly nothing inside, not including tissue, such as leaf epidermis, or extremely maximized single cell with full field occupied) area in the image:
<br/>
<br/>
Your effective area (area of the cell) is very important to the accuracy of the result. While we do have the automatic cell segment algorithm, but we cannot promise the area is absolutely accurate. First please find a channel where the total cell area are most contrast to blank area. Please understand I am talking about a "dirty" channel or bright field (we tested animal cell actin and it works); highly specific dye targeting small spots of proteins is a bad choice for cell coverage. To obtain accurate results, please use your "dirtiest" channel to run the pipeline to test/optimize whether the total cell segmentation is accurate:
<br/>
[https://colab.research.google.com/github/phylars/ILEE_CSK/blob/ipynb/2D_mode_calculating_effective_area.ipynb]
<br/>
<br/>
If it looks approximately ideal, you can use the general pipeline without any concern. IF NOT, we specifically generate another pipeline for manual input of effective area, which several recommend choices:
<br/>
[https://github.com/phylars/ILEE_CSK/blob/ipynb/ilee_2d_mode_with_custom_effective_area.ipynb]
<br/>
<br/>
You may:
<br/>
(1) use the test/optimize Colab pipeline to determine a your preferred parameter setting and generate the area excel file in the end, as input of ILEE 2D pipeline with manual input of effective area.
<br/>
(2) mixing measured areas of different parameters from out pipeline by editting the excel file, but please maintain the excel structure. 
<br/>
(3) If you understand the mechanism here, you can find your own way to measure effective area (such as ImageJ/FIJI, etc.).
<br/>
<br/>
NOTE: please make sure the order of images samples in the area excel file is NEVER changed, and this order MUST be same as your final result output (otherwise your data does not make sense). We do not check whether these two orders match during the pipeline, so you have to be cautious.
<br/>

<br/>
Github:

https://github.com/phylars/ILEE_CSK

Documentation:

https://phylars.github.io/ILEE_CSK/

