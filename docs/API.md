# API

------

## General batch-processing pipeline

<br/>

**opt_k2 (folder_path, target_channel = None):**

This functions calculates the optimal K2 for a batch of samples (raw TIFF files) to be analyzed. Please include all samples from different groups to be compared into the same document.

Parameters:

folder_path : string, a local path of your image samples. The final document name should contain an "/" in the end, for example, folder_path = 'D:/experiment/all my samples/'

target_channel : int, optional. <u>It should be an provided int for most of the situations.</u> If your raw TIFF file does not have the channel dimension, please use default.

Returns:

K2: float, your estimated optimal K2 for a total batch of sample for comparison

<br/>

**analyze_document_2D (folder_path, obj_channel, k2, pixel_size = 1):**

Process and analyze all the raw tiff samples in your document by 2D mode and output the result as a table (dataframe). All images are automatically transformed to 12 bit to be compatible. In 2D mode, K1 = 2.5.

Parameters:

folder_path : string, a local path of your image samples. The final document name should contain an "/" in the end, for example, folder_path = 'D:/experiment/all my samples/'

obj_channel : int. Your channel index of cytoskeleton fluoresence in the 4D array read from raw TIFF

k2 : number. The universal K2 you set. Recommend use the optimal K2 calculated by the function opt_k2

pixel_size : number, optional. The physical lenghth unit of your pixel, by micrometer (μm). The default is 1. If you set it into the real pixel size, the unit "PU" in your output table should be considered as μm.

Returns:

Pandas dataframe, all cytoskeleton indices. You can process it inside python IDE or export it as excel file. 

<br/>

**analyze_document_3D (folder_path, obj_channel, k2, xy_unit, z_unit, pixel_size = 1, single_k = True, use_GPU = False):**

Process and analyze all the raw tiff samples in your document in 3D mode and output the result as a table (dataframe). All images are automatically transformed to 12 bit to be compatible. Please note that it can be slow if the MATLAB-based GPU acceleration is not available. In 3D mode, K1 = 10^((log10(2.5)+log10(K2))/2).

Parameters

folder_path : string, a local path of your image samples. The final document name should contain an "/" in the end, for example, folder_path = 'D:/experiment/all my samples/'

obj_channel : int, your channel index of cytoskeleton fluoresence in the 4D array read from raw TIFF.

k2 : number, the universal K2 you set. Recommend use the optimal K2 calculated by the function opt_k2.

xy_unit : number, the unit size of x- and y-axis by μm. Please check using ImageJ.

z_unit : number, the unit size of z-axis by μm, which is also your step size of stack imaging. Please check using ImageJ.

pixel_size : number, optional. The physical lenghth unit (on x- and y-axis) of your voxel, by micrometer (μm). The default is 1. If you set it into the real pixel size, it should be equal to the parameter xy_unit, and the unit "PU" in your output table should be considered as μm.

single_k : bool, optional. Whether to use single-K mode (use only K1) to save time. The default is True.

use_GPU : bool, optional. Whether to activate GPU acceleration. We highly recommend you turn it on if you have MATLAB and compatible GPU. The default is False.

Returns

pandas dataframe, all cytoskeleton indices. You can process it inside python IDE or export it as excel file.

<br/>

------

## Processing individual sample

<br/>

**rearrange_image_dimension (img, target_channel = None):**

Transform your image array of any demension order (most likely x(height)-y(width)-z(stack)-channel or channel-x-y-z structure, depending on your confocal microscope) into z-x-y array with your object channel.This function depends on sorting, so your channel number must be lower than your stack number, and your stack number (z-resolution) must be lower than your x and y resolution.

Parameters

img : 3D or 4D array, the image read from a tiff file by skimage.io

target_channel : int, your cytoskeleton fluorescence channel index. If img is 3D array (only one channel, which is your cytoskeleton), target_channel must be None.

Returns

3D array, the image of your object channel extracted from TIFF-generated 4D array, <u>in z-x-y format</u>. 

<br/>

**ILEE_2d (img, k2, k1 = 2.5, pL_type = 'pL_8', gauss_dif = True):**

In 2D mode, generate a difference image (raw_image - threshold_image) as ILEE output for downstream analysis or visualization.

Parameters

img : float 2D array. An image of z-axis maximum projection of your target channel or the 2D "raw" image. Values should be transformed into float type with 0-4095 dynamic range. For example, if you have a 12 bit int-type image, as a variable img, then it should be used as img.astype('float'); if you have a 16 bit int-type image, then you should transform img into img2 = img.astype('float')/65535*4095.

k2 : a number; should be your pre-determined optimal K2 for comparison purpose. 

k1 : The default is 2.5, which can identify one-pixel-size filament.

pL_type : string, Laplacian matrix type. The default and recommended is 'pL_8'. 

gauss_dif : Bool, default is True; whether the difference image is slightly gaussian-blurred for smoother render.

Returns

float 2D array, the rendered 2D difference image. 

<br/>

**analyze_actin_2d_standard (img, img_dif, pixel_size = 1):**

Compute all cytoskeleton indices using 2d raw and difference image as input.

Parameters

img : float 2D array. The input 2D image. 

img_dif : float 2D array. The difference image (raw_image - threshold_image) generated by ILEE_2d.

pixel_size : number, optional. The physical lenghth unit of your pixel, by micrometer (μm). The default is 1. If you set it into the real pixel size, the unit "PU" in your output table should be considered as μm.

Returns

pandas.DataFrame, all cytoskeleton indices of input sample: occupancy, linear_density, skewness, cv, Diameter_tdt, Diameter_sdt, sev_act, branching_act, anisotropy.

<br/>

**ILEE_3d (img, xy_unit, z_unit, k1, k2, single_k1 = False, use_matlab = False, use_matlabGPU = False, gauss_dif = True):**

In 3D mode, generate a difference image (raw_image - threshold_image) as ILEE output for downstream analysis or visualization.

Parameters

img : 3D array. The input 3D image. Values should be transformed into float type with 0-4095 dynamic range. For example, if you have a 12 bit int-type image, as a variable img, then it should be used as img.astype('float'); if you have a 16 bit int-type image, then you should transform img into img2 = img.astype('float')/65535*4095. 

xy_unit : number. The unit size of voxel on x- and y-axis by μm. Please check using ImageJ.

z_unit : number. The unit size of voxel on z-axis by μm, which is also your step size of stack imaging. Please check using ImageJ.

k1 : number. ILEE K1. for 3D mode we recommend K1 = 10^((log10(2.5)+log10(K2))/2).

k2 : number. ILEE K2.

single_k1 : bool, optional. Whether to use single-K mode (use only K1) to save time. The default is False.

use_matlab : bool, optional. Whether to use MATLAB for ILEE. The default is False. 

use_matlabGPU : bool, optional. Whether to activate GPU_acceleration in MATLAB. The default is False. If use_matlab is False, this makes no difference. 

gauss_dif : Bool, optional. Whether the difference image is gaussian-blurred slightly for smoother render. Default is True. 

Returns

float 3D array, the rendered 3D difference image, the same shape as input.  

<br/>

**analyze_actin_3d_standard (img, img_dif_ori, xy_unit, z_unit, oversampling_for_bundle = True, pixel_size = 1):**

Compute all cytoskeleton indices using 3d raw and difference image as input.

Parameters

img : 3D array. Input 3D image

img_dif_ori : 3D array. The difference image to define the cytoskeleton components.

xy_unit : number. The unit size of voxel on x- and y-axis by μm. Please check using ImageJ.

z_unit : number. The unit size of voxel on z-axis by μm, which is also your step size of stack imaging. Please check using ImageJ.

oversampling_for_bundle : bool, optional. Whether to use oversampled image to calculated diameter indices. Recommend to turn on for accuracy. The default is True.

pixel_size : number, optional. The physical lenghth unit (on x- and y-axis) of your voxel, by micrometer (μm). The default is 1. If you set it into the real pixel size, it should be equal to the parameter xy_unit, and the unit "PU" in your output table should be considered as μm.

Returns
-------
pandas.DataFrame, all cytoskeleton indices of input sample: occupancy, linear_density, skewness, cv, Diameter_tdt, Diameter_sdt, sev_act, branching_act, anisotropy.

