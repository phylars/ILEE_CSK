# For users without Python experience

ILEE_CSK is a python library and therefore we need to setup your python environment to process your data. Don't be afraid, should not be difficult.

First, we need to install your python environment. For novice, we only recommend the one-step simplest way: [Anaconda](https://www.anaconda.com/). Please find the Anaconda installer that applies to your computer at https://www.anaconda.com/products/individual, and install it. Anaconda has a python environment has most of the libraries we need, so you do not need to install other libraries. 

Anaconda has a well-established Python IDE (Integrated development environment) [Spyder](https://www.spyder-ide.org/) optimized for scientific computation and visualization. Please just regard it as a powerful calculator in your computer. When you installed our library,  you need to open Spyder to compute your samples. You can either search Spyder in your apps or open Anaconda Navigator to open it. Typically it looks like this:

![spyder.png](https://github.com/phylars/image_host/blob/main/ILEE/spyder.png?raw=true)

Since you need to install ILEE_CSK, you need to use commend lines. <u>Please note that you should **NOT** use your Windows or Apple command lines, but instead, please open Anaconda Prompt and enter the install commands you wants.</u> This makes your intended libraries installed in Anaconda environment particularly. 

You can either search Anaconda Prompt in your apps or open Anaconda Navigator to open it. It should look like this:

![cmd.png](https://github.com/phylars/image_host/blob/main/ILEE/cmd.png?raw=true)

![cmd1.png](https://github.com/phylars/image_host/blob/main/ILEE/cmd1.png?raw=true)

This is the interface where you can install any library. We use use it for library installing purposes, which should not be confounded with Spyder.

# Installing ILEE_CSK

ILEE_CSK is available at pipy...

```sh
pip install ILEE_CSK
```

ILEE_CSK supports MATLAB-based GPU acceleration for native-3D image processing. If you have MATLAB license as well as CUDA GPU, and wish to proceed native-3D mode. Please download......

# Installing napari for 3D visualization

The [napari](https://napari.org/) is a multi-dimensional image viewer for python. It is a outstanding tool for confocal image visualization, especially for visualizing 3D dnapariata structure. We just recommend you directly follow their documentation to install it: https://napari.org/.

# Cytoskeleton indices computation

## Batch processing all samples in a document

Confocal images of cytoskeletons are generally taken as 3D stacks (Native 3D structure). Hereby we support batch processing of files of  3D stack. If your source images are already 2D, please check our API to process it individually.

### Generating TIFF raw images by FIJI (ImageJ). 

[FIJI](https://fiji.sc/) is a java-based free-to-use image processing GUI software with optimized for scientific use. If you don't have it, please download it and congratulations -- you have found a great treasure. Put all of your confocal image files in the same document. Open FIJI, select process->batch->convert

![fiji_batch.png](https://github.com/phylars/image_host/blob/main/ILEE/fiji_batch.png?raw=true)

Input is the document of your raw confocal-format files; output is wherever you want. Press convert, and we just get what we need for ILEE. Please note that output format MUST BE tiff, interpolation MUST be None, and scale factor MUST BE 1. This make sure your tiff is literally raw data. 

### Find your channel index and K2

You need to learn the data structure of the raw TIFF files you generated to finish parameter setting. Please get the path of one of your generated TIFF file and read it by skimage.io, such as:

```python
image_example_path = 'D:/my_experiments/generated_TIFF_images/sample_001.tiff'
from skimage import io
img = io.imread(image_example_path)
```

Generally, img should be a 4D array, and one of the dimension with the smallest size among the four dimensions should be your fluorescence channel. (For Python novice, you can check the Spyder variable explorer on the column "size".) For example, when I use a Olympus confocal system and I gets a array with dimension size of (25, 3, 800, 800), which means (z, channel, x, y). Difference confocal system have these numbers occur at different order. Since you know how did you take the image by confocal, you can understand what these number means easily. <u>NOTE: Our batch processing pipeline requires your channel to be the dimension with the lowest size, and the second lowest to be z-axis resolution.</u> This applies to almost all of the cytoskeleton confocal images, because the number of channels is generally no greater than 5, with maximum 4 fluorescence channels plus a transmission channel, and z-stack number is usually greater than 5 to cover the total cell. If you have a abnormal data structure that does not apply, we wish to you inform us about what the structure is and why it is necessary, so that we may improve our algorithm. 

You need to figure out which of your channel is the cytoskeleton fluorescence channel you want to analyze. Generally, the channel order will be the same as your confocal setting, and be aware that Python index start from 0 instead of 1. For example, if I used channel 1 to see GFP-labeled actin, channel 2 to see chloroplast far-red light, and an additional transmission channel to see bright field, my objective index is supposed to be 0 (for channel 1). However, to clearly verify which one is your objective channel, you can try visualize their projected image:

```python
import matplotlib.pyplot as plt
img = img[:,0,:,:] #You can change this "0" to the other channel indexes you want to check.
img = np.amax(img, axis=0) #In my case, my first dimension (index 0) is z-axis, if your z-axis is another dimension, please change accordingly.
#run this three line together:
plt.figure(figsize = (18, 18))
plt.imshow(img)
plt.show()
```

Another way is directly see the 3D voxel image by napari:

```python
import napari
with napari.gui_qt():
    viewer = napari.Viewer()
    for i in range(3): #This 3 is the total number of your channels. Change accordingly.
        viewer.add_image(img[:,i,:,:])   
```

 

When you learn your data structure, you can start to evaluate your K2 for the batch using the function **opt_k2**. For example, if my index of cytoskeleton fluorescence channel is 0, you can do like such:

```python
import ILEE_CSK
folder = 'D:/my_experiments/generated_TIFF_images/'
the_optimal_K2 = ILEE_CSK.opt_k2(folder_path = folder, target_channel = 0)
```

### Batch processing by 2D mode

After you determine the optimal K2 for you sample pool. Using the same fluorescence channel index, you can finish batch processing by just one step:

```python
results = ILEE_CSK.analyze_document_2D (folder_path = folder, obj_channel = 0, k2 = the_optimal_K2)
```

Please read the API to learn about the parameters. 

The "results" will be a table in format of Pandas dataframe, which you can check and process within python IDE. You may also want to directly output it as an Excel file, which can be done as such:

```python
result.to_excel('D:/my_experiments/generated_TIFF_images/analysis_result.xlsx')
```

### Batch processing by 3D mode on CPU

Batch processing on 3D mode use a similar pipeline function that is easy to use. We recommend no less than 16G memory, otherwise it may fail. You need to manually check your voxel size on x/y- and z-axis by μm, if you don't remember. A simple way is to check by FIJI (ImageJ); after open one of you samples, select image->Show Info to get it. In the example above, my voxel size is 0.132μm on x/y-axis and 0.5μm on z-axis. 

![fiji_info.png](https://github.com/phylars/image_host/blob/main/ILEE/fiji_info.png?raw=true)

You can use batch processing function on 3D mode. I will use my image information as an example:

```python
results = ILEE_CSK.analyze_document_3D (folder_path = folder,
                                        obj_channel = 0, 
                                        k2 = the_optimal_K2, 
                                        xy_unit = 0.132, 
                                        z_unit = 0.5,     
                                        pixel_size = 0.132, 
                                        single_k = True, 
                                        use_GPU = False)
```

The setting above uses CPU mode, which is very slow (5~7min per 25x800x800 image). Therefore, we highly recommend to turn on single_k mode to save half of the time with acceptable sacrifice at accuracy on thick filaments.

### Batch processing by 3D mode on GPU

