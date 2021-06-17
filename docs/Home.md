# ILEE_CSK

The eukaryotic cytoskeleton plays essential roles in cell signaling, trafficking, and motion. We aim to provide the best cytoskeleton segmentation-analysis algorithm and platform. 

ILEE_CSK is a Python library that provides a platform for unguided and automated quantitative analysis of cytoskeletal images at high accuracy, which support both 2D and native 3D data structure. ILEE means Implicit Laplacian of Enhanced Edge, which is a local thresholding algorithm that detects and uses the edge of cytoskeleton filaments as the local baseline for fine thresholding. Our robust algorithm captures and computes native pixel/voxel value, first-order derivative (gradient), and second-order derivative (Laplacian) of the image to automatedly generate a threshold image for cytoskeleton segmentation. 

![scheme.PNG](https://github.com/phylars/image_host/blob/main/ILEE/scheme.PNG?raw=true)

The binary image of cytoskeleton fraction, together with the raw image is processed to compute more than 10 cytoskeleton indices, which are the matrices that measure the cytoskeletal features, of 5 classes: density, bundling, severing, branching, and directionality. The generated data can be used to determine and discover critical cytoskeletal feature of your own biological samples. 

In general, ILEE_CSK provides following state-of-art features to make outstanding performance:

* Support of native 3D mode, avoiding image projection that losses a considerable amount of information perpendicular to z-axis. 

* Unguided algorithm, avoiding human-bias related to subjective judgment in manual global thresholding, the currently most popular approach.

* High accuracy, the top among the kind. Please read our paper.

* Most complete cytoskeleton index list, with many novel indices developed by ourselves. 

# Presentation

Visualization of results processed by ILEE.

![demo1.PNG](https://github.com/phylars/image_host/blob/main/ILEE/demo1.PNG?raw=true)

Comparison of ILEE and other classic image segmentation approaches.

![demo2.PNG](https://github.com/phylars/image_host/blob/main/ILEE/demo2.PNG?raw=true)

![demo3.PNG](https://github.com/phylars/image_host/blob/main/ILEE/demo3.PNG?raw=true)
