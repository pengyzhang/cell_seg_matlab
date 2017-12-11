This is the training dataset for the digital pathology segmentation challenge. In this 
dataset there are 15 images. For each image, the dataset contains three files: 

imageXX.png -- the original image 
imageXX_mask.txt -- the manual segmentation results as labeled masks in ASCII 
					text file format. 
imageXX_poly.png -- the rendering of the segmentation results on the image

imageXX_poly.png is provided for visualization of the segmentation results only. 
You should use imageXX.png and imageXX_mask.txt files to train your segmentation 
algorithms. 

The 15 images are extracted as rectangular tiles from whole slide tissue images at the 
highest resolution. The image capture magnifications (objective powers) at which the 
source whole slide images were obtained are as follows: 

Image Name		Magnification
------------	-------------
image01.png 		20x
image02.png			20x
image03.png 		20x
image04.png 		20x
image05.png			20x
image06.png			20x
image07.png			20x
image08.png			40x
image09.png			40x
image10.png			40x
image11.png			40x
image12.png			40x
image13.png			40x
image14.png			40x
image15.png			40x

The format of the ASCII labeled mask files is as follows. Each labeled mask file is an 
array of integer values. The array has the same resolution as the image -- i.e., the 
same width and height. Each array element stores a value between 0 and N, where N is the 
number of segmented nuclei. A value of 0 in an array element means the corresponding 
pixel in the image tile is not part of a nucleus. A non-zero value means the corresponding 
pixel in the image tile is part of a nucleus. All the pixels that belong to the same 
nucleus are labeled with the same value. The mask array is organized as following in the 
mask file: 

width height    
pixel_label_id  // pixel (0,0)
pixel_label_id  // pixel (1,0)
pixel_label_id  // pixel (2,0)
…
pixel_label_id  // pixel (width-1,height-1)

For example, assume two nuclei were segmented in an image of 5x4 pixels:

00000
01100
11200
02220

The mask file would have the following content:

5 4
0
0
0
0
0
0
1
1
0
0
1
1
2
0
0
0
2
2
2
0
