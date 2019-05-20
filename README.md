Intima-media Segmentation ReadMe

The repository contains the Matlab implementations of the parallel boundary detection algorithm proposed in:  

Zhou, Y.; Cheng, X.; Xu, X. & Song, E.  
Dynamic programming in parallel boundary detection with application to ultrasound intima-media segmentation  
Medical image analysis, Elsevier, 2013, 17, 892-906  
(https://www.sciencedirect.com/science/article/abs/pii/S1361841513000832)


--------------------------------- Usage -------------------------------------

Simply running "LDLD.m" will execute the algorithm on the attached synthetic image.  

The "dld_dp" folder contains the c++ files for the dynamic programming algorithm (compiled into dld_dp.mexw64 and dld_dp.mexw32).  

"LDLD.m" - main demo file.

"imt_snake_model.m" - implements the final curve refinement step.  

"scale_multi.m" - calculates the edge map.  


-------------------------------- Contact ----------------------------------

If you have any questions, please contact:

Yuan Zhou  

zhouyuanzxcv@gmail.com

