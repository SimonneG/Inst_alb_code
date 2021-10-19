# Inst_alb_code
Scripts to select cer-pol colocated couples
selection_albedos.py takes two arguments : -f POLDER filename (which incudes the path to the file) and -c CERES date which helps to choose the right CERES files the script calls. 
For each POLDER orbit, the timespan of the orbit is identified (eg from 9:45 am to 11:05am) and the corresponding 1h CERES orbits are called in a loop. In this example, three CERES orbits are called : from 9:00 am to 10:00 am, from 10:00 am to 11:00 am and from 11:00 am to 12:00 am. 
Within the loop, all the polder superpixels within a ceres pixel are averaged (up to 36 values) to a single value, then the albedo value of the ceres pixel is printed along the averaged polder value, the coordinates of the ceres pixel and averaged coordinates of the polder superpixels. Other parameters are also printed (hour, solar zenith angle at the time of measurement). 

