**What are the software requirements?**

1. Python IDE
2. Anaconda 3
3. Microsoft Excel

**How to set up the python environment and install packages?**

1. Create a new python project with environment as conda.
![Conda Project](https://github.com/stccenter/COVID-19/blob/master/analysis/China%20-%20NTL%20AQI/conda.png)
2. Download python files (*.py) and requirements.txt file from repository and place it in under python project.
3. Execute below comments in terminal to install required python packages:<br/>
            1. conda install -c conda-forge  hdf4<br/>
            2. conda install -c conda-forge  pandas<br/>
            3. conda install -c conda-forge  basemap<br/>
            4. conda install -c conda-forge pyhdf<br/>
            5. pip install -r requirements.txt<br/>

**Where to download the data?**

Nighttime light: [https://ladsweb.modaps.eosdis.nasa.gov/search/order/3/VNP46A1--5000/2020-02-01..2020-04-03/DB/World](https://ladsweb.modaps.eosdis.nasa.gov/search/order/3/VNP46A1--5000/2020-02-01..2020-04-03/DB/World)

Air Quality Index: [http://www.cnemc.cn/](http://www.cnemc.cn/)

**How to get the results?**

Run the below scripts to generate nighttime light and AQI data. Install required packages for the scripts.

1. VIIRS-VNP46A1.py

This script is used to calculate the daily nighttime light radiances for the entire target region. Change following variables accordingly.

    1. data_infolder– root directory for input files
    2. output_dir – root directory for output files
    3. alldate – study period

2. dnb-monthly-VNP.py

This script is used to calculate the monthly mean of nighttime light over the target region. Change the following variables accordingly.

    1. input_dir – root directory for input files
    2. output_dir – root directory for output files
    3. alldate – study period

3. dnb-difference-VNP.py

This script is used to calculate the difference of nighttime light between different months. Change the following variables accordingly.

    1. input_dir – root directory for input files
    2. output_dir – root directory for output files
    3. region – study area
    4. pre_nc_filename – monthly mean of nighttime light in China during former month
    5. post_nc_filename – monthly mean of nighttime light in China during latter month

**Results**

1. NTL-2019.xlsx contains the monthly mean of nighttime light in all the provinces of mainland China in 2019.
2. Nighttime_light_data.xlsx contains monthly mean nighttime light of provinces and entire country of China.
3. AQI_data.xlsx contains daily and monthly mean of AQI in provinces and entire country of China

**Tutorial Video**

[<img src="https://github.com/stccenter/COVID-19/blob/master/analysis/China%20-%20NTL%20AQI/Youtube%20screenshot.png" width="60%">](https://www.youtube.com/watch?v=pWMRpG_2-Fs&t=7s)
