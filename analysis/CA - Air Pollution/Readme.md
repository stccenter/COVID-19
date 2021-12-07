## **Research Background**

The mitigation policies or lockdown measures imposed by local and national gov-ernments to control the spread of COVID-19 have created drastic changes in people's livelihood, economy, and environment. We investigated the impacts of the COVID-19 crisis on the atmospheric environment by comprehensively understanding the spatio-temporal patterns of air pollutants in California, USA. The investigation results provided knowledge about humans and environmental re-sponses to the pandemic.

**Prerequisite**

1. Python IDE (Visual Studio Code, Pycharm, or anything of your choice)
2. Python 2.7 or above
3. Anaconda 3

**Steps**

* Step 1: Download the project related materials using the link:

* Step 2: The downloaded folder has below folder:
    1. Air Pollutants Data - This folder has satellite NO2 observations data that is splitted into three periods - pre, peri, and post
    2. Air Quality Analytical Tool
   
* Step 3: Go to the terminal. 

* Step 4: Create a new conda environment.

**For Windows and Mac**
```
conda create -n env-analysis-no2
```

* Step 5: Activate the conda environment.

**For Windows and Mac**
```
conda activate env-analysis-no2
```

* Install python packages


1. Download python file (OMI_static_ca.py) and requirements.txt file from repository and place it in python project.
2. Execute below comment to install required python packages.<br/>
    pip install -r requirements.txt

**How to get the results?**

Run the script OMI_statitic_ca.py to calculate periodical (pre, peri, and post) means of 2020 and 2015-2019, and their differences. Change the below variables accordingly.

1. period – either pre, peri, or post
2. infolder – root directory for input files
3. output_dir – root directory for output files
4. period1 – former period
5. period2 – latter period
6. path – root directory for input files

**Results**

1. california_counties_covid_env_data.xlsx contains daily average concentration for each pollutant. Seven-day moving average and standard error are calculated in using excel.

* Useful links

1. Ground-based observations of air pollutants: https://www.epa.gov/outdoor-air-quality-data/download-daily-data

2. Satellite NO2 observations: [https://disc.gsfc.nasa.gov/datasets/OMNO2d_003/summary](https://disc.gsfc.nasa.gov/datasets/OMNO2d_003/summary)

3. National highways: https://catalog.data.gov/dataset/tiger-line-shapefile-2016-nation-u-s-primary-roads-national-shapefile


**Author**

Qian Liu<br>
Email: qliu6@gmu.edu


**Tutorial Video**

[<img src="https://github.com/stccenter/COVID-19/blob/master/analysis/CA%20-%20Air%20Pollution/Screenshots/Screenshot%20for%20video.png" width="60%">](https://www.youtube.com/watch?v=hwQF3_ZJSJY)

