**What are the software requirements?**

1. Python IDE
2. Anaconda 3
3. Microsoft Excel

**How to set up the python environment and install packages?**
1. Create a new python project with environment as conda
![Conda Project](https://github.com/stccenter/COVID-19/blob/master/analysis/CA%20-%20Air%20Pollution/conda.png)
2. Download python file (OMI_static_ca.py) and requirements.txt file from repository and place it in python project.
3. Execute below comment to install required python packages.<br/>
    pip install -r requirements.txt

**Where to download the data?**

Ground-based observations of air pollutants: https://www.epa.gov/outdoor-air-quality-data/download-daily-data

Satellite NO2 observations: [https://disc.gsfc.nasa.gov/datasets/OMNO2d_003/summary](https://disc.gsfc.nasa.gov/datasets/OMNO2d_003/summary)

National highways: https://catalog.data.gov/dataset/tiger-line-shapefile-2016-nation-u-s-primary-roads-national-shapefile

**How to get the results?**

Run the script OMI_statitic_ca.py to calculate periodical (pre, peri, and post) means of 2020 and 2015-2019, and their differences. Change the below variables accordingly.

1. period – either pre, peri, or post
2. infolder – root directory for input files
3. output_dir – root directory for output files
4. period1 – former period
5. period2 – latter period
6. path – root directory for input files

Results:

1. california_counties_covid_env_data.xlsx contains daily average concentration for each pollutant. Seven-day moving average and standard error are calculated in using excel.
