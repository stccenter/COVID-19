<!--
 * @Author: your name
 * @Date: 2020-09-24 15:17:10
 * @LastEditTime: 2020-10-08 16:18:58
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: \Environmental Data - COVID-19\readme.md
-->
## What are the software requirements?

1. Python IDE
2. Anaconda 3

## How to set up the python environment and install packages?

1. Create a new python project with environment as conda.
![Conda Project](https://github.com/stccenter/COVID-19/blob/master/analysis/CA%20-%20Air%20Pollution/conda.png)
2. Download all python files (*.py) and requirements.txt file from repository and place it in python project.
3. Execute below comment to install required python packages.<br/>
	pip install -r requirements.txt

## Where to download the data?

1. For temperature, humidity, environmental condition – https://disc.gsfc.nasa.gov/datasets/M2T1NXSLV_5.12.4/summary?keywords=MERRA2_400.tavg1_2d_slv_Nx
2. For precipitation rate - https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGHHE_06/summary?keywords=IMERG
3. For nighttime light radiance, human activities, community distributions, human-gathering levels - https://ladsweb.modaps.eosdis.nasa.gov/
4. For concentration of air pollutants - https://disc.gsfc.nasa.gov/datasets/OMNO2d_003/summary?keywords=omi
5. For air quality index, air pollution concentration - http://data.cma.cn/

## How to get the results?**

1. Run script hourly2daily_humidity.py to convert hourly humidity data to daily. Change below variables accordingly.

        infolder – root directory for input files

        output_dir – root directory for output files

2. Run script hourly2daily_IMERG.py to convert hourly IMERG data to daily. Change below variables accordingly.

        infolder – root directory for input files

        output_dir – root directory for output files

        data_range – study period

3. Run script hourly2daily_temperature.py to convert temperature data to daily. Change below variables accordingly.

        infolder – root directory for input files

        output_dir – root directory for output files

## Results:
Find mean value statistics corresponding with regions below
1. Nighttime light radiance - https://github.com/stccenter/COVID-19/tree/master/analysis/nightlight
2. Air Quality - https://github.com/stccenter/COVID-19/tree/master/analysis/nightlight 
3. Precipitation - https://github.com/stccenter/COVID-19-Data/tree/master/Environmental%20factors 
4. Temperature - https://github.com/stccenter/COVID-19-Data/tree/master/Environmental%20factors 
5. Humidity - https://github.com/stccenter/COVID-19-Data/tree/master/Environmental%20factors 
6. NO2 TVCD - https://github.com/stccenter/COVID-19-Data/tree/master/Environmental%20factors 
 
Find global Distribution maps below
1. Nighttime light radiance - https://covid19datagmu.s3.us-east-2.amazonaws.com/NightTImeLight/NightTimeLight.zip
2. Precipitation<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/precipitation/daily/daily_precipitation_JAN_2020.zip<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/precipitation/daily/daily_precipitation_FEB_2020.zip<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/precipitation/daily/daily_precipitation_MAR_2020.zip<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/precipitation/daily/daily_precipitation_APR_2020.zip
3. Temperature<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/temperature_humidity/daily/Temperature/daily_MEAN_TEMP_JAN_2020.zip<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/temperature_humidity/daily/Temperature/daily_MEAN_TEMP_FEB_2020.zip<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/temperature_humidity/daily/Temperature/daily_MEAN_TEMP_MAR_2020.zip
4. Humidity<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/temperature_humidity/daily/+humidity/daily_MEAN_JAN_2020.zip<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/temperature_humidity/daily/+humidity/daily_MEAN_FEB_2020.zip<br/>
	https://covid19datagmu.s3.us-east-2.amazonaws.com/temperature_humidity/daily/+humidity/daily_MEAN_MAR_2020.zip

## Tutorial Video
[<img src="https://github.com/stccenter/COVID-19/blob/master/analysis/Environmental%20Data%20-%20COVID-19/youtube%20screen%20shot.png" width="60%">](https://youtu.be/Y_a-hs5n0oo)