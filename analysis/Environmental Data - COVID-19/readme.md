**What are the software requirements?**

1. Python IDE

**Where to download the data?**

1. For temperature, humidity, environmental condition –

https://disc.gsfc.nasa.gov/ datasets/M2T1NXSLV_5.12.4/ summary?keywords=MERRA2_ 400.tavg1_2d_slv_Nx

1. For precipitation rate - https://disc.gsfc.nasa.gov/ datasets/GPM_3IMERGHHE_06/ summary?keywords=IMERG
2. Fornighttime light radiance, human activities, community distributions, human-gathering levels - https://ladsweb.modaps.eosdis. nasa.gov/
3. For concentration of air pollutants - https://disc.gsfc.nasa.gov/ datasets/OMNO2d_003/ summary?keywords=omi
4. For air quality index, air pollution concentration - http://data.cma.cn/

**How to get the results?**

1. Run script hourly2daily_humidity.py to convert hourly humidity data to daily. Change below variables accordingly.

infolder – root directory for input files

output_dir – root directory for output files

1. Run script hourly2daily_IMERG.py to convert hourly IMERG data to daily. Change below variables accordingly.

infolder – root directory for input files

output_dir – root directory for output files

data_range – study period

1. Run script hourly2daily_temperature.py to convert temperature data to daily. Change below variables accordingly.

infolder – root directory for input files

output_dir – root directory for output files

**Results** :

|
 | **Derived Values** | **Downloading Path** |
| --- | --- | --- |
| Mean value statistics corresponding with regions | Nighttime light radiance | https://github.com/stccenter/COVID- 19/tree/master/analysis/nightlight |
| Air Quality | https://github.com/stccenter/COVID- 19/tree/master/analysis/nightlight |
| Precipitation | https://github.com/stccenter/COVID- 19-Data/tree/master/Environmental% 20factors |
| Temperature | https://github.com/stccenter/COVID- 19-Data/tree/master/Environmental% 20factors |
| Humidity | https://github.com/stccenter/COVID- 19-Data/tree/master/Environmental% 20factors |
| NO2 TVCD | https://github.com/stccenter/COVID- 19-Data/tree/master/Environmental% 20factors |
| Global Distribution maps | Nighttime light radiance | https://covid19datagmu.s3.us-east-2. amazonaws.com/NightTImeLight/ NightTimeLight.zip |
| Precipitation | https://covid19datagmu.s3.us-east-2. amazonaws.com/precipitation/daily/ daily_precipitation_JAN_2020.zip https://covid19datagmu.s3.us-east-2. amazonaws.com/precipitation/daily/ daily_precipitation_FEB_2020.zip https://covid19datagmu.s3.us-east-2. amazonaws.com/precipitation/daily/ daily_precipitation_MAR_2020.zip https://covid19datagmu.s3.us-east-2. amazonaws.com/precipitation/daily/ daily_precipitation_APR_2020.zip |
| Temperature | https://covid19datagmu.s3.us-east-2. amazonaws.com/temperature_ humidity/daily/Temperature/daily_ MEAN_TEMP_JAN_2020.zip https://covid19datagmu.s3.us-east-2. amazonaws.com/temperature_ humidity/daily/Temperature/daily_ MEAN_TEMP_FEB_2020.zip https://covid19datagmu.s3.us-east-2. amazonaws.com/temperature_ humidity/daily/Temperature/daily_ MEAN_TEMP_MAR_2020.zip |
| Humidity | https://covid19datagmu.s3.us-east-2. amazonaws.com/temperature_ humidity/daily/+humidity/daily_ MEAN_JAN_2020.zip https://covid19datagmu.s3.us-east-2. amazonaws.com/temperature_ humidity/daily/+humidity/daily_ MEAN_FEB_2020.zip https://covid19datagmu.s3.us-east-2. amazonaws.com/temperature_ humidity/daily/+humidity/daily_ MEAN_MAR_2020.zip |
