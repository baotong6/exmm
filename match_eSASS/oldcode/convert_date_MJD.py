from datetime import datetime
from astropy.time import Time

def convert_to_mjd(time_str):
    # 将时间字符串转换为datetime对象
    dt = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S.%f')
    
    # 使用astropy将datetime对象转换为MJD
    mjd = Time(dt).mjd
    
    return mjd