#!python
import psycopg2
from scipy.spatial import distance
import pandas as pd
import datetime
import argparse

ildatacols=['OmaTunniste','OmaIta','OmaPohjoinen','Kunta','siteid','aika','vuosi','kk','paiva','longitude','latitude','t_mean',\
            't_max','t_min','rainfall','radiation','hpa','lamposumma_v','rainfall_v']
            
def compare_distance(x):
    return x[0]

#Find closest point, euclidian distance
def find_closest_point(x,y,ls):
    ls1 = [(distance.euclidean([x,y],[x1,y1]),(x1,y1)) for (x1,y1,z) in ls]
    ls1 = sorted(ls1,key=compare_distance)
    return ls1[0]

#column 'aika' = str(year)+str(month)+str(day)
def date_to_aika(date):
    year=str(date.year)
    month=str(date.month)
    day=str(date.day)
    if date.month < 10:
        month=str(0)+month
    if date.day < 10:
        day=str(0)+day
    date=year+month+day
    return date

def add_aika_vuosi_kk_paiva(ls):
    ls1 = [(x,y,date_to_aika(date),date.year,date.month,date.day,x,y,temp_avg,temp_max,temp_min,prec,global_rad,vapour_press) for\
               (x,y,date,x,y,temp_avg,temp_max,temp_min,prec,global_rad,vapour_press)  in ls]
    return ls1

def add_fill_missing_column(df,col_num,name,initial_value):
    df.insert(col_num,name,initial_value)
    return df

def calculate_annual_rainfall(df):
    annual_rainfall = 0.0
    for i in range(len(df.index)):
        day = df.loc[i,'paiva']
        month = df.loc[i,'kk']
        if day == 1 and month == 1:
            annual_rainfall=0.0
        rainfall = df.loc[i,'rainfall']
        annual_rainfall=annual_rainfall+rainfall
        df.loc[i,'rainfall_v'] = annual_rainfall
    return df
        
def calculate_day_degree(df):
    day_degree=0.0
    dd_limit=5.0
    for i in range(len(df.index)):
        day = df.loc[i,'paiva']
        month = df.loc[i,'kk']
        if day == 1 and month == 1:
            day_degree=0.0
        t_mean=df.loc[i,'t_mean']
        dd_tmp=t_mean-dd_limit
        day_degree = day_degree + max(0,dd_tmp)
        df.loc[i,'lamposumma_v']=day_degree
    return df
    
def fetch_ildata(host,database,user,password,query): 
    conn = psycopg2.connect(host=host,database=database,user=user,password=password) 
    cur = conn.cursor()
    #cur.execute('SELECT version()')                                                                      
    #version=cur.fetchone()
    #print(version)
    cur.execute(query)
    ls = cur.fetchall()
    conn.close()
    return ls

def location_to_ildata(x,y,user,password):
    n=x+1000
    s=x-1000
    #day before yesterday
    dayday=(datetime.date.today()-datetime.timedelta(2)).strftime('%Y-%m-%d')
    #SQL syntax requires the use of single quotes (i.e. the "'"-character)
    daydaydate="'"+dayday+"'"+"::date"
    #Fetch first a couple of days of data
    ls1 = fetch_ildata(host="lukedb1.ns.luke.fi",database='weather',user=user,password=password,
                       query="SELECT x,y,date FROM grid_day WHERE date >="+daydaydate)
    t = find_closest_point(x,y,ls1)
    x1=t[1][0]
    y1=t[1][1]
    #With the closest point select all data  
    ls = fetch_ildata(host="lukedb1.ns.luke.fi",database='weather',user=user,password=password,
                      query="SELECT x,y,date,x,y,temp_avg,temp_max,temp_min,prec,global_rad,vapour_press FROM grid_day WHERE x="+str(x1)+" AND y="+str(y1)+" AND date > '1961-01-01'::date")
    return ls

def write_ls_to_excel(file_name,ls,column_ls):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')
    df = pd.DataFrame(ls)
    df.columns=column_ls
    df.to_excel(writer,sheet_name='ILData')
    writer.save()

def write_df_to_excel(file_name,sheet_name,df):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')
    df.to_excel(writer,sheet_name=sheet_name)
    writer.save()

def write_df_to_csv(file_path,df,sep=';'):
    df.to_csv(file_path,sep)

def read_and_write_weather_data(x,y,loc,file_path,user,passwd):
    ls = location_to_ildata(x,y,user,passwd)
    ls1 = add_aika_vuosi_kk_paiva(ls)
    df = pd.DataFrame(ls1)
    df = add_fill_missing_column(df,0,'OmaTunniste',1)
    df = add_fill_missing_column(df,3,'Kunta',loc)
    df = add_fill_missing_column(df,4,'siteid',9999)
    df = add_fill_missing_column(df,len(df.columns),'lamposumma_v',0)
    df = add_fill_missing_column(df,len(df.columns),'rainfall_v',0)
    df.columns=ildatacols
    df['OmaIta']=x
    df['OmaPohjoinen']=y
    df=calculate_day_degree(df)
    df=calculate_annual_rainfall(df)
    df.to_csv(file_path,sep=';',index=False)
    excel_file = file_path.replace('.csv','.xlsx')
    write_df_to_excel(excel_file,'ILData',df)
    return df

#Vihti research site: (x=3356593 AND y=6703104)
def test(x=3356593,y=6703104):
    ls = location_to_ildata(x,y)
    ls1 = add_aika_vuosi_kk_paiva(ls)
    df = pd.DataFrame(ls1)
    df = add_fill_missing_column(df,0,'OmaTunniste',1)
    df = add_fill_missing_column(df,3,'Kunta','Vihti')
    df = add_fill_missing_column(df,4,'siteid',9999)
    df = add_fill_missing_column(df,len(df.columns),'lamposumma_v',0)
    df = add_fill_missing_column(df,len(df.columns),'rainfall_v',0)
    df.columns=ildatacols
    df['OmaIta']=x
    df['OmaPohjoinen']=y
    df=calculate_day_degree(df)
    df=calculate_annual_rainfall(df)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-x',type=int,dest='x',help='East coordinate')
    parser.add_argument('-y',type=int,dest='y',help='North coordinate')
    parser.add_argument('-l',type=str,dest='l',default='xxxx',help='Location name')
    parser.add_argument('-f',type=str,dest='f',help='Output file name (type CSV)')
    parser.add_argument('-u','--user',type=str,dest='u',help='User name (Luke network) for weather database')
    parser.add_argument('-p','--passwd',type=str,dest='p',help='Password for weather database')
    args=parser.parse_args()
    #print(args.x,args.y,args.f,args.u,args.p)
    read_and_write_weather_data(args.x,args.y,args.l,args.f,args.u,args.p)
