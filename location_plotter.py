import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()



xlsfile=pd.ExcelFile('locations.xlsx')
df=xlsfile.parse()
df.Latitude=(df.Latitude**3)*(7.96749E-05)+(df.Latitude**2)*(0.000361054)+(df.Latitude**1)*(0.678924168)-24.18847156



BBox = (-180,   180,
         -90, 90)

ruh_m = plt.imread('cold_map.png')

fig, ax = plt.subplots(figsize = (8,7))

#ax.scatter(df.Longitude, df.Latitude, zorder=1, alpha= 1, c='r', s=(2**(df.Heat/500))*2)

#ax2 = sns.scatterplot(x="Longitude", y="Latitude", hue="Heat", size="Heat", sizes=(100, 200), hue_norm=(1000, 2000), data=df)
cmap = sns.cubehelix_palette(dark=0.8, light=0.3, as_cmap=True)
ax2 = sns.scatterplot(x="Longitude", y="Latitude", hue="Heat", hue_norm=(1000, 2500),data=df, palette=cmap)

ax.set_title('Coldest Cities in Canada, United States and Europe')
ax.set_xlim(BBox[0],BBox[1])
ax.set_ylim(BBox[2],BBox[3])

ax.imshow(ruh_m, zorder=0, extent = BBox, aspect= 'equal')

for i in df.index:
    #ax.annotate(str(df.Location[i]), (df.Longitude[i], df.Latitude[i]))
    ax.annotate(str(i+1), (df.Longitude[i], df.Latitude[i]))