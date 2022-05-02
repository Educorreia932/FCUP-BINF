# Assignment 4

## How to run

```sh
python similarity.py
```

## Exercises 

### Exercise 1

```py
<Protein>: <Score>
YP_009724390.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]: 6722
QHO60594.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]: 6722
QHR63300.2 spike glycoprotein [Bat coronavirus RaTG13]: 6523
QVT76606.1 spike protein [Pangolin coronavirus]: 6237
AVP78042.1 spike protein [Bat SARS-like coronavirus]: 5281
AAP13441.1 S protein [SARS coronavirus Urbani]: 5145
AAP41037.1 spike glycoprotein [SARS coronavirus Tor2]: 5142
AAU04646.1 spike glycoprotein [Civet SARS CoV 007/2004]: 5141
YP_009047204.1 spike protein [Middle East respiratory syndrome-related coronavirus]: 1190
AMK59677.1 S protein [Human coronavirus OC43]: 976
```

### Exercise 3

#### Scores matrix

```py
                  YP_009724390.1      QHO60594.1      AAP13441.1      QHR63300.2      AAU04646.1      AMK59677.1      AAP41037.1  YP_009047204.1      QVT76606.1      AVP78042.1
YP_009724390.1	            6722            6722            5145            6523            5141             976            5142            1190            6237            5281
    QHO60594.1	            6722            6722            5145            6523            5141             976            5142            1190            6237            5281
    AAP13441.1	            5145            5145            6632            5190            6516             934            6629            1113            5159            4925
    QHR63300.2	            6523            6523            5190            6704            5183             948            5187            1117            6283            5281
    AAU04646.1	            5141            5141            6516            5183            6624             937            6513            1106            5152            4897
    AMK59677.1	             976             976             934             948             937            7248             933            1150             945             788
    AAP41037.1	            5142            5142            6629            5187            6513             933            6632            1115            5156            4922
YP_009047204.1	            1190            1190            1113            1117            1106            1150            1115            7124            1114            1024
    QVT76606.1	            6237            6237            5159            6283            5152             945            5156            1114            6696            5308
    AVP78042.1	            5281            5281            4925            5281            4897             788            4922            1024            5308            6557
```

#### Mismatches matrix

```py
                  YP_009724390.1      QHO60594.1      AAP13441.1      QHR63300.2      AAU04646.1      AMK59677.1      AAP41037.1  YP_009047204.1      QVT76606.1      AVP78042.1
YP_009724390.1	               0               0             299              33             302             957             300             938              97             252
    QHO60594.1	               0               0             299              33             302             957             300             938              97             252
    AAP13441.1	             299             299               0             292              20             951               1             923             295             298
    QHR63300.2	              33              33             292               0             296             962             293             943              88             254
    AAU04646.1	             302             302              20             296               0             951              21             924             297             303
    AMK59677.1	             957             957             951             962             951               0             951             966             958             963
    AAP41037.1	             300             300               1             293              21             951               0             923             296             299
YP_009047204.1	             938             938             923             943             924             966             923               0             943             930
    QVT76606.1	              97              97             295              88             297             958             296             943               0             254
    AVP78042.1	             252             252             298             254             303             963             299             930             254               0
```

## Exercise 5.1

```py
<Protein>: <Score>
YP_009724390.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]: 6722
QHO60594.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]: 6722
QHR63300.2 spike glycoprotein [Bat coronavirus RaTG13]: 6523
QVT76606.1 spike protein [Pangolin coronavirus]: 6237
AVP78042.1 spike protein [Bat SARS-like coronavirus]: 5281
AAP13441.1 S protein [SARS coronavirus Urbani]: 5145
AAP41037.1 spike glycoprotein [SARS coronavirus Tor2]: 5142
AAU04646.1 spike glycoprotein [Civet SARS CoV 007/2004]: 5141
YP_009047204.1 spike protein [Middle East respiratory syndrome-related coronavirus]: 1393
AMK59677.1 S protein [Human coronavirus OC43]: 1164
```

## Exercise 5.3

#### Scores matrix

```py
                  YP_009724390.1      QHO60594.1      AAP13441.1      QHR63300.2      AAU04646.1      AMK59677.1      AAP41037.1  YP_009047204.1      QVT76606.1      AVP78042.1
YP_009724390.1	            6722            6722            5145            6523            5141            1164            5142            1393            6237            5281
    QHO60594.1	            6722            6722            5145            6523            5141            1164            5142            1393            6237            5281
    AAP13441.1	            5145            5145            6632            5190            6516            1118            6629            1350            5159            4928
    QHR63300.2	            6523            6523            5190            6704            5183            1158            5187            1327            6283            5281
    AAU04646.1	            5141            5141            6516            5183            6624            1131            6513            1356            5152            4900
    AMK59677.1	            1164            1164            1118            1158            1131            7248            1118            1317            1132            1158
    AAP41037.1	            5142            5142            6629            5187            6513            1118            6632            1352            5156            4925
YP_009047204.1	            1393            1393            1350            1327            1356            1317            1352            7124            1323            1275
    QVT76606.1	            6237            6237            5159            6283            5152            1132            5156            1323            6696            5308
    AVP78042.1	            5281            5281            4928            5281            4900            1158            4925            1275            5308            6557
```

#### Mismatches matrix

```py
                  YP_009724390.1      QHO60594.1      AAP13441.1      QHR63300.2      AAU04646.1      AMK59677.1      AAP41037.1  YP_009047204.1      QVT76606.1      AVP78042.1
YP_009724390.1	               0               0             299              33             302             458             300             674              97             252
    QHO60594.1	               0               0             299              33             302             458             300             674              97             252
    AAP13441.1	             299             299               0             292              20             469               1             675             295             296
    QHR63300.2	              33              33             292               0             296             459             293             676              88             254
    AAU04646.1	             302             302              20             296               0             467              21             674             297             301
    AMK59677.1	             458             458             469             459             467               0             469             432             459             455
    AAP41037.1	             300             300               1             293              21             469               0             675             296             297
YP_009047204.1	             674             674             675             676             674             432             675               0             679             671
    QVT76606.1	              97              97             295              88             297             459             296             679               0             254
    AVP78042.1	             252             252             296             254             301             455             297             671             254               0
```    