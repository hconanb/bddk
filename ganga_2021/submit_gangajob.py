import os
# '12197024',
# '12197422',
# '11198022',
# Old MC

typelist_run2 = [
'11198006',
# # '11198005',
# #### August 24th 2021
# #### January 12 2022
#
# # '11198400',
# # '11198410',
# # '11198401',
# # '13198040',
# # '13198400',
# # '13198200',
# # '13198600',
#
# #### August 25th 2021
# #### January 13th 2022
#
# # '12197400',
# # '12197401',
# # '12197410',
# # #### August 26th 2021
# # '12197023',
# # '11196414',
# # '11196413',
# # '11196019',
# # '12197008',
# # '11198007',
# #### August 30th 2021
# "12197423",
# "12197045",
# "11198023",
# #### September 16th 2021
]
pollist = ['Up','Down']
yearlist_run2 = ['2016','2017','2018']

for type in typelist_run2:
    for year in yearlist_run2:
        for pol in pollist:
            os.system("ganga gangajob_BDDK.py --TYPE {0} --YEAR {1} --POL {2} --TESTING no".format(type,year,pol))

# for year in yearlist_run2:
#     for pol in pollist:
#         os.system(f"ganga bddk/gangajob_BDDK.py --TYPE Data --YEAR {year} --POL {pol} --TESTING no")


os.system("ganga -i")
