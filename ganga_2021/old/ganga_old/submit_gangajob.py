import os

yearlist_run1 = ['2011','2012']
typelist_run1 = [
'12197022',
# '11198020',
'12197420',
# '11196420',
# '11196410',
# '12197220',
# '11196210',
# '11196200',
# '11196620',
]
yearlist_run2 = ['2016']
typelist_run2 = [
# # '11198000',
# # '11198005',
# # '11198400',
# # '11198401',
# '11198410',
# # # "12197400",
# "12197401",
# # "11196000",
# # "13198040",
# # "13198200",
# # "13198400",
# # "13198600",
# "11198030",
# "12197009",
]

for year in yearlist_run1:
    for type in typelist_run1:
        os.system("ganga bddk/gangajob_BDDK.py --TYPE {0} --YEAR {1} --TESTING no".format(type,year))
#
# for year in yearlist_run2:
#     for type in typelist_run2:
#         os.system("ganga bddk/gangajob_BDDK.py --TYPE {0} --YEAR {1} --TESTING no".format(type,year))

# os.system("ganga bddk/gangajob_BDDK.py --TYPE Data --YEAR 2016 --TESTING yes")
# os.system("ganga bddk/gangajob_BDDK.py --TYPE Data --YEAR 2017 --TESTING no")
# os.system("ganga bddk/gangajob_BDDK.py --TYPE Data --YEAR 2018 --TESTING no")

os.system("ganga -i")
