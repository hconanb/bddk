import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *

def saveplot(canvas, name):
    now = datetime.datetime.now()
    tag = name.split("_")[-1]
    if not os.path.exists(f'plots/{now.month}_{now.day}_pdfs/{tag}/'):
        os.makedirs(f'plots/{now.month}_{now.day}_pdfs/{tag}/')
    if not os.path.exists(f'plots/{now.month}_{now.day}_pngs/{tag}/'):
        os.makedirs(f'plots/{now.month}_{now.day}_pngs/{tag}/')
    canvas.SaveAs(f"plots/{now.month}_{now.day}_pdfs/{tag}/{name}.pdf")
    print(f"Saved: plots/{now.month}_{now.day}_pdfs/{tag}/{name}.pdf")
    canvas.SaveAs(f"plots/{now.month}_{now.day}_pngs/{tag}/{name}.png")
    print(f"Saved: plots/{now.month}_{now.day}_pdfs/{tag}/{name}.png")



def savesplitplot(list, split1, split2, atag):
    canvas = ROOT.TCanvas("c1","c1")
    canvas.Divide(split1, split2)
    name = list[0].GetTitle()
    for index in range(len(list)):
        canvas.cd(index+1)
        list[index].Draw()
    now = datetime.datetime.now()
    tag = name.split("_")[-1]
    if not os.path.exists(f'plots/{now.month}_{now.day}_pdfs/{tag}/'):
        os.makedirs(f'plots/{now.month}_{now.day}_pdfs/{tag}/')
    if not os.path.exists(f'plots/{now.month}_{now.day}_pngs/{tag}/'):
        os.makedirs(f'plots/{now.month}_{now.day}_pngs/{tag}/')
    canvas.SaveAs(f"plots/{now.month}_{now.day}_pdfs/{tag}/{atag}.pdf")
    print(f"Saved: plots/{now.month}_{now.day}_pdfs/{tag}/{atag}.pdf")
    canvas.SaveAs(f"plots/{now.month}_{now.day}_pngs/{tag}/{atag}.png")
    print(f"Saved: plots/{now.month}_{now.day}_pdfs/{tag}/{atag}.png")
