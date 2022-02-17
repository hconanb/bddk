
import ROOT

def autocropXaxis( histlist ):
    """Remove zero bins (on every histogram in the list) on left and right of x-axis automatically, by setting same range for all."""

    if len(histlist) == 0:
        print ("Why do you ask me to crop no histograms?")
        return

    firstbin = 1
    nbins = histlist[0].GetNbinsX()
    lastbin = nbins

    while ( all(h.GetBinContent(firstbin) == 0 for h in histlist) ):
        firstbin += 1
        if firstbin == nbins:
            return

    #leave a single empty bin to make it clear things cut off
    if firstbin > 1:
        firstbin -= 1

    while ( all(h.GetBinContent(lastbin) == 0 for h in histlist) ):
        lastbin -= 1

    #leave a single empty bin to make it clear things cut off
    if lastbin < histlist[0].GetNbinsX():
        lastbin += 1

    for h in histlist:
        h.GetXaxis().SetRange( firstbin, lastbin )
def histomax( *hists ):
    try: ##works if you give it a list of histograms
        hlist = [ h for sub in hists for h in sub ]
        return max( [ h.GetBinContent(h.GetMaximumBin()) for h in hlist] )
    except:
        ## works if you give it individual hists as args
        return max( [ h.GetBinContent(h.GetMaximumBin()) for h in hists] )
def graphmax( gr, noError = False ):
    m = 0.0
    x = ROOT.Double()
    y = ROOT.Double()
    for i in range(gr.GetN()):
        gr.GetPoint( i, x, y )

        e = gr.GetErrorYhigh( i )

        if noError and float(y) > m:
            m = float(y)
        elif  float(y) + e > m:
            m = float(y)+e

    return m
def graphmin( gr, noError = False ):
    m = 1e12
    x = ROOT.Double()
    y = ROOT.Double()
    for i in range(gr.GetN()):
        gr.GetPoint( i, x, y )
        e = gr.GetErrorYlow( i )
        if noError and float(y) < m:
            m = float(y)
        elif (not noError) and float(y) - e < m:
            m = float(y)-e

    return m
class compPlot():
    """Helper class that plots histograms for comparison. Auto setting of axis range"""
    def __init__(self,name='compPlot'):
        self.pt = ROOT.TCanvas(name,name)


    def save(self, path):
        self.pt.SaveAs(path)


    def setref(self, h ):
        ## helps keep interface same as ratio plot
        self.hlist = [h]
        self.href = h

    def addhists(self, *h):
        try:
            self.hlist.extend(h)
        except:
            self.hlist = list(h)


    def reset(self):
        self.hlist = None

    def draw(self, opts='', runAutoCrop = False, axismult = 1.2):
        self.pt.cd()

        if runAutoCrop:
            autocropXaxis( self.hlist )

        maxbin = max( [ h.GetBinContent(h.GetMaximumBin()) for h in self.hlist] )

        self.hlist[0].Draw(opts)
        self.hlist[0].GetYaxis().SetRangeUser(0.0, axismult*maxbin)
        for h in self.hlist[1:]:
            h.Draw("SAME"+opts)
class splitPlot():
    """Helper class that creates a canvas with pads above and below with no margin in between"""
    def __init__(self, fr = 1./3., name='splitPlot'):
        basew = ROOT.gStyle.GetCanvasDefW()
        baseh = ROOT.gStyle.GetCanvasDefH()

        bmargin = ROOT.gStyle.GetPadBottomMargin()
        tmargin = ROOT.gStyle.GetPadTopMargin()

        ##frame height in relative coord
        frameh = 1. - bmargin - tmargin

        ##bottom frame height over top frame
        self.frameratio = fr

        ##frameh*frameratio is the new "relative" length of the canvas
        newscale = 1 + frameh*self.frameratio
        ##actual pixel size of new canvas
        newh = baseh*newscale

        ##new margins -- fractions of the whole canvas
        newbmargin = bmargin / newscale
        newtmargin = tmargin / newscale

        ##where the pads split
        ybreak = ( bmargin + frameh*self.frameratio )/newscale
        self.padratio = ybreak/(1-ybreak)

        ##how much of the pad does the margin take?
        actualtmargin = newtmargin/(1-ybreak)
        actualbmargin = newbmargin/ybreak

        ## the label ratio is not the padratio, but the bottom pad size relative to the default
        self.tlabelratio = baseh / ( newh*(1-ybreak) )
        self.blabelratio = baseh / ( newh*ybreak )

        self.c = ROOT.TCanvas(name, name,basew, int(newh) )
        self.c.cd()
        self.pb = ROOT.TPad("pb","pb",0,0,1,ybreak)
        self.pb.Draw()
        self.c.cd()
        self.pt = ROOT.TPad("pt","pt",0,ybreak,1,1)
        self.pt.Draw()

        self.pt.cd()
        self.pt.SetBottomMargin(0.02)
        self.pt.SetTopMargin(actualtmargin)
        self.pt.DrawFrame(0,0,1,1)


        self.pb.cd()
        self.pb.SetBottomMargin(actualbmargin)
        self.pb.SetTopMargin(0.03)
        self.pb.DrawFrame(0,0,1,1)


    def save(self, path):

        self.c.SaveAs(path)
class ratioPlot(splitPlot):
    """Helper class to create a canvas with histograms above and ratio plot below"""

    def setref(self, href):
        self.href = href.Clone("href")

    def addhists(self, *h):
        try:
            self.hlist.extend(h)
        except:
            self.hlist = list(h)


    def reset(self):
        self.href = None
        self.hlist = None

    def draw(self, topopts='', bottomopts='', runAutoCrop = False, axismult = 1.2):


        ## Notes on options:
        ## top and bottomopts are used for draw styles for histograms on top and bottom frames
        ## runAutoCrop will use the autocrop method from this utils file to crop empty bins on the xaxis
        ## axismult controls the top frame Y axis scale relative  to the maximum bin. if you have a legend, you can set this to
        ## ( 1 - topmargin - bottommargin ) /( legend_y1 - bottommargin - dy ), legend_y1 is where the legend starts in NDC coords and where dy is how far in NDC coordinates above the max bin the legend should start

        self.pt.cd()
        self.href.Draw(topopts)

        maxbin = max( [ h.GetBinContent(h.GetMaximumBin()) for h in self.hlist] + [self.href.GetBinContent(self.href.GetMaximumBin()),] )

        self.href.GetXaxis().SetLabelSize(0)
        self.href.GetXaxis().SetTitleSize(0)
        self.href.GetYaxis().SetRangeUser(1e-6, maxbin*axismult)
        self.href.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()*self.tlabelratio)
        self.href.GetYaxis().SetTitleSize(ROOT.gStyle.GetLabelSize()*self.tlabelratio)

        if runAutoCrop:
            autocropXaxis( [self.href,] + self.hlist )

        for h in self.hlist:
            h.Draw("SAME"+topopts)

        self.href.Draw("SAME"+topopts)

        self.pb.cd()

        self.ratiolist = [ h.Clone(h.GetName() + "clone") for h in self.hlist ]
        for h in self.ratiolist:
            #h.GetXaxis().SetRange( self.href.GetXaxis().GetFirst(), self.href.GetXaxis().GetLast() )
            h.Divide(self.href)

        if len(self.ratiolist) > 0:
            maxratio = max( [ h.GetBinContent(h.GetMaximumBin()) for h in self.ratiolist ] )
            minratio = min( [ h.GetBinContent(h.GetMinimumBin()) for h in self.ratiolist ] )
        else:
            maxratio = 2
            minratio = 0

        up = max( maxratio - 1., 1 - minratio )
        upratio = round( up,1 )
        if upratio < up:
            upratio += 0.1

        downratio = min(upratio,1.)

        self.bottomframe = self.href.Clone("pbframe")
        self.bottomframe.GetYaxis().SetRangeUser( 1- downratio, 1+upratio)
        self.bottomframe.GetYaxis().SetNdivisions(202,False)
        for i in range(1, self.href.GetNbinsX()+1):
            self.bottomframe.SetBinContent(i,1)
            self.bottomframe.SetBinError(i,0)

        #self.bottomframe.SetLineColor(1)
        self.bottomframe.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
        self.bottomframe.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()*self.blabelratio)
        self.bottomframe.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize("Y")*(1+self.padratio)/(2*self.padratio))
        self.bottomframe.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*self.blabelratio)
        self.bottomframe.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*self.blabelratio)
        self.bottomframe.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset()/self.blabelratio)
        self.bottomframe.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength()*self.blabelratio)
        self.bottomframe.GetYaxis().SetTitle("Ratio")
        self.bottomframe.GetYaxis().CenterTitle()
        self.bottomframe.Draw("HIST")

        for h in self.ratiolist:
            h.Draw("SAME" + bottomopts)
class residualPlot(splitPlot):
    """Helper class to create a canvas with histograms above and residual plot below"""
    def __init__(self, usePull = False ):
        self.usePull = usePull
        splitPlot.__init__(self)

    def setref(self, ref):
        self.ref = ref.Clone("ref")

    def addhists(self, *h):
        try:
            self.hlist.extend(h)
        except:
            self.hlist = list(h)

    def reset(self):
        self.ref = None
        self.hlist = None

    def draw(self, topopts='', bottomopts='', runAutoCrop = False):

        self.pt.cd()

        self.ref.Draw(topopts)
        if runAutoCrop:
            autocropXaxis( [self.href,] + self.hlist )

        try:
            maxbin = max( [ h.GetBinContent(h.GetMaximumBin()) for h in self.hlist] + [self.ref.GetBinContent(self.ref.GetMaximumBin()),] )
        except:
            ##in case of function
            maxbin = max( [ h.GetBinContent(h.GetMaximumBin()) for h in self.hlist] + [self.ref.GetMaximum(),] )

        self.ref.GetXaxis().SetLabelSize(0)
        self.ref.GetYaxis().SetRangeUser(1e-6,maxbin*1.2)

        for h in self.hlist:
            h.Draw("SAME"+topopts)

        self.ref.Draw("SAME"+topopts)


        self.pb.cd()

        self.reslist = [ h.Clone(h.GetName() + "clone") for h in self.hlist ]
        self.pullplots = []
        for h in self.reslist:
            #h.GetXaxis().SetRange( self.href.GetXaxis().GetFirst(), self.href.GetXaxis().GetLast() )
            h.Add(self.ref,-1.0)

            if self.usePull:
                self.pullplots.append( ROOT.TH1F("pull"+h.GetName(),"pull",50,-5,5) )
                for b in range(1,h.GetNbinsX()+1):
                    h.SetBinContent( b, h.GetBinContent(b)/ h.GetBinError(b))
                    self.pullplots[-1].Fill( h.GetBinContent(b))
                    h.SetBinError( b, 1.0)

        maxres = max( [ h.GetBinContent(h.GetMaximumBin()) for h in self.reslist ] )
        minres = min( [ h.GetBinContent(h.GetMinimumBin()) for h in self.reslist ] )

        up = max( abs(maxres), abs(minres) )

        self.bottomframe = self.ref.Clone("pbframe")
        self.bottomframe.GetYaxis().SetRangeUser( -up, up)
        self.bottomframe.GetYaxis().SetNdivisions(202,False)
        for i in range(1, self.ref.GetNbinsX()+1):
            self.bottomframe.SetBinContent(i,0)
            self.bottomframe.SetBinError(i,0)

        #self.bottomframe.SetLineColor(1)
        self.bottomframe.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
        self.bottomframe.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()/self.padratio)
        self.bottomframe.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize("Y")*(1+self.padratio)/(2*self.padratio))
        self.bottomframe.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/self.padratio)
        self.bottomframe.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()/self.padratio)
        self.bottomframe.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset()*self.padratio)
        self.bottomframe.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength()/self.padratio)
        if self.usePull:
            self.bottomframe.GetYaxis().SetTitle("Pull")
        else:
            self.bottomframe.GetYaxis().SetTitle("Residual")
        self.bottomframe.GetYaxis().CenterTitle()
        self.bottomframe.Draw("HIST")

        for h in self.reslist:
            h.Draw("SAME"+bottomopts)
from math import sqrt
def efferr( nr, ntot ):
    return nr/ntot, sqrt( nr*((ntot-nr)**2) + (nr**2)*(ntot-nr) ) / ntot**2
