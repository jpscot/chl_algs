''' Module for MERIS MPH Chl algorithm

Maximum Peak Height of chlorophyll for MERIS

Mark William Matthews, Daniel Odermatt
Remote Sensing of Environment 156 (2015) 374–382

author: jpscott6, 2018.05.30

'''

#==========================================================================================================================================

def calc_MPH_chl(rhos620, rhos664, rhos681, rhos709, rhos753, rhos885):

    from numpy import exp

    #wav6 = 620
    #wav7 = 664
    #wav8 = 681
    #wav9 = 709
    #wav10 = 753
    #wav14 = 885

#    if (rhos620 <= 0.0):
        #if(Rrs620 <= 0.0):
        #if(Rrs620 <= 0.0 || Rrs664 <= 0.0 || Rrs681 <= 0.0 || Rrs709 <= 0.0):
#        chl_mph = float('NaN')
        #l2rec->l1rec->flags[ip] |= PRODFAIL

#    else:

    if (rhos681 > rhos709):
        wavmax0 = 681
        Rmax0 = rhos681
    else:
        wavmax0 = 709
        Rmax0 = rhos709

    if (Rmax0 > rhos753):
        wavmax1 = wavmax0
        Rmax1 = Rmax0
    else:
        wavmax1 = 753
        Rmax1 = rhos753

    #normalised difference vegetation index
    ndvi = (rhos885 - rhos664) / (rhos885 + rhos664)

    #sun-induced phycocyanin absorption fluorescence
    sipf = rhos664 - rhos620 - (rhos681 - rhos620) * (664 - 619) / (681 - 619)

    #sun induced chlorophyll fluorescence
    sicf = rhos681 - rhos664 - (rhos709 - rhos664) * (681 - 664) / (709 - 664)

    #backscatter and absorption induced reflectance
    bair = rhos709 - rhos664 - (rhos885 - rhos664) * (709 - 664) / (885 - 664)

    mph0 = Rmax0 - rhos664 - (rhos885 - rhos664) * (wavmax0 - 664) / (885 - 664)
    mph1 = Rmax1 - rhos664 - (rhos885 - rhos664) * (wavmax1 - 664) / (885 - 664)

    if (wavmax1 != 753):

        if (sicf >= 0 or sipf <= 0 or bair <= 0.002):
            #immersed eukaryotes
            #R.Healy had 3.02e3 instead of 4.02e3 in the paper
            chl_mph = 5.24e9 * mph0**4 - 1.95e8 * mph0**3 + 2.46e6 * mph0**2 + 4.02e3 * mph0 + 1.97

        else:
            #cyanobacteria, immersed or floating
            chl_mph = 22.44 * exp(35.79 * mph1)

    else:

        if (mph1 >= 0.02 or ndvi >= 0.2):

            if (sicf < 0 and sipf > 0):
                #cyanobacteria, immersed or floating
                chl_mph = 22.44 * exp(35.79 * mph1)

            else:
                #floating vegetation
                chl_mph = float('NaN')

        else:
            #immersed eukaryotes
            #R.Healy had 3.02e3 instead of 4.02e3 in the paper
            chl_mph = 5.24e9 * mph0**4 - 1.95e8 * mph0**3 + 2.46e6 * mph0**2 + 4.02e3 * mph0 + 1.97

#    if (chl_mph < 0.):
#        chl_mph = 0.

    return chl_mph

#==========================================================================================================================================

def calc_MPH_chl_865(rhos620, rhos664, rhos681, rhos709, rhos753, rhos865):

    from numpy import exp

    #wav6 = 620
    #wav7 = 664
    #wav8 = 681
    #wav9 = 709
    #wav10 = 753
    #wav14 = 865

#    if (rhos620 <= 0.0):
        #if(Rrs620 <= 0.0):
        #if(Rrs620 <= 0.0 || Rrs664 <= 0.0 || Rrs681 <= 0.0 || Rrs709 <= 0.0):
#        chl_mph = float('NaN')
        #l2rec->l1rec->flags[ip] |= PRODFAIL

#    else:

    if (rhos681 > rhos709):
        wavmax0 = 681
        Rmax0 = rhos681
    else:
        wavmax0 = 709
        Rmax0 = rhos709

    if (Rmax0 > rhos753):
        wavmax1 = wavmax0
        Rmax1 = Rmax0
    else:
        wavmax1 = 753
        Rmax1 = rhos753

    #normalised difference vegetation index
    ndvi = (rhos865 - rhos664) / (rhos865 + rhos664)

    #sun-induced phycocyanin absorption fluorescence
    sipf = rhos664 - rhos620 - (rhos681 - rhos620) * (664 - 619) / (681 - 619)

    #sun induced chlorophyll fluorescence
    sicf = rhos681 - rhos664 - (rhos709 - rhos664) * (681 - 664) / (709 - 664)

    #backscatter and absorption induced reflectance
    bair = rhos709 - rhos664 - (rhos865 - rhos664) * (709 - 664) / (865 - 664)

    mph0 = Rmax0 - rhos664 - (rhos865 - rhos664) * (wavmax0 - 664) / (865 - 664)
    mph1 = Rmax1 - rhos664 - (rhos865 - rhos664) * (wavmax1 - 664) / (865 - 664)

    if (wavmax1 != 753):

        if (sicf >= 0 or sipf <= 0 or bair <= 0.002):
            #immersed eukaryotes
            #R.Healy had 3.02e3 instead of 4.02e3 in the paper
            chl_mph = 5.24e9 * mph0**4 - 1.95e8 * mph0**3 + 2.46e6 * mph0**2 + 4.02e3 * mph0 + 1.97

        else:
            #cyanobacteria, immersed or floating
            chl_mph = 22.44 * exp(35.79 * mph1)

    else:

        if (mph1 >= 0.02 or ndvi >= 0.2):

            if (sicf < 0 and sipf > 0):
                #cyanobacteria, immersed or floating
                chl_mph = 22.44 * exp(35.79 * mph1)

            else:
                #floating vegetation
                chl_mph = float('NaN')

        else:
            #immersed eukaryotes
            #R.Healy had 3.02e3 instead of 4.02e3 in the paper
            chl_mph = 5.24e9 * mph0**4 - 1.95e8 * mph0**3 + 2.46e6 * mph0**2 + 4.02e3 * mph0 + 1.97

#    if (chl_mph < 0.):
#        chl_mph = 0.

    return chl_mph
