''' Module for Gilerson Chl algorithms (rrs-based and rhos-based)

author: jpscott6, 2018.05.08

'''

#==========================================================================================================================================

def calc_kd_rhos(rhos_442, rhos_490, rhos_620, rhos_665, rhos_865):

    #Kd calculation with 865 correction applied
    kd = (0.7 * ((((rhos_620 + rhos_665) / 2.) - rhos_865) / (((rhos_442 + rhos_490) / 2.) - rhos_865)))

    #apply calibration
    kd = (4.0 * kd) - 0.69

    return kd

#==========================================================================================================================================

def calc_chl_rhos(rhos_442, rhos_490, rhos_620, rhos_665, rhos_709, rhos_865):

    #Gilerson chl Rho_s
    if (rhos_665 - rhos_865) == 0.0:
        return float('nan')

    else:
        rel_chl = (rhos_709 - rhos_865) / (rhos_665 - rhos_865)

        if rel_chl <= 19.30/35.75:
            rel_chl = (19.30/35.75)

        chl = ((35.75 * rel_chl) - 19.30) ** 1.124
    
    #calculate and apply a correction in clear water based on kd where kd_min_threshold = 0.31
    kd_min_threshold = 0.31
    kd = calc_kd_rhos(rhos_442, rhos_490, rhos_620, rhos_665, rhos_865)            

    if (kd < kd_min_threshold) & (kd > 0):
        chl = float('nan')

    return chl

#==========================================================================================================================================

def calc_kd_rrs(rrs_442, rrs_490, rrs_620, rrs_665):

    #Kd calculation
    kd = (0.7 * (((rrs_620 + rrs_665) / 2.) / ((rrs_442 + rrs_490) / 2.)))

    #apply calibration
    kd = (4.0 * kd) - 0.69

    return kd

#==========================================================================================================================================

def calc_chl_rrs(rrs_442, rrs_490, rrs_620, rrs_665, rrs_709):

    #Gilerson chl rrs
    if rrs_665 == 0.0:
        return float('nan')

    else:
        rel_chl = (rrs_709 / rrs_665)

        if rel_chl <= 19.30/35.75:
            rel_chl = (19.30/35.75)

        chl = ((35.75 * rel_chl) - 19.30) ** 1.124
    
    #calculate and apply a correction in clear water based on kd where kd_min_threshold = 0.31
    kd_min_threshold = 0.31
    kd = calc_kd_rrs(rrs_442, rrs_490, rrs_620, rrs_665)            

    if (kd < kd_min_threshold) & (kd > 0):
        chl = float('nan')
 
    return chl

#==========================================================================================================================================

def calc_chl_rrs_2band(rrs_665, rrs_709):

    #Gilerson 2010, advanced 2-band chl via rrs
    if rrs_665 == 0.0:
        return float('nan')

    else:
        rel_chl = (rrs_709 / rrs_665)

        if rel_chl <= 19.30/35.75:
            rel_chl = (19.30/35.75) #keeps chl from becoming complex if base is negative

        chl = ((35.75 * rel_chl) - 19.30) ** 1.124
 
    return chl

#==========================================================================================================================================

def calc_chl_rhos_2band(rhos_665, rhos_709, rhos_865):

    #Gilerson chl Rho_s
    if (rhos_665 - rhos_865) == 0.0:
        return float('nan')

    else:
        rel_chl = (rhos_709 - rhos_865) / (rhos_665 - rhos_865)

        if rel_chl <= 19.30/35.75:
            rel_chl = (19.30/35.75)

        chl = ((35.75 * rel_chl) - 19.30) ** 1.124

    return chl
