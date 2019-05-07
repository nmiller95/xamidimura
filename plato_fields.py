from numpy import polyfit, polyval, array, linspace, pi, interp 
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.modeling.projections import Sky2Pix_HammerAitoff 
from astropy.modeling.projections import Pix2Sky_HammerAitoff 
from matplotlib.path import Path

# Declares nominal locations of centre of south and north fields
c_SPF = SkyCoord(253, -30, unit='deg',frame='galactic')
c_NPF = SkyCoord(65, 30, unit='deg',frame='galactic')


# Maps image pixels to galactic coordinates 
Sky2Pix = Sky2Pix_HammerAitoff
Pix2Sky = Pix2Sky_HammerAitoff

tcal = Table.read('plato_info/coords_gal.dat',format='ascii',names=['l','b','x','y'])

xd,yd = Sky2Pix.evaluate(((tcal['l'] + 180) % 360) -180,tcal['b'])
px = polyfit(xd,tcal['x'],1)
py = polyfit(yd,tcal['y'],1)
qx = polyfit(tcal['x'],xd,1)
qy = polyfit(tcal['y'],yd,1)

tspf = Table.read('plato_info/SPF_gal_pxl.dat',format='ascii',names=['x','y'])
xd = polyval(qx,tspf['x'])
yd = polyval(qy,tspf['y'])
c_SPF_edge = SkyCoord(array(Pix2Sky.evaluate(xd,yd)).T, unit='deg',frame='galactic')

tnpf = Table.read('plato_info/NPF_gal_pxl.dat',format='ascii',names=['x','y'])
xd = polyval(qx,tnpf['x'])
yd = polyval(qy,tnpf['y'])
c_NPF_edge = SkyCoord(array(Pix2Sky.evaluate(xd,yd)).T, unit='deg',frame='galactic')


def calculate(coords):
    """
    Calculates relative location of given target to nominal PLATO fields
    
    :param coords: Skycoord object of target
    
    :returns: in_npf (bool), in_spf (bool), closest (field, "NPF" or "SPF"),
        d_cen (dist to centre of closest field in deg), d_edg (dist to edge of
        closest field in deg)
    
    """
    # Recalculates target coords in terms of PLATO field images
    xd,yd = Sky2Pix.evaluate([((coords.galactic.l.deg + 180) % 360) -180],[coords.galactic.b.deg])
    xpxl = polyval(px,xd)
    ypxl = polyval(py,yd)

    path_spf = Path(array([tspf['x'],tspf['y']]).T)
    path_npf = Path(array([tnpf['x'],tnpf['y']]).T)

    # Determines if target is in NPF or SPF or neither
    if path_npf.contains_points(array([xpxl,ypxl]).T):
        in_npf = True
    else: in_npf = False
    if path_spf.contains_points(array([xpxl,ypxl]).T):
        in_spf = True
    else: in_spf = False

    # Interpolates edge coordinates onto a finer grid
    phigrid = linspace(0,2*pi,3601)

    phi_spf = c_SPF_edge.position_angle(c_SPF)
    l_fine_spf = interp(phigrid,phi_spf,c_SPF_edge.l,period=2*pi)
    b_fine_spf = interp(phigrid,phi_spf,c_SPF_edge.b,period=2*pi)
    c_spf_fine = SkyCoord(l_fine_spf, b_fine_spf, unit='deg',frame='galactic')

    phi_npf = c_NPF_edge.position_angle(c_NPF)
    l_fine_npf = interp(phigrid,phi_npf,c_NPF_edge.l,period=2*pi)
    b_fine_npf = interp(phigrid,phi_npf,c_NPF_edge.b,period=2*pi)
    c_npf_fine = SkyCoord(l_fine_npf, b_fine_npf, unit='deg',frame='galactic')

    d_NPF = c_NPF.separation(coords)
    d_SPF = c_SPF.separation(coords)

    # Calculates distances to centre/edges of closest field
    if d_NPF < d_SPF:
        closest = "NPF"
        d_cen   = d_NPF.deg
        d_edg   = min(c_npf_fine.separation(coords)).deg
    else:
        closest = "SPF"
        d_cen   = d_SPF.deg
        d_edg   = min(c_spf_fine.separation(coords)).deg
    
    return (in_npf, in_spf, closest, d_cen, d_edg)