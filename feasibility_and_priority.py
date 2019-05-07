"""
Functions for Xamidimura schedulers.

"""

############################################# IMPORTS ##############################################

# main modules for calculations
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

# logging module & config
import logging
logging.basicConfig(filename='log_files/feasibility_and_priority.log', 
                    filemode='w',
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    level=logging.INFO)

# custom modules
import plato_fields
from connect_database import match_target_name, match_target_id, update_priority_info

########################################## FUNCTIONS ############################################
def airmass_calc(Z):
    """
    Optical pathlength [Krisciunas & Schaefer 1991]
    
    :param Z: Zenith distance of the sky position in radians
    
    :returns: Optical pathlength in air masses
    
    """
    return (1-0.96*(np.sin(Z)**2))**(-0.5) #eq3

def sky_brightness_change(phase, sep, Z, Zm, k, vzen):
    """
    Change in sky brightness caused by the Moon [Krisciunas & Schaefer 1991]
    
    :param phase: Phase angle of the Moon (angular distance between Earth
                    and Sun as seen by the Moon)
    :param sep: Moon/sky seperation in degrees
    :param Z: Zenith distance of sky position in degrees
    :param Zm: Zenith distance of Moon in degrees
    :param k: Extinction coefficient
    :param vzen: Zenith brightness at location, in mag/arcsec^2
    
    :returns: Change in sky brightness in mag/arcsec^2
    
    """
    sep_r = np.deg2rad(sep*u.deg)
    Z = np.deg2rad(Z*u.deg)
    Zm = np.deg2rad(Zm*u.deg)
    
    # Zenith brightness in nL
    Bzen = 34.08*np.exp(20.7233-0.92104*vzen) #eq1
    # Illuminance of moon outside of atmosphere
    illum = 10**(-0.4*(3.84+0.026*abs(phase)+4e-9*phase**4)) #eq20
    # Rayleigh scattering
    scat_r = 10**5.36*(1.06+np.cos(sep_r)**2) #eq17
    # Mie scattering
    if sep >= 10:
        scat_m = 10**(6.15-sep/40.) #eq18
    else:
        #scat_m = 6.2e7*sep**(-2) #eq19
        return(-20)
    # Total scattering
    scat = scat_r + scat_m #eq16
    
    # Model surface brightness of Moon
    Bmoon = scat*illum*10**(-0.4*k*airmass_calc(Zm))*(1-10**(-0.4*k*airmass_calc(Z))) #eq15
    # Dark nighttime sky brightness as a function of zenith distance
    B0 = Bzen*10**(-0.4*k*(airmass_calc(Z)-1))*airmass_calc(Z) #eq2 
    # Change in V band sky brightness caused by moonlight
    delv = -2.5*np.log10((Bmoon+B0)/B0) #eq22
    return delv

def daylight_test(location, time):
    """
    Boolean report of darkness
    
    :param location: astroplan.observer object
    :param time: astropy.Time object
    
    :returns: True for dark (sun below -12), False for light
    
    """
    return location.is_night(time)

def make_altaz_frame(location, time):
    """
    Turns location and time into altaz frame
    
    :param location: astroplan.observer object
    :param time: astropy.Time object
    
    :returns: Altaz frame
    
    """
    return location.altaz(time)

def get_eclipse_type_from_phase(target_name, phase, curs):
    """
    General method for calculating eclipse type from any phase
    
    :param target_name: String name of target as appears in target_info table
    :param phase: Float value for phase (0-1)
    :param curs: Database cursor
    
    :returns: 1 for primary, 2 for secondary, 0 for not-in-eclipse
    
    """
    info = match_target_name(target_name,'target_info',curs)
    width_1 = info[0][11]
    width_2 = info[0][12]
    phase_2 = info[0][13]
    
    primary_start = 1-0.6*width_1
    primary_end   = 0.6*width_1
    secondary_start = phase_2-0.6*width_2
    secondary_end = phase_2+0.6*width_2
    
    if (phase>=primary_start and phase<=1) or (phase<=primary_end and phase>=0):
        return 1
    elif (phase>=secondary_start and phase<=phase_2) or (phase<=secondary_end and phase>=phase_2):
        return 2
    else:
        return 0

def get_bin_number(target_name, phase, curs):
    """
    Calculates which phase bin an observation falls in
    
    :param target_name: String name of target as appears in target_info table
    :param phase: Float value for phase (0-1)
    :param curs: Database cursor
    
    :returns: Integer bin number
    
    """
    ecl_type = get_eclipse_type_from_phase(target_name, phase, curs)
    info = match_target_name(target_name,'target_info',curs)
    target_id = info[0][0]
    pri_info = match_target_id(target_id,'priority_table',curs)
    
    if ecl_type == 1:
        width_1 = info[0][11]
        primary_start = 1-0.6*width_1
        primary_end   = 0.6*width_1
        # Calculates parameters of the phase bins
        no_bins = pri_info[0][2]
        ecl_len = (primary_end - primary_start)%1
        bin_width = ecl_len/no_bins
        bin_start = primary_start
        # Returns the bin number for the current phase
        for n in range(0,no_bins):
            bin_end = bin_start + bin_width
            if bin_end>= 1:
                bin_end -= 1
                if (phase>=bin_start and phase<1) or (phase<bin_end and phase>=0):
                    return n
            elif (phase>=bin_start and phase<bin_end):
                return n
            bin_start += bin_width
            if bin_start>1: bin_start-=1
    elif ecl_type == 2:
        width_2 = info[0][12]
        phase_2 = info[0][13]
        secondary_start = phase_2-0.6*width_2
        secondary_end = phase_2+0.6*width_2
        # Calculates parameters of the phase bins
        no_bins = pri_info[0][5]
        ecl_len = (secondary_end - secondary_start)%1
        bin_width = ecl_len/no_bins
        bin_start = secondary_start
        # Returns the bin number for the current phase
        for n in range(0,no_bins):
            bin_end = bin_start + bin_width
            if phase>=bin_start and phase<bin_end:
                return n
            bin_start += bin_width

def increment_completeness(target_name, phase, curs, conn):
    """
    Increments the completeness value for particular phase bin by one
    
    :param target_name: String name of target as appears in target_info table
    :param phase: Float value for phase (0-1)
    :param curs: Database cursor
    :param conn: Database connection
    
    """
    ecl_type = get_eclipse_type_from_phase(target_name, phase, curs)
    bin_no   = get_bin_number(target_name, phase, curs)
    if bin_no == None:
        logging.warning('No bin number found for {}'.format(target_name))
        return 0
    
    info = match_target_name(target_name,'target_info',curs)
    target_id = info[0][0]
    pri_info = match_target_id(target_id,'priority_table',curs)
    
    if ecl_type == 1:
        completeness = pri_info[0][4]
        col_name = 'COMPLETENESS_1'
    elif ecl_type == 2:
        completeness = pri_info[0][7]
        col_name = 'COMPLETENESS_2'
    else:
        return 0
    
    # Increments the bin value by one
    bin_value = int(completeness[2*bin_no:2*(bin_no+1)],16)
    bin_value = hex(bin_value + 1)
    if len(bin_value) == 3: # is of format '0x[value]'
        completeness = completeness[:2*bin_no+1] + bin_value[2] + completeness[2*(bin_no+1):]
    elif len(bin_value) == 4:
        completeness = completeness[:2*bin_no] + bin_value[2:4] + completeness[2*(bin_no+1):]
    else:
        print('Invalid entry - number of observations in bin exceeds 256')
        pass
    update_priority_info(target_id, col_name, completeness, curs, conn)


################################### TARGET CLASS --- PARENT #####################################
class Target:
    def __init__(self, name, ra, dec, mag, moon_tol):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.mag = mag
        self.moon_tol = moon_tol
    
    @property
    def coords(self):
        """
        Converts hmsdms RA and Dec into skycoord object
        
        :returns: RA/Dec skycoord object
        
        """
        return SkyCoord(self.ra,self.dec,unit=(u.hourangle,u.deg))
     
    def altaz_transform(self, frame):
        """
        Transforms RA/Dec into AltAz format
        
        :param frame: Altaz frame

        :returns: AltAz skycoord object
        
        """
        return self.coords.transform_to(frame)
    
    def airmass(self, frame):
        """
        Modified sec(z) airmass calculation
        
        :param frame: Altaz frame
        
        :returns: Airmass
        
        """
        Z = np.deg2rad(90*u.deg - self.altaz_transform(frame).alt)
        return (1/np.cos(Z)) - 0.010*(1/np.cos(Z) - 1)**2

    def secz(self, frame):
        """
        Classic sec(z) airmass calculation
        
        :param frame: Altaz frame
        
        :returns: Airmass
        
        """
        Z = np.deg2rad(90*u.deg - self.altaz_transform(frame).alt)
        return 1/np.cos(Z)
              
    def sky_brightness(self, location, time, k, vzen, frame):
        """
        Sky brightness due to moonlight calculation
        
        :param location: astroplan.observer object
        :param time: astropy.Time object
        :param k: Extinction coefficient
        :param vzen: Zenith brightness at location, in mag/arcsec^2
        :param frame: Altaz frame
        
        :returns: Sky brightness in visual magnitudes
        
        """
        moon = location.moon_altaz(time)
        Zm = 90 - moon.alt/u.deg
        phase = np.rad2deg(location.moon_phase(time))/u.deg
        Z = 90 - self.altaz_transform(frame).alt/u.deg
        sep = moon.separation(self.coords)/u.deg
        return vzen+sky_brightness_change(phase,sep,Z,Zm,k,vzen)
    
    ############################### feasibility ###############################
    
    def airmass_test(self, frame, airmass_limit):
        """
        Boolean test for acceptable airmass
        
        :param frame: Altaz frame
        :param airmass_limit: Maximum airmass tolerated
        
        :returns: True if airmass acceptable, False otherwise
        
        """
        airmass = self.airmass(frame)
        if airmass >= 1 and airmass <= airmass_limit:    return True
        else:    return False
            
    def sky_contrast_test(self, location, time, k, vzen, frame):
        """
        Boolean test for acceptable sky contrast
        
        :param location: astroplan.observer object
        :param time: astropy.Time object
        :param k: Extinction coefficient
        :param vzen: Zenith brightness at location, in mag/arcsec^2
        :param frame: Altaz frame
        
        :returns: True if sky contrast acceptable, False otherwise
        
        """
        contrast = self.sky_brightness(location, time, k, vzen, frame)-self.mag
        if self.moon_tol < contrast:    return True
        else:    return False
    
    def obs_feas(self, location, time, k, vzen, frame, airmass_limit):
        """
        Reports on current observational feasibility for ANY target type
        
        :param location: astroplan.observer object
        :param time: astropy.Time object
        :param k: Extinction coefficient
        :param vzen: Zenith brightness at location, in mag/arcsec^2
        :param frame: Altaz frame
        :param airmass_limit: Maximum airmass tolerated
        
        :returns: Boolean score for darkness, moon sky contrast and airmass
        
        """
        dark = daylight_test(location, time)
        moon_dim = self.sky_contrast_test(location, time, k, vzen, frame)
        airmass_ok = self.airmass_test(frame, airmass_limit)
        return dark, moon_dim, airmass_ok
    
    
    ############################## priority ###############################
    
    @property
    def plato_score(self):
        """
        Calculates PLATO field proximity score
        
        :returns: Score (float, 0.1-1.0)
        
        """
        plato = plato_fields.calculate(self.coords)
        logging.debug(f"PLATO field calculation: {plato}")
        if (plato[0] or plato[1]) and plato[4]>10:      return 1.0
        if (plato[0] or plato[1]) and 5<plato[4]<10:    return 0.9
        if (plato[0] or plato[1]) and plato[4]<5:       return 0.7
        if not(plato[0] or plato[1]) and plato[4]<5:    return 0.5
        if not(plato[0] or plato[1]) and 5<plato[4]<10: return 0.3
        if not(plato[0] or plato[1]) and plato[4]>10:   return 0.1
    
    def airmass_score(self, location, frame):
        """
        Calculates airmass score
        
        :param location: astroplan.observer object
        :param frame: Altaz frame
        
        :returns: Score (float, 0.0-1.0)
        
        """
        current_airmass = self.airmass(frame)
        if current_airmass < 1:
            return 0
        max_elevation = 90*u.deg - abs(location.location.lat-self.coords.dec)
        Z = np.deg2rad(90*u.deg - max_elevation)
        best_airmass = (1/np.cos(Z)) - 0.010*(1/np.cos(Z) - 1)**2
        return max(((2-current_airmass)/(2-best_airmass)),0)    
    
    def __repr__(self):
        return "Target('{}')".format(self.name)
    
    
###################################### ECLIPSING BINARY CLASS #######################################
class EclipsingBinary(Target):
    def __init__(self, name, ra, dec, mag, moon_tol, t0, period, 
                 width_1, phase_2, width_2):
        super().__init__(name, ra, dec, mag, moon_tol)
        self.t0 = t0
        self.period = period 
        self.width_1 = width_1
        self.phase_2 = phase_2 
        self.width_2 = width_2
        self.primary = 0
        self.secondary = 0
        self.phase = None
        self.score = 0
    
    def primary_eclipse(self, time):
        """
        Boolean report of primary eclipse status
        
        :param time: astropy.Time object
        
        :returns: True if primary eclipse is occurring, False otherwise
        
        """
        current_time = time.jd
        current_phase = ((current_time-self.t0)/self.period)%1
        phase_start = 1-0.6*self.width_1
        phase_end   = 0+0.6*self.width_1
        if current_phase >= phase_start and current_phase <= 1.:
            return True, current_phase
        elif current_phase <= phase_end and current_phase >= 0.:
            return True, current_phase
        else:
            return False, current_phase
    
    def secondary_eclipse(self, time):
        """
        Boolean report of secondary eclipse status
        
        :param time: astropy.Time object
        
        :returns: True if secondary eclipse is occurring, False otherwise
        
        """
        current_time = time.jd
        current_phase = ((current_time-self.t0)/self.period)%1
        phase_start = self.phase_2-0.6*self.width_2
        phase_end   = self.phase_2+0.6*self.width_2
        if current_phase >= phase_start and current_phase <= phase_end:
            return True, current_phase
        else:
            return False, current_phase
    
    ############################### feasibility ###############################
    
    def eclipse_status(self, time):
        """
        Reports current eclipsing status of target
        
        :param time: astropy.Time object
       
        :returns 1 if primary eclipse, 2 if secondary eclipse, 0 if out of eclipse
        
        """
        self.primary, self.phase = self.primary_eclipse(time)
        self.secondary, self.phase = self.secondary_eclipse(time)
        if self.primary and self.secondary:
            logging.warning("You broke physics again >_>")
            pass
        elif self.primary:    return 1
        elif self.secondary:  return 2
        else:                 return 0
        
    ############################## priority ###############################
    
    def phase_score(self, time, curs):
        """
        Calculation of phase score for current bin
        
        :param time: astropy.Time object
        :param curs: Database cursor
        
        :returns: Phase score (prioritises complete eclipse coverage)
        
        """
        ecl_type = self.eclipse_status(time)
        bin_no   = get_bin_number(self.name, self.phase, curs)
        
        info = match_target_name(self.name,'target_info',curs)
        target_id = info[0][0]
        pri_info = match_target_id(target_id,'priority_table',curs)
        
        if ecl_type == 1:
            no_bins = pri_info[0][2]
            obs_per_bin = pri_info[0][3]
            completeness = pri_info[0][4]
        elif ecl_type == 2:
            no_bins = pri_info[0][5]
            obs_per_bin = pri_info[0][6]
            completeness = pri_info[0][7]
        else:
            return 0
        
        bin_value = int(completeness[2*bin_no:2*(bin_no+1)],16)
        gradient = 0.5 / (no_bins-3)#2) # highest score when 1 observation away from complete
        intercept = 0.5 # unobserved score = 0.5
        if int(bin_value) > obs_per_bin-1:    phase_score = 0.0 #obs_per_bin
        elif int(bin_value) == obs_per_bin-1: phase_score = 0.5 #obs_per_bin
        else: phase_score = int(bin_value)*gradient + intercept
        return(phase_score)
    
    def mean_priority(self, location, time, frame, curs):
        """
        Calculates each priority and returns average
        
        :param location: astroplan.observer object
        :param time: astropy.Time object
        :param frame: Altaz frame
        :param curs: Database cursor
        
        """
        priority_airmass = self.airmass_score(location, frame)
        priority_plato   = self.plato_score
        priority_phase   = self.phase_score(time, curs)
        weight_airmass   = 0.5
        weight_plato     = 0.5
        weight_phase     = 2
        self.score = priority_airmass**weight_airmass * priority_plato**weight_plato *  priority_phase**weight_phase
        return(self.score)
    
    def __repr__(self):
        return "EclipsingBinary('{}')".format(self.name)