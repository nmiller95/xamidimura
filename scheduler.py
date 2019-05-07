"""
Simplest scheduler for Xamidimura telescopes.

1. Adjust preferences if necessary
2. Run file

"""

############################################# IMPORTS ##############################################

# astropy and astroplan modules for bulk of work
from astropy.time import Time
from astroplan import Observer

# logging module & config
import logging
logging.basicConfig(filename='log_files/scheduler.log', 
                    filemode='w',
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    level=logging.INFO)

# custom modules
from connect_database import connect_to, get_table_into_pandas
import feasibility_and_priority as f


####################################### SPECIFY PREFERENCES HERE #####################################

vzen = 21.7 # Zenith brightnesss at site (mag/arcsec^2)
k = 0.172 # Extinction coefficient at site
airmass_limit = 2.0 # Maximum allowable airmass
moon_tol = 5. # No. magnitudes of contrast acceptable between target & sky
del_pri = 0.2*0.2 # Change in priority required to mandate change of target

site = Observer.at_site('SAAO')
time = Time.now()
dbconn, dbcurs = connect_to()


############################################## DATABASE ##############################################

def make_eb_objects(conn, moon_tol):
    """
    Reads in database as pandas table
    
    :param conn:        Database connection
    :param moon_tol:    Smallest allowable sky/target brightness difference (mag)
    
    :returns:           List of targets as EclipsingBinary objects
    
    """
    t = get_table_into_pandas('target_info',conn)
    target_list = []
    for i in range(0,len(t)):
        name = t.loc[i][1]
        ra = t.loc[i][2]
        dec = t.loc[i][3]
        mag = t.loc[i][5]
        t0 = t.loc[i][7]
        period = t.loc[i][8]
        width_1 = t.loc[i][11]
        phase_2 = t.loc[i][13]
        width_2 = t.loc[i][12]
        temp = f.EclipsingBinary(name, ra, dec, mag, moon_tol, t0, period, width_1, phase_2, width_2)
        target_list.append(temp)
    logging.info("{} targets read from database".format(len(target_list)))
    return target_list


def observe_now(targets, time=None, current_target=None):
    """
    Main scheduler function to select best target at given time
    
    :param targets:         List of targets as EclipsingBinary objects
    :param time:            Astropy.time object
    :param current_target:  Target currently observed, as EclipsingBinary object
    
    :returns:               Best target (as EclipsingBinary object)
    
    """
    # UPDATES TIME TO CURRENT TIME, INTIALISATION
    if not time:
        time = Time.now()
    status, ecls = 0, 0
    frame = f.make_altaz_frame(site, time)
    
    # UPDATES PRIORITY SCORE FOR CURRENT TARGET, IF SPECIFIED AND FEASIBLE
    best_target, previous_best = None, 0
    if current_target:
        eclipsing = current_target.eclipse_status(time)
        if eclipsing:
            dark, moon_dim, airmass_ok = current_target.obs_feas(site, time, k, vzen, frame, airmass_limit)
            if dark and moon_dim and airmass_ok:
                mean_priority = current_target.mean_priority(site, time, frame, dbcurs)
                current_target.score = mean_priority
                previous_best = mean_priority
                best_target = current_target
    # IF CURRENT TARGET NOT FEASIBLE, REMOVE IT
            else: current_target = None
        else: current_target = None
    
    # ITERATES THROUGH LIST OF TARGETS
    for temp in targets:
    # FEASIBILITY MODEL - DETERMINES WHICH TARGETS ARE FEASIBLE
        eclipsing = temp.eclipse_status(time)
        if eclipsing:
            ecls += 1
            dark, moon_dim, airmass_ok = temp.obs_feas(site, time, k, vzen, frame, airmass_limit)
            if dark and moon_dim and airmass_ok:
                status += 1
    # SCORING MODEL - CALCULATES PRIORITY SCORES
                mean_priority = temp.mean_priority(site, time, frame, dbcurs)
                temp.score = mean_priority
    # SELECTION MODEL - CHOOSES BEST TARGET
                if mean_priority > previous_best:
                    best_target = temp
                    previous_best = mean_priority
    if best_target and current_target and best_target != current_target:
        if best_target.score < del_pri+current_target.score:
            best_target = current_target
    
    # LOGGING AND RETURNING BEST TARGET
    if best_target and best_target.score == 0.: # Will not observe if score=0
        best_target = None
        logging.info("Best target has score 0. Not observing")
    if status == 0:
        logging.info('{} targets eclipsing but zero are feasible'.format(ecls))
    else:
        logging.info('{} targets eclipsing, {} are feasible'.format(ecls,status))
    # EXTRA OUTPUT TO SCREEN
    print(best_target.name)
    return best_target


############################################## MAIN ##############################################

def main():
    targets = make_eb_objects(dbconn, moon_tol)
    current_target = None
    
    # Loops until user terminates program
    while True:
        cmd = input("Press enter to observe now, q to quit: ")
        if cmd == 'q':
            break
        else:
            current_target = observe_now(targets, current_target=current_target)
    
if __name__=="__main__":
    main()