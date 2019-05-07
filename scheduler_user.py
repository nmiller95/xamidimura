"""
Comprehensive scheduler for Xamidimura telescopes with optional user interface.

1. Adjust preferences if necessary, choose interface = True or False
2. Run file
3. Follow instructions from interface
4. See readme for info on functions if not using interface

"""

############################################# IMPORTS ##############################################

# basic modules
import numpy as np
import sys
import matplotlib.pyplot as plt

# astropy and astroplan modules for bulk of work
from astropy.time import Time
from astroplan import Observer
import astropy.units as u

# modules for displaying output
from prettytable import PrettyTable
from astroplan.plots import plot_sky
from astroplan import FixedTarget

# logging module & config
import logging
logging.basicConfig(filename='log_files/scheduler_user.log', 
                    filemode='w',
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    level=logging.INFO)

# custom modules
from connect_database import connect_to, get_table_into_pandas
import feasibility_and_priority as f


####################################### SPECIFY PREFERENCES HERE #####################################

# RUN INTERFACE OR JUST ACCESS FUNCTIONS
interface = True

site = Observer.at_site('SAAO')
time = Time.now()

vzen = 21.7 # Zenith brightnesss at site (mag/arcsec^2)
k = 0.172 # Extinction coefficient at site
airmass_limit = 2.0 # Maximum allowable airmass
moon_tol = 5. # No. magnitudes of contrast acceptable between target & sky
del_pri = 0.2*0.2 # Change in priority required to mandate change of target


############################################## DATABASE ##############################################

dbconn, dbcurs = connect_to()
t = get_table_into_pandas('target_info',dbconn)

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

    
######################################### FUNCTIONS ############################################
    
def why_not(target_name, target_list=target_list, time=time):
    """
    Report of feasibility status for specified target
    
    :param target_name: String target name as it appears in the database. 
                        Can also use EclipsingBinary.name
    :param target_list: List of targets as EclipsingBinary objects
    :param time:        Astropy.time object
    
    :returns:           Eclipsing status, darkness, moon brightness contrast,
                        and airmass feasibility scores.
    
    """
    frame = f.make_altaz_frame(site, time)
    for item in target_list:
        if item.name == target_name:
            eclipsing = item.eclipse_status(time)
            dark, moon_dim, airmass_ok = item.obs_feas(site, time, k, vzen, frame, airmass_limit)
            print('Feasibility status report for {}\n'.format(item.name))
            print('Eclipsing status:             {}'.format(eclipsing))
            print('Suitable darkness:            {}'.format(dark))
            print('Moon brightness contrast:     {}'.format(moon_dim))
            print('Acceptable airmass:           {}\n'.format(airmass_ok))
            return (eclipsing, dark, moon_dim, airmass_ok)


def observe_now(time=time, targets=target_list, current_target=None):
    """
    Main scheduler function to select best target at given time
    
    :param time:            Astropy.time object
    :param targets:         List of targets as EclipsingBinary objects
    :param current_target:  Target currently observed, as EclipsingBinary object
    
    :returns:               Best target (as EclipsingBinary object)
    
    """
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
    for temp in target_list:
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
    return (best_target, status)


def feasible_now(time=time, targets=target_list):
    """
    Returns all feasible eclipsing targets at a given time
    
    :param time:            Astropy.time object
    :param targets:         List of targets as EclipsingBinary objects
    
    :returns:               List of targets as EclipsingBinary objects
    
    """
    frame = f.make_altaz_frame(site, time)
    status, ecls = 0, 0
    feasibles = []
    # ITERATES THROUGH LIST OF TARGETS
    for temp in target_list:
    # FEASIBILITY MODEL - DETERMINES WHICH TARGETS ARE FEASIBLE
        eclipsing = temp.eclipse_status(time)
        if eclipsing:
            ecls += 1
            dark, moon_dim, airmass_ok = temp.obs_feas(site, time, k, vzen, frame, airmass_limit)
            if dark and moon_dim and airmass_ok:
                status += 1
                feasibles.append(temp)
    # LOGGING AND RETURNING BEST TARGET
    if status == 0:
        logging.info('{} targets eclipsing but zero are feasible'.format(ecls))
    else:
        logging.info('{} targets eclipsing, {} are feasible'.format(ecls,status))
    return feasibles


def eclipse_info(target_name, time=time, display=True):
    """
    Report of eclipse times and airmasses for targets *currently in eclipse*
    
    :param target_name: String target name as it appears in the database. 
                        Can also use EclipsingBinary.name
    :param time:        Astropy.time object
    
    :returns:           Eclipsing status, eclipsing start time, eclipsing
                        end time, altitude at eclipse start and at end, 
                        priority score
    
    """
    frame = f.make_altaz_frame(site, time)
    for item in target_list:
        if item.name == target_name:
            eclipsing = item.eclipse_status(time)
            if eclipsing:
                # eclipse start and end time
                current_phase = item.phase
                if eclipsing == 1:
                    phase_start = 1 - 0.6*item.width_1
                    phase_end = 0 + 0.6*item.width_1
                    if current_phase < 0.1:
                        phase_del_minus = abs(phase_start - current_phase+1)
                        phase_del_plus  = (phase_end - current_phase)
                    else:
                        phase_del_minus = abs(phase_start - current_phase)
                        phase_del_plus  = (phase_end+1 - current_phase)
                elif eclipsing == 2:
                    phase_start = item.phase_2 - 0.6*item.width_2
                    phase_end = item.phase_2 + 0.6*item.width_2
                    phase_del_minus = abs(phase_start - current_phase)
                    phase_del_plus  = abs(phase_end - current_phase)
                time_start = time - (phase_del_minus*item.period)*u.day
                time_end = time + (phase_del_plus*item.period)*u.day
                # altitude and airmasses
                frame_start = f.make_altaz_frame(site, time_start)
                frame_end = f.make_altaz_frame(site, time_end)
                alt_start = item.altaz_transform(frame_start).alt
                alt_end   = item.altaz_transform(frame_end).alt
                alt_now   = item.altaz_transform(frame).alt
                airmass_now   = item.airmass(frame)
                # priority scores
                priority_score = round(float(item.mean_priority(site, time, frame, dbcurs)),ndigits=3)
                if display:
                    # for use in detailed eclipse info function
                    return(
                    print('Status report for {}\n'.format(item.name)),
                    print('Eclipsing status:           {}'.format(eclipsing)),
                    print('Current phase:              {}'.format(round(float(item.phase),ndigits=3))),
                    print('Current airmass:            {}'.format(round(float(airmass_now),ndigits=3))),
                    print('Current priority score:     {}'.format(round(float(priority_score),ndigits=3))),
                    print('Eclipse started:            {}'.format(time_start)),#.isot[11:19])),
                    print('Eclipse ends:               {}'.format(time_end)),#.isot[11:19])),
                    print('Altitude @ start:           {} deg'.format(round(float(alt_start/u.deg),ndigits=3))),
                    print('Current altitude:           {} deg'.format(round(float(alt_now/u.deg),ndigits=3))),
                    print('Altitude @ end:             {} deg\n'.format(round(float(alt_end/u.deg),ndigits=3)))
                    )
                else:
                    # for use in tonight's best function
                    return (eclipsing, time_start, time_end, alt_start, alt_end, priority_score)
            else:
                return None


def tonights_best(time=time, targets=target_list):
    """
    List of tonight's best targets
    
    :param time: Astropy.time object
    :param targets: List of targets as EclipsingBinary objects
    
    :returns: target_array - List of targets as eclipsing binary objects
              info         - Eclipsing status, eclipsing start time,  eclipsing 
                              end time, altitude at eclipse start and at end
    
    """
    twi_pm = site.twilight_evening_astronomical(time,which='nearest')
    twi_am = site.twilight_morning_astronomical(time,which='next')
    times = np.linspace(twi_pm.jd, twi_am.jd, 10)
    target_array, info = [], []
    # LOOP THROUGH TIMES AND COLLECT FEASIBLE TARGETS
    for i in times:
        i = Time(i,format='jd')
        feas_targets = feasible_now(i, targets)
        for target in feas_targets:
            if target not in target_array:
                target_array.append(target)
                info.append(eclipse_info(target.name,time=i,display=False))
    return(target_array, info)


def make_sky_plot(target=None, target_array=None, time=time):
    """
    Print position of target(s) over duration of the night
    
    :param target_array: List or single EclipsingBinary object
    :param time: Astropy.time object
    
    """
    twi_pm = site.twilight_evening_astronomical(time,which='nearest')
    twi_am = site.twilight_morning_astronomical(time,which='next')
    obstimes = Time(twi_pm.iso) + np.linspace(0,9,10)*((Time(twi_am.iso) - Time(twi_pm.iso)) /10)
    plt.figure(figsize=(12,6))
    if target:
        targ = FixedTarget(coord=target.coords, name=target.name)
        plot_sky(targ, site, obstimes)
    elif target_array:
        for j in range(0,len(target_array)):
            targ = FixedTarget(coord=target_array[j].coords, name=target_array[j].name)
            plot_sky(targ, site, obstimes)
            plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()


######################################### MAIN ############################################

def main():
    print("""SCHEDULER FOR XAMIDIMURA \n
          Shortcuts: \n
          t --- generates table of tonight's best targets\n 
          d --- detailed information on a specific target\n
          b --- returns current single best target\n
          q --- quit
    """)
    # INITIALISE
    cmd = None
    recognised_cmds = ['t','d','b','q']
    target_match = []
    for item in target_list:
        target_match.append(item.name)
        
    # MAIN LOOP
    while cmd is not 'q':
        print('Enter new command or press q to quit')
        cmd = input(">>> ")
        if cmd not in recognised_cmds:
            print('Command not recognised, try again')
        
        # TONIGHT'S BEST TARGETS
        if cmd == 't':
            print('*** Tonight\'s best targets --- between astronomical twilights for given date ***')
            print('Which date to use? Format YYYY-MM-DD. Default: current date')
            date = input(">>> ")
            if date == '': time = Time.now()
            else:          time = Time(date) + 12*u.hour
            target_array, info = tonights_best(time)
            # CREATE TABLE OUTPUT
            tab = PrettyTable(['Target name', 'EclType', 'Time start', 'Time end', 'Alt. start', 'Alt. end', 'Priority'])
            for j in range(0,len(target_array)):
                # MAKES EACH ROW OF TABLE FROM INFO LIST
                name = target_array[j].name
                eclty = info[j][0]
                tstart = info[j][1].isot[11:19]
                tend = info[j][2].isot[11:19]
                altstart = info[j][3] 
                altend = info[j][4] 
                priority = info[j][5]
                tab.add_row([name,eclty,tstart,tend,altstart,altend,priority])
            print('Tonight\'s best for {}. All times UTC'.format(str(time)[0:10]))
            print(tab)
            print('')
        
        # DETAILED TARGET INFORMATION
        if cmd == 'd':
            print('*** Detailed target information for specific time ***')
            print('Which time to use? Format YYYY-MM-DD HH:MM:SS (UTC). Default: current time')
            date = input(">>> ")
            if date == '': time = Time.now()
            else:          time = Time(date)
            print('Type target name as it appears in the database.')
            print('Press \'l\' for full target list')
            target = input(">>> ")
            while target not in target_match:
                if target in target_match: break
                elif target == 'l':        print(target_match)
                elif target == 'q':        sys.exit()
                else: print('Invalid target. Try again or press \'l\' for full target list')
                target = input(">>> ")
            result = eclipse_info(target, time)
            if not result:
                print('Target not eclipsing at given time. Returning feasibility report instead')
                why_not(target,time=time)
                print('')
        
        # BEST TARGET TO OBSERVE AT GIVEN TIME
        if cmd == 'b':
            print('*** Best eclipsing target at specific time ***')
            print('Which time to use? Format YYYY-MM-DD HH:MM:SS (UTC). Default: current time')
            date = input(">>> ")
            if date == '': time = Time.now()
            else:          time = Time(date)
            #result = observe_now(time)
            result, status = observe_now(time)
            if result: print('{}\n'.format(result.name))
            else:      print('No targets eclipsing at {}'.format(time.isot))
    
    # QUIT INTERFACE
    if cmd == 'q':
        sys.exit(0)
    

if __name__ == "__main__":
    if interface:
        main()
    else:
        pass