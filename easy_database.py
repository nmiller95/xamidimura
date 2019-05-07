"""
Functions to manipulate target info table and priority table in database.

1. Run file
2. Consider using the following lines before using other functions:
    dbconn, dbcurs = connect_database.connect_to()
    t = targets_to_dataframe(dbconn)
    s = priorities_to_dataframe(dbconn)

"""

############################################# IMPORTS ##############################################

import connect_database
import settings_and_error_codes as set_err_codes
from astropy.table import Table


############################# DISPLAY AS DATAFRAMES #############################

def targets_to_dataframe(conn):
    """
    Displays target info table as pandas dataframe
    
    :param curs: Database cursor
    
    :returns: Dataframe
    
    """
    return connect_database.get_table_into_pandas('target_info',conn)


def priorities_to_dataframe(conn):
    """
    Displays priority table as pandas dataframe
    
    :param curs: Database cursor
    
    :returns: Dataframe
    
    """
    return connect_database.get_table_into_pandas('priority_table',conn)


############################# TARGET INFO TABLE #############################

def wipe_target_table(curs, conn):
    """
    Wipes target info table
    
    :param curs: Database cursor
    :param conn: Database connection
    
    :returns: Empty target info table in database
    
    """
    connect_database.remove_table_if_exists(curs, set_err_codes.TARGET_INFORMATION_TABLE)
    curs.execute('CREATE TABLE '+set_err_codes.TARGET_INFORMATION_TABLE+\
                   ' (TAR_ID TEXT, TAR_NAME TEXT NOT NULL, RA TEXT NOT NULL, \
                   DEC TEXT NOT NULL, TAR_TYPE TEXT NOT NULL, MAGNITUDE REAL, \
                   SPEC_TYPE TEXT, T_0 REAL, PERIOD REAL, DEPTH1 REAL, \
                   DEPTH2 REAL, WIDTH1 REAL, WIDTH2 REAL, PHASE2 REAL, NOTES TEXT, \
                   MAG_SOURCE TEXT, PRIMARY KEY (TAR_ID ASC));')
    conn.commit()


def filter_targets(file_name="database/master_targets.csv"):
    """
    Reads and filters new targets from file
    
    :param file_name: String filename of csv file containing targets
    
    :returns: List of rows ready to add to database
    
    """
    tab = Table.read(file_name,format="csv")
    tab = tab.filled(-99999.)
    new_data = []
    for i in tab:
        mag    = i['MAGNITUDE'] >= 6. and i['MAGNITUDE'] <= 13.
        width  = i['WIDTH1'] <= 0.05 and i['WIDTH2'] <= 0.05
        period = i['PERIOD'] >= 5. and i['PERIOD'] <= 100.
        depth  = i['DEPTH1'] >= 0.1 and i['DEPTH2'] >= 0.05
        dec    = int(i['DEC'][0:3]) < 30 
        if mag and width and period and depth and dec:
            new_data.append(list(i))
    print("Targets filtered from original {} to {}".format(len(tab),len(new_data)))
    return new_data


def fill_target_table(new_data, curs, conn, overwrite=False):
    """
    Fills target info table with new data
    
    :param new_data: Table or list or list of lists
    :param curs: Database cursor
    :param conn: Database connection
    
    :returns: Filled target info table in database
    
    """
    for i in new_data:
        connect_database.add_target_to_database(list(i), curs, conn, overwrite_exsiting = overwrite)
    conn.commit()


############################# PRIORITY TABLE #############################

def wipe_priority_table(curs, conn):
    """
    Wipes priority table
    
    :param curs: Database cursor
    :param conn: Database connection
    
    :returns: Empty priority table in database
    
    """
    connect_database.remove_table_if_exists(curs, set_err_codes.PRIORITY_TABLE)
    curs.execute('CREATE TABLE '+set_err_codes.PRIORITY_TABLE+' (PRIORITY_ID \
                 INTEGER, TAR_ID TEXT UNIQUE NOT NULL, NUMBER_BINS_1 INTEGER, \
                 NUMBER_OBS_PER_BIN_1 INTEGER, COMPLETENESS_1 TEXT, \
                 NUMBER_BINS_2 INTEGER, NUMBER_OBS_PER_BIN_2 INTEGER, \
                 COMPLETENESS_2 TEXT,URGENCY INTEGER, PRIMARY KEY (PRIORITY_ID), \
                 FOREIGN KEY(TAR_ID) REFERENCES '+set_err_codes.TARGET_INFORMATION_TABLE+'(TAR_ID));')
    conn.commit()


def fill_priority_table(target_table, curs, conn, no_bins=5, no_obs=4):
    """
    Fills priority table with defaults
    
    :param target_table: Target info table in pandas dataframe format
    :param no_bins: Integer, number of bins per eclipse per target
    :param no_obs: Integer, number of observations per bin to require
    :param curs: Database cursor
    :param conn: Database connection
    
    :returns: Filled priority table in database
    
    """
    tar_ids = list(target_table['TAR_ID'])
    no_bins_1 = [no_bins]*len(target_table)
    noobs_bin_1 = [no_obs]*len(target_table)
    complete_1 = ['00'*no_bins]*len(target_table)
    no_bins_2 = [no_bins]*len(target_table)
    noobs_bin_2 = [no_obs]*len(target_table)
    complete_2 = ['00'*no_bins]*len(target_table)
    urgency = [0]*len(target_table)
    for i in range(0,len(target_table)):
        item = [tar_ids[i], no_bins_1[i], noobs_bin_1[i], complete_1[i], 
                no_bins_2[i], noobs_bin_2[i], complete_2[i], urgency[i]]
        connect_database.add_priority_to_database(item, curs, conn)


def change_target_urgency(target_name, value, curs, conn):
    """
    Changes urgency of single target
    
    :param target_name: String, target name
    :param value: Integer, value to set urgency to
    :param curs: Database cursor
    :param conn: Database connection
    
    :returns: Updated priority table in database
    
    """
    urgs = [0,1,2,3,4]
    if value not in urgs:
        raise ValueError('Not a valid urgency score')
    else:
        target_info = connect_database.match_target_name(target_name,'target_info',curs)
        if target_info == []:
            raise ValueError('Target not in database')
        else:
            target_id = target_info[0][0]
            connect_database.update_priority_info(target_id, 'URGENCY', value, curs, conn)