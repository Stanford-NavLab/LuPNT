import argparse

# import gnss_lib_py as glp
from datetime import datetime, timedelta
import numpy as np

# def parsing_routine(start_time, end_time, dt, constellations, sv_ids, verb = False):
#     """Retrieve high precision ephemeris of Gnss constellations with gnss-lib-py (https://gnss-lib-py.readthedocs.io/en/latest/index.html)

#     Args:
#         start_time (datetime.datetime): start time of parsing in UTC
#         end_time (datetime.datetime): end time of parsing in UTC
#         dt (float): time interval (seconds)
#         constellations (list): Gnss constellations to parse
#         sv_ids (list): list of SV ids to parse
#         verb (bool, optional): enables verbose. Defaults to False.

#     Returns:
#         NavData(): obj with precise ephemeris
#     """
#     #generate all the time instances we care about
#     t0 = start_time
#     times = []
#     gnss_ids = []
#     sv_ids_list = []

#     #the total length of each list should be len(time list)*len(constellations)*len(sv_ids)
#     while t0 <= end_time :
#         for constellation in constellations:
#             for sv_id in sv_ids:
#                 gnss_ids.append(constellation)
#                 sv_ids_list.append(sv_id)
#                 times.append(t0)
#         #next time step
#         t0 = t0 + timedelta(seconds=dt)

#     #must convert to gps milliseconds (ASSUMING THE DATETIMES ARE IN UTC)
#     gps_millis = glp.datetime_to_gps_millis(times)
#     #create NaVData() object
#     lupnt_sp3 = glp.NavData()
#     lupnt_sp3['gps_millis'] = gps_millis
#     lupnt_sp3['gnss_id'] = np.asarray(gnss_ids, dtype=object)
#     lupnt_sp3['sv_id'] = sv_ids_list
#     lupnt_sp3['raw_pr_m'] = 0

#     #get all the data and return it to save it as csv
#     lupnt_sp3_w_sv_states = glp.add_sv_states_precise(lupnt_sp3, download_directory = 'data/ephemeris/gnsslibpy', verbose = verb)

#     return lupnt_sp3_w_sv_states

# def main():
#     parser = argparse.ArgumentParser(
#                     prog='Ephemeris Parsing ft. gnss lib py',
#                     description='Gathers necessary ephemeris files of Gnss constellations in a CSV file for LuPNT usage')
#     #add arguments to parse
#     parser.add_argument('start', type= datetime.fromisoformat, help ='ISOformat - YYYY-MM-DD:HH:mm:ss')
#     parser.add_argument('end', type= datetime.fromisoformat, help = 'ISOformat - YYYY-MM-DD:HH:mm:ss')
#     #if not specified, these arguments have defaults
#     parser.add_argument('--dt', type= float, default= 60*5, help = 'Time interval (default is 5 minutes)')
#     parser.add_argument('--constellation', '--c', choices = ['gps', 'galileo', 'glonass', 'beidou', 'qzss', 'irnss', 'sbas'], nargs="+", default=['gps'])
#     parser.add_argument('--svids', '--s', nargs = "+", default= range(1,33), help = 'If no ids specified, we are grabbing within range(1,33)')
#     parser.add_argument('--outfile', '--o', type = str, default= 'data/gnss_ephemeris.csv')
#     #parse the arguments provided (includes defaults)
#     args = parser.parse_args()
#     start_time = args.start
#     end_time = args.end
#     dt = args.dt
#     constellations = args.constellation
#     sv_ids = args.svids
#     outfile =  args.outfile
#     #run gnss lib py functions
#     lupnt_sp3 = parsing_routine(start_time, end_time, dt, constellations, sv_ids)
#     #get csv file
#     lupnt_sp3.to_csv(outfile)

# if __name__ == '__main__':
#     main()
