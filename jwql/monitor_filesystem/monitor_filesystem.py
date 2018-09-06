#! /usr/bin/env python

"""
This module is meant to monitor and gather statistics of the filesystem
that hosts data for the ``jwql`` application. This will answer
questions such as the total number of files, how much disk space is
being used, and then plot these values over time.

Authors
-------

    - Misty Cracraft
    - Gray Kanarek

Use
---

    This module can be executed from the command line:

    ::

        python monitor_filesystem.py

    Alternatively, it can be called from scripts with the following
    import statements:

    ::
        from monitor_filesystem import MonitorFilesystem


    Required arguments (in a ``config.json`` file):
    ``outputs`` - The path to the output files needs to be in a
    ``config.json`` file in the ``utils`` directory.


Dependencies
------------

    The user must have a configuration file named ``config.json``
    placed in the ``utils`` directory.

Notes
-----

    The ``MonitorFilesystem.monitor`` function queries the filesystem,
    calculates the statistics and saves the output file(s) in the
    directory specified in the ``config.json`` file. Every time this function
    is called, the associated plots will be updated, and the embeddable
    components will be written out to the appropriate .html and .js files.
"""

from collections import defaultdict
import datetime
import logging
import numpy as np
import os
import subprocess

from jwql.bokeh_templating import BokehTemplate
from jwql.logging.logging_functions import configure_logging, log_info, log_fail
from jwql.permissions.permissions import set_permissions
from jwql.utils.utils import filename_parser
from jwql.utils.utils import get_config

script_dir = os.path.dirname(os.path.abspath(__file__))

fig_formats = """
Figure:
    tools: 'pan,box_zoom,reset,wheel_zoom,save'
    x_axis_type: 'datetime'
    x_axis_label: 'Date'
    sizing_mode: 'stretch_both'
Line:
    line_width: 2
    x: 'dates'
Circle:
    x: 'dates'
Square:
    x: 'dates'
Triangle:
    x: 'dates'
Diamond:
    x: 'dates'
Asterisk:
    x: 'dates'
XGlyph:
    x: 'dates'
"""

class MonitorFilesystem(BokehTemplate):
    
    def pre_init(self):
        self._embed = True
        
        #App design
        self.format_string = fig_formats
        self.interface_file = os.path.join(script_dir, "monitor_filesystem_interface.yaml")
        
        # Get path, directories and files in system and count files in all directories
        self.settings = get_config()
        self.filesystem = self.settings['filesystem']
        self.outputs_dir = os.path.join(self.settings['outputs'], 'monitor_filesystem')
        
        self.statistics = defaultdict(list)
        self.results = defaultdict(list)
        self.sizes = defaultdict(list)
        
        self.types_k = ['circle', 'diamond', 'square', 'triangle', 
                        'asterisk'] + ['x']*6
        self.types_y = ['fits', 'uncal', 'cal', 'rate', 'rateint', 
                        'i2d', 'nrc', 'nrs', 'nis', 'mir', 'fgs']
        self.types_c = ['black', 'red', 'blue', 'green', 'orange', 'purple', 
                        'midnightblue', 'springgreen', 'darkcyan', 
                        'dodgerblue', 'darkred']
        self.types_l = ['Total FITS files', 'Uncalibrated FITS files',
                        'Calibrated FITS files', 'Rate FITS files',
                        'Rateints FITS files', 'I2D FITS files',
                        'NIRCam FITS files', 'NIRSpec FITS files',
                        'NIRISS FITS files', 'MIRI FITS files',
                        'FGS FITS files']
        
    post_init = None
    
    @log_fail
    @log_info
    def monitor(self):
        """
        Monitoring script to inventory the JWST filesystem, save file 
        statistics, and generate plots.
        """
    
        # Begin logging
        logging.info('Beginning filesystem monitoring.')
    
        # re-initialize dictionaries for output
        results_dict = defaultdict(int)
        size_dict = defaultdict(float)
        
        # Walk through all directories recursively and count files
        logging.info('Searching filesystem...')
        for dirpath, dirs, files in os.walk(self.filesystem):
            results_dict['file_count'] += len(files)  # find number of all files
            for filename in files:
                file_path = os.path.join(dirpath, filename)
                if filename.endswith(".fits"):  # find total number of fits files
                    results_dict['fits_files'] += 1
                    size_dict['size_fits'] += os.path.getsize(file_path)
                    suffix = filename_parser(filename)['suffix']
                    results_dict[suffix] += 1
                    size_dict[suffix] += os.path.getsize(file_path)
                    detector = filename_parser(filename)['detector']
                    instrument = detector[0:3]  # first three characters of detector specify instrument
                    results_dict[instrument] += 1
                    size_dict[instrument] += os.path.getsize(file_path)
        logging.info('{} files found in filesystem'.format(results_dict['fits_files']))
    
        # Get df style stats on file system
        out = subprocess.check_output('df {}'.format(self.filesystem), shell=True)
        outstring = out.decode("utf-8")  # put into string for parsing from byte format
        parsed = outstring.split(sep=None)
    
        # Select desired elements from parsed string
        stats = {
                'total': int(parsed[8]),  # in blocks of 512 bytes
                'used': int(parsed[9]),
                'available': int(parsed[10]),
                'percent_used': parsed[11],
                'file_count': results_dict.pop('file_count'),
                'timestamp': datetime.datetime.now().isoformat(sep='T', timespec='auto')  # get date of stats
                }
        
        #store results & sizes in the appropriate dictionaries
        for key, val in results_dict.items():
            self.results[key].append(val)
        for key, val in size_dict.items():
            self.sizes[key].append(val)
        for key, val in stats.items():
            self.statistics[key].append(val)
    
        # set up output file and write stats
        statsfile = os.path.join(self.outputs_dir, 'statsfile.txt')
        with open(statsfile, "a+") as f:
            f.write("{timestamp} {file_count:15d} {total:15d} {available:15d} {used:15d} {percent_used}\n".format(**stats))
        set_permissions(statsfile)
        logging.info('Saved file statistics to: {}'.format(statsfile))
    
        output_stub = "{fits_files} {uncal} {cal} {rate} {rateints} {i2d} {nrc} {nrs} {nis} {mir} {gui}\n"
        # set up and read out stats on files by type
        filesbytype = os.path.join(self.outputs_dir, 'filesbytype.txt')
        with open(filesbytype, "a+") as f2:
            f2.write(output_stub.format(**results_dict))
        set_permissions(filesbytype, verbose=False)
        logging.info('Saved file statistics by type to {}'.format(filesbytype))
    
        # set up file size by type file
        sizebytype = os.path.join(self.outputs_dir, 'sizebytype.txt')
        with open(sizebytype, "a+") as f3:
            f3.write(output_stub.format(**size_dict))
        set_permissions(sizebytype, verbose=False)
        logging.info('Saved file sizes by type to {}'.format(sizebytype))
    
        logging.info('Filesystem statistics calculation complete.')
        
        #Update the plots based on new information
        self.update_plots()
        
    def update_plots(self):
        """
        Update the ColumnDataSource objects for the filesystem monitor plots.
        """
        
        logging.info("Beginning filesystem statistics monitor plot updates")
        
        # We'll ensure that all the statistics are in the correct formats
        dates = np.array(self.statistics['timestamp'], dtype='datetime64')
        
        self.refs['source_filecount'].data = {
                'dates': dates,
                'filecount': np.array(self.statistics['file_count'], dtype=float)
                }
        
        self.refs['source_stats'].data = {
                'dates': dates,
                'systemsize': np.array(self.statistics['total'], dtype=float) / (1024.**3),
                'freesize': np.array(self.statistics['available'], dtype=float) / (1024.**3),
                'usedsize': np.array(self.statistics['used'], dtype=float) / (1024.**3)
                }
        
        self.refs['source_files'].data = {
                'dates': dates,
                'fits': np.array(self.results['fits_files'], dtype=int),
                'uncal': np.array(self.results['uncal'], dtype=int),
                'cal': np.array(self.results['cal'], dtype=int),
                'rate': np.array(self.results['rate'], dtype=int),
                'rateint': np.array(self.results['rateint'], dtype=int),
                'i2d': np.array(self.results['i2d'], dtype=int),
                'nrc': np.array(self.results['nrc'], dtype=int),
                'nrs': np.array(self.results['nrs'], dtype=int),
                'nis': np.array(self.results['nis'], dtype=int),
                'mir': np.array(self.results['mir'], dtype=int),
                'fgs': np.array(self.results['gui'], dtype=int)
            }
        
        self.refs['source_sizes'].data = {
                'dates': dates,
                'fits': np.array(self.sizes['fits_files'], dtype=float) / (1024.**3),
                'uncal': np.array(self.sizes['uncal'], dtype=float) / (1024.**3),
                'cal': np.array(self.sizes['cal'], dtype=float) / (1024.**3),
                'rate': np.array(self.sizes['rate'], dtype=float) / (1024.**3),
                'rateint': np.array(self.sizes['rateint'], dtype=float) / (1024.**3),
                'i2d': np.array(self.sizes['i2d'], dtype=float) / (1024.**3),
                'nrc': np.array(self.sizes['nrc'], dtype=float) / (1024.**3),
                'nrs': np.array(self.sizes['nrs'], dtype=float) / (1024.**3),
                'nis': np.array(self.sizes['nis'], dtype=float) / (1024.**3),
                'mir': np.array(self.sizes['mir'], dtype=float) / (1024.**3),
                'fgs': np.array(self.sizes['gui'], dtype=float) / (1024.**3)
            }
        
        # Write scripts out to files
        for name in ['filecount', 'system_stats', 'filecount_type', 'size_type']:
            script, div = self.embed('fig_'+name)
            div_outfile = os.path.join(self.outputs_dir, "{}_component.html".format(name))
            with open(div_outfile, 'w') as f:
                f.write(div)
                f.close()
            set_permissions(div_outfile)
    
            script_outfile = os.path.join(self.outputs_dir, "{}_component.js".format(name))
            with open(script_outfile, 'w') as f:
                f.write(script)
                f.close()
            set_permissions(script_outfile)
            
            logging.info('Saved components files: {}_component.html and {}_component.js'.format(name, name))

        logging.info('Filesystem statistics plot updates complete.')
        


if __name__ == '__main__':

    # Configure logging
    module = os.path.basename(__file__).strip('.py')
    configure_logging(module)

    monitor = MonitorFilesystem()
    monitor.monitor()