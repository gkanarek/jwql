"""This module is home to a suite of MAST queries that gather bulk properties
of available JWST data for JWQL

Authors
-------

    Joe Filippazzo
    Gray Kanarek

Use
---

    To get an inventory of all JWST files do:
    ::

        from jwql.monitor_mast import monitor_mast
        inventory, keywords = monitor_mast.jwst_inventory()
"""

import logging
import os

from astroquery.mast import Mast
import pandas as pd

from jwql.logging.logging_functions import configure_logging, log_info, log_fail
from jwql.permissions.permissions import set_permissions
from jwql.utils.utils import get_config, JWST_DATAPRODUCTS, JWST_INSTRUMENTS
from jwql.bokeh_templating import BokehTemplate

script_dir = os.path.dirname(os.path.abspath(__file__))

fig_formats = """
Donut:
    label: ['instrument', 'dataproduct']
    values: 'files'
    text_font_size: '12pt'
    hover_text: 'files'
    name: "JWST Inventory"
    plot_width: 600
    plot_height: 600
"""
    

def instrument_inventory(instrument, dataproduct=JWST_DATAPRODUCTS,
                         add_filters=None, add_requests=None,
                         caom=False, return_data=False):
    """Get the counts for a given instrument and data product

    Parameters
    ----------
    instrument: str
        The instrument name, i.e. ['NIRISS','NIRCam','NIRSpec','MIRI','FGS']
    dataproduct: sequence, str
        The type of data product to search
    add_filters: dict
        The ('paramName':'values') pairs to include in the 'filters' argument
        of the request e.g. add_filters = {'filter':'GR150R'}
    add_requests: dict
        The ('request':'value') pairs to include in the request
        e.g. add_requests = {'pagesize':1, 'page':1}
    caom: bool
        Query CAOM service
    return_data: bool
        Return the actual data instead of counts only

    Returns
    -------
    int, dict
        The number of database records that satisfy the search criteria
        or a dictionary of the data if `return_data=True`
    """
    filters = []

    # Make sure the dataproduct is a list
    if isinstance(dataproduct, str):
        dataproduct = [dataproduct]

    # Make sure the instrument is supported
    if instrument.lower() not in [ins.lower() for ins in JWST_INSTRUMENTS]:
        raise TypeError('Supported instruments include:', JWST_INSTRUMENTS)

    # CAOM service
    if caom:

        # Declare the service
        service = 'Mast.Caom.Filtered'

        # Set the filters
        filters += [{'paramName': 'obs_collection', 'values': ['JWST']},
                    {'paramName': 'instrument_name', 'values': [instrument]},
                    {'paramName': 'dataproduct_type', 'values': dataproduct}]

    # Instruent filtered service
    else:

        # Declare the service
        service = 'Mast.Jwst.Filtered.{}'.format(instrument.title())

    # Include additonal filters
    if isinstance(add_filters, dict):
        filters += [{"paramName": name, "values": [val]}
                    for name, val in add_filters.items()]

    # Assemble the request
    params = {'columns': 'COUNT_BIG(*)',
              'filters': filters,
              'removenullcolumns': True}

    # Just get the counts
    if return_data:
        params['columns'] = '*'

    # Add requests
    if isinstance(add_requests, dict):
        params.update(add_requests)

    response = Mast.service_request_async(service, params)
    result = response[0].json()

    # Return all the data
    if return_data:
        return result

    # Or just the counts
    else:
        return result['data'][0]['Column1']


def instrument_keywords(instrument, caom=False):
    """Get the keywords for a given instrument service

    Parameters
    ----------
    instrument: str
        The instrument name, i.e. ['NIRISS','NIRCam','NIRSpec','MIRI','FGS']
    caom: bool
        Query CAOM service

    Returns
    -------
    pd.DataFrame
        A DataFrame of the keywords
    """
    # Retrieve one dataset to get header keywords
    sample = instrument_inventory(instrument, return_data=True, caom=caom,
                                  add_requests={'pagesize': 1, 'page': 1})
    data = [[i['name'], i['type']] for i in sample['fields']]
    keywords = pd.DataFrame(data, columns=('keyword', 'dtype'))

    return keywords


def jwst_inventory(instruments=JWST_INSTRUMENTS,
                   dataproducts=['image', 'spectrum', 'cube'],
                   caom=False):
    """Gather a full inventory of all JWST data in each instrument
    service by instrument/dtype

    Parameters
    ----------
    instruments: sequence
        The list of instruments to count
    dataproducts: sequence
        The types of dataproducts to count
    caom: bool
        Query CAOM service
    plot: bool
        Return a pie chart of the data

    Returns
    -------
    astropy.table.table.Table
        The table of record counts for each instrument and mode
    """
    logging.info('Searching database...')
    # Iterate through instruments
    inventory, keywords = [], {}
    for instrument in instruments:
        ins = [instrument]
        for dp in dataproducts:
            count = instrument_inventory(instrument, dataproduct=dp, caom=caom)
            ins.append(count)

        # Get the total
        ins.append(sum(ins[-3:]))

        # Add it to the list
        inventory.append(ins)

        # Add the keywords to the dict
        keywords[instrument] = instrument_keywords(instrument, caom=caom)

    logging.info('Completed database search for {} instruments and {} data products.'.
                 format(instruments, dataproducts))

    # Make the table
    all_cols = ['instrument']+dataproducts+['total']
    table = pd.DataFrame(inventory, columns=all_cols)

    # Melt the table
    table = pd.melt(table, id_vars=['instrument'],
                    value_vars=dataproducts,
                    value_name='files', var_name='dataproduct')

    return table, keywords

class MastMonitor(BokehTemplate):
    
    instruments = JWST_INSTRUMENTS
    dataproducts = ['image', 'spectrum', 'cube']
    
    def pre_init(self):
        self._embed = True
        
        #App design
        self.format_string = None
        self.interface_file = os.path.join(script_dir, "monitor_mast_interface.yaml")
        
        self.settings = get_config()
        self.output_dir = self.settings['outputs']
        
        #Have to initialize DataFrames for plotting
        self.monitor(initialize=True)
        
        
    post_init = None
    
    @log_info
    @log_fail
    def monitor(self, initialize=False):
        """Tabulates the inventory of all JWST data products in the MAST
        archive and generates plots.
        """
        logging.info('Beginning database monitoring.')
        self.jwst_table, self.jwst_keywords = jwst_inventory(self.instruments,
                                                             self.dataproducts,
                                                             caom=False)
        self.caom_table, self.caom_keywords = jwst_inventory(self.instruments,
                                                             self.dataproducts,
                                                             caom=True)
        if not initialize:
            self.update_plots()
    
    def update_plots(self):
        # Determine plot location and names
        html_stub = 'database_monitor_{}_component.html'
        js_stub = 'database_monitor_{}_component.js'
        
        self.refs["fig_jwst"].data = self.jwst_table
        self.refs["fig_caom"].data = self.caom_table

        # Save the plots as components
        for name in ["jwst", "caom"]:
            script, div = self.embed("fig_"+name)
            html_name = html_stub.format(name)
            js_name = js_stub.format(name)
            div_outfile = os.path.join(self.output_dir, 'monitor_mast', 
                                       html_name)
            with open(div_outfile, 'w') as f:
                f.write(div)
                f.close()
            set_permissions(div_outfile)

            script_outfile = os.path.join(self.output_dir, 'monitor_mast', 
                                          js_name)
            with open(script_outfile, 'w') as f:
                f.write(script)
                f.close()
            set_permissions(script_outfile)

            logging.info('Saved Bokeh components files: {} and {}'.format(html_name, 
                                                                          js_name))

if __name__ == '__main__':

    # Configure logging
    module = os.path.basename(__file__).strip('.py')
    configure_logging(module)

    # Run the monitors
    monitor = MastMonitor()
    monitor.monitor()
