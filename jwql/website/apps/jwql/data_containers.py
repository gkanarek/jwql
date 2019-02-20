"""Various functions to collect data to be used by the ``views`` of the
``jwql`` app.

This module contains several functions that assist in collecting and
producing various data to be rendered in ``views.py`` for use by the
``jwql`` app.

Authors
-------

    - Lauren Chambers
    - Matthew Bourque

Use
---

    The functions within this module are intended to be imported and
    used by ``views.py``, e.g.:

    ::
        from .data_containers import get_proposal_info
"""

import glob
import os

from astropy.io import fits
from astropy.time import Time
from astroquery.mast import Mast
from bokeh.embed import components
from bokeh.plotting import figure
import numpy as np

from .forms import MnemonicSearchForm, MnemonicQueryForm
from jwql.jwql_monitors import monitor_cron_jobs
from jwql.utils.constants import MONITORS
from jwql.utils.preview_image import PreviewImage
from jwql.utils.utils import get_config, filename_parser
from jwql.utils.engineering_database import query_mnemonic_info, query_single_mnemonic

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
FILESYSTEM_DIR = os.path.join(get_config()['jwql_dir'], 'filesystem')
PREVIEW_IMAGE_FILESYSTEM = os.path.join(get_config()['jwql_dir'], 'preview_images')
THUMBNAIL_FILESYSTEM = os.path.join(get_config()['jwql_dir'], 'thumbnails')
PACKAGE_DIR = os.path.dirname(__location__.split('website')[0])
REPO_DIR = os.path.split(PACKAGE_DIR)[0]


def webpage_template_data():
    """An example data container function for the webpage template.

    Define variables to pass to the webpage_template view, including
    the output of a Bokeh plot to embed.

    Returns
    -------
    jwst_launch_date : int
        The current expected JWST launch date
    plot_data : list
        A list containing the JavaScript and HTML content for the
        Bokeh plot
    """
    # Define a single variable to pass
    jwst_launch_date = 2021

    # Define a basic bokeh plot showing value as a function of time.
    # Start by defining the data for the plot.
    list_of_dates = ['2018-02-11 00:00:00.000', '2018-02-12 00:00:00.000',
                     '2018-02-13 00:00:00.000', '2018-02-14 00:00:00.000',
                     '2018-02-15 00:00:00.000', '2018-02-16 00:00:00.000',
                     '2018-02-17 00:00:00.000']
    date = Time(list_of_dates, format='iso').datetime
    avg_temp = [33, 36, 42, 43, 56, 44, 34]

    # Build the plot with Bokeh
    p1 = figure(tools='pan,box_zoom,reset,wheel_zoom,save', x_axis_type='datetime',
                title='Baltimore Winter Temperatures', x_axis_label='Date',
                y_axis_label='Average Temperature (Degrees F)')
    p1.line(date, avg_temp, line_width=1, line_color='blue', line_dash='dashed')
    p1.circle(date, avg_temp, color='blue')

    # Save out the JavaScript and HTML
    script, div = components(p1)
    plot_data = [div, script]

    return jwst_launch_date, plot_data

def data_trending():
    """Container for Miri datatrending dashboard and components

    Returns
    -------
    variables : int
        nonsense
    dashboard : list
        A list containing the JavaScript and HTML content for the
        dashboard
    """

    import jwql.instrument_monitors.miri_monitors.data_trending.dashboard as dash

    dashboard = dash.data_trending_dashboard()

    variables = {
        "var1": 10,
        "var2": 20
    }

    return variables, dashboard


def get_acknowledgements():
    """Returns a list of individuals who are acknowledged on the
    ``about`` page.

    The list is generated by reading in the contents of the ``jwql``
    ``README`` file.  In this way, the website will automatically
    update with updates to the ``README`` file.

    Returns
    -------
    acknowledgements : list
        A list of individuals to be acknowledged.
    """

    # Locate README file
    readme_file = os.path.join(REPO_DIR, 'README.md')

    # Get contents of the README file
    with open(readme_file, 'r') as f:
        data = f.readlines()

    # Find where the acknowledgements start
    for i, line in enumerate(data):
        if 'Acknowledgments' in line:
            index = i

    # Parse out the list of individuals
    acknowledgements = data[index + 1:]
    acknowledgements = [item.strip().replace('- ', '').split(' [@')[0].strip() for item in acknowledgements]

    return acknowledgements


def get_all_proposals():
    """Return a list of all proposals that exist in the filesystem.

    Returns
    -------
    proposals : list
        A list of proposal numbers for all proposals that exist in the
        filesystem
    """

    proposals = glob.glob(os.path.join(FILESYSTEM_DIR, '*'))
    proposals = [proposal.split('jw')[-1] for proposal in proposals]
    proposals = [proposal for proposal in proposals if len(proposal) == 5]

    return proposals


def get_dashboard_components():
    """Build and return dictionaries containing components and html
    needed for the dashboard.

    Returns
    -------
    dashboard_components : dict
        A dictionary containing components needed for the dashboard.
    dashboard_html : dict
        A dictionary containing full HTML needed for the dashboard.
    """

    output_dir = get_config()['outputs']
    name_dict = {'': '',
                 'monitor_mast': 'Database Monitor',
                 'database_monitor_jwst': 'JWST',
                 'database_monitor_caom': 'JWST (CAOM)',
                 'monitor_filesystem': 'Filesystem Monitor',
                 'filecount_type': 'Total File Counts by Type',
                 'size_type': 'Total File Sizes by Type',
                 'filecount': 'Total File Counts',
                 'system_stats': 'System Statistics'}

    # Exclude monitors that can't be saved as components
    exclude_list = ['monitor_cron_jobs']

    # Run the cron job monitor to produce an updated table
    monitor_cron_jobs.status(production_mode=True)

    # Build dictionary of components
    dashboard_components = {}
    for dir_name, subdir_list, file_list in os.walk(output_dir):
        monitor_name = os.path.basename(dir_name)
        if monitor_name not in exclude_list:
            dashboard_components[name_dict[monitor_name]] = {}
            for fname in file_list:
                if 'component' in fname:
                    full_fname = '{}/{}'.format(monitor_name, fname)
                    plot_name = fname.split('_component')[0]

                    # Get the div
                    html_file = full_fname.split('.')[0] + '.html'
                    with open(os.path.join(output_dir, html_file)) as f:
                        div = f.read()

                    # Get the script
                    js_file = full_fname.split('.')[0] + '.js'
                    with open(os.path.join(output_dir, js_file)) as f:
                        script = f.read()
                    dashboard_components[name_dict[monitor_name]][name_dict[plot_name]] = [div, script]

    # Add HTML that cannot be saved as components to the dictionary
    with open(os.path.join(output_dir, 'monitor_cron_jobs', 'cron_status_table.html')) as f:
        cron_status_table_html = f.read()
    dashboard_html = {}
    dashboard_html['Cron Job Monitor'] = cron_status_table_html

    return dashboard_components, dashboard_html


def get_edb_components(request):
    """Return dictionary with content needed for the EDB page.

    Parameters
    ----------
    request : HttpRequest object
        Incoming request from the webpage

    Returns
    -------
    edb_components : dict
        Dictionary with the required components

    """
    mnemonic_name_search_result = {}
    mnemonic_query_result = {}
    mnemonic_query_result_plot = None

    # If this is a POST request, we need to process the form data
    if request.method == 'POST':

        if 'mnemonic_name_search' in request.POST.keys():
            mnemonic_name_search_form = MnemonicSearchForm(request.POST,
                                                           prefix='mnemonic_name_search')

            if mnemonic_name_search_form.is_valid():
                mnemonic_identifier = mnemonic_name_search_form['search'].value()
                if mnemonic_identifier is not None:
                    mnemonic_name_search_result = query_mnemonic_info(mnemonic_identifier)

            # create forms for search fields not clicked
            mnemonic_query_form = MnemonicQueryForm(prefix='mnemonic_query')

        elif 'mnemonic_query' in request.POST.keys():
            mnemonic_query_form = MnemonicQueryForm(request.POST, prefix='mnemonic_query')

            # proceed only if entries make sense
            if mnemonic_query_form.is_valid():
                mnemonic_identifier = mnemonic_query_form['search'].value()
                start_time = Time(mnemonic_query_form['start_time'].value(), format='iso')
                end_time = Time(mnemonic_query_form['end_time'].value(), format='iso')

                if mnemonic_identifier is not None:
                    mnemonic_query_result = query_single_mnemonic(mnemonic_identifier, start_time,
                                                                  end_time)
                    mnemonic_query_result_plot = mnemonic_query_result.bokeh_plot()

            # create forms for search fields not clicked
            mnemonic_name_search_form = MnemonicSearchForm(prefix='mnemonic_name_search')

    else:
        mnemonic_name_search_form = MnemonicSearchForm(prefix='mnemonic_name_search')
        mnemonic_query_form = MnemonicQueryForm(prefix='mnemonic_query')

    edb_components = {'mnemonic_query_form' : mnemonic_query_form,
                      'mnemonic_query_result' : mnemonic_query_result,
                      'mnemonic_query_result_plot' : mnemonic_query_result_plot,
                      'mnemonic_name_search_form' : mnemonic_name_search_form,
                      'mnemonic_name_search_result': mnemonic_name_search_result}

    return edb_components


def get_expstart(rootname):
    """Return the exposure start time (``expstart``) for the given
    group of files.

    The ``expstart`` is gathered from a query to the
    ``astroquery.mast`` service.

    Parameters
    ----------
    rootname : str
        The rootname of the observation of interest (e.g.
        ``jw86700006001_02101_00006_guider1``).

    Returns
    -------
    expstart : float
        The exposure start time of the observation (in MJD).
    """

    return 5000.00


def get_filenames_by_instrument(instrument):
    """Returns a list of paths to files that match the given
    ``instrument``.

    Parameters
    ----------
    instrument : str
        The instrument of interest (e.g. `FGS`).

    Returns
    -------
    filepaths : list
        A list of full paths to the files that match the given
        instrument.
    """

    # Query files from MAST database
    # filepaths, filenames = DatabaseConnection('MAST', instrument=instrument).\
    #     get_files_for_instrument(instrument)

    # Find all of the matching files in filesytem
    # (TEMPORARY WHILE THE MAST STUFF IS BEING WORKED OUT)
    instrument_match = {'FGS': 'guider',
                        'MIRI': 'mir',
                        'NIRCam': 'nrc',
                        'NIRISS': 'nis',
                        'NIRSpec': 'nrs'}
    search_filepath = os.path.join(FILESYSTEM_DIR, '*', '*.fits')
    filepaths = [f for f in glob.glob(search_filepath) if instrument_match[instrument] in f]

    return filepaths


def get_filenames_by_proposal(proposal):
    """Return a list of filenames that are available in the filesystem
    for the given ``proposal``.

    Parameters
    ----------
    proposal : str
        The five-digit proposal number (e.g. ``88600``).

    Returns
    -------
    filenames : list
        A list of filenames associated with the given ``proposal``.
    """

    filenames = sorted(glob.glob(os.path.join(
        FILESYSTEM_DIR, 'jw{}'.format(proposal), '*')))
    filenames = [os.path.basename(filename) for filename in filenames]

    return filenames


def get_filenames_by_rootname(rootname):
    """Return a list of filenames available in the filesystem that
    are part of the given ``rootname``.

    Parameters
    ----------
    rootname : str
        The rootname of interest (e.g. ``jw86600008001_02101_00007_guider2``).

    Returns
    -------
    filenames : list
        A list of filenames associated with the given ``rootname``.
    """

    proposal = rootname.split('_')[0].split('jw')[-1][0:5]
    filenames = sorted(glob.glob(os.path.join(
        FILESYSTEM_DIR,
        'jw{}'.format(proposal),
        '{}*'.format(rootname))))
    filenames = [os.path.basename(filename) for filename in filenames]

    return filenames


def get_header_info(file):
    """Return the header information for a given ``file``.

    Parameters
    ----------
    file : str
        The name of the file of interest.

    Returns
    -------
    header : str
        The primary FITS header for the given ``file``.
    """

    dirname = file[:7]
    fits_filepath = os.path.join(FILESYSTEM_DIR, dirname, file)
    header = fits.getheader(fits_filepath, ext=0).tostring(sep='\n')

    return header


def get_image_info(file_root, rewrite):
    """Build and return a dictionary containing information for a given
    ``file_root``.

    Parameters
    ----------
    file_root : str
        The rootname of the file of interest.
    rewrite : bool
        ``True`` if the corresponding JPEG needs to be rewritten,
        ``False`` if not.

    Returns
    -------
    image_info : dict
        A dictionary containing various information for the given
        ``file_root``.
    """

    # Initialize dictionary to store information
    image_info = {}
    image_info['all_jpegs'] = []
    image_info['suffixes'] = []
    image_info['num_ints'] = {}

    preview_dir = os.path.join(get_config()['jwql_dir'], 'preview_images')

    # Find all of the matching files
    dirname = file_root[:7]
    search_filepath = os.path.join(FILESYSTEM_DIR, dirname, file_root + '*.fits')
    image_info['all_files'] = glob.glob(search_filepath)

    for file in image_info['all_files']:

        # Get suffix information
        suffix = os.path.basename(file).split('_')[4].split('.')[0]
        image_info['suffixes'].append(suffix)

        # Determine JPEG file location
        jpg_dir = os.path.join(preview_dir, dirname)
        jpg_filename = os.path.basename(os.path.splitext(file)[0] + '_integ0.jpg')
        jpg_filepath = os.path.join(jpg_dir, jpg_filename)

        # Check that a jpg does not already exist. If it does (and rewrite=False),
        # just call the existing jpg file
        if os.path.exists(jpg_filepath) and not rewrite:
            pass

        # If it doesn't, make it using the preview_image module
        else:
            if not os.path.exists(jpg_dir):
                os.makedirs(jpg_dir)
            im = PreviewImage(file, 'SCI')
            im.output_directory = jpg_dir
            im.make_image()

        # Record how many integrations there are per filetype
        search_jpgs = os.path.join(preview_dir, dirname, file_root + '_{}_integ*.jpg'.format(suffix))
        num_jpgs = len(glob.glob(search_jpgs))
        image_info['num_ints'][suffix] = num_jpgs

        image_info['all_jpegs'].append(jpg_filepath)

    return image_info


def get_instrument_proposals(instrument):
    """Return a list of proposals for the given instrument

    Parameters
    ----------
    instrument : str
        Name of the JWST instrument

    Returns
    -------
    proposals : list
        List of proposals for the given instrument
    """

    service = "Mast.Jwst.Filtered.{}".format(instrument)
    params = {"columns": "program",
              "filters": []}
    response = Mast.service_request_async(service, params)
    results = response[0].json()['data']
    proposals = list(set(result['program'] for result in results))

    return proposals


def get_preview_images_by_instrument(inst):
    """Return a list of preview images available in the filesystem for
    the given instrument.

    Parameters
    ----------
    inst : str
        The instrument of interest (e.g. ``NIRCam``).

    Returns
    -------
    preview_images : list
        A list of preview images available in the filesystem for the
        given instrument.
    """

    # Make sure the instrument is of the proper format (e.g. "Nircam")
    instrument = inst[0].upper() + inst[1:].lower()

    # Query MAST for all rootnames for the instrument
    service = "Mast.Jwst.Filtered.{}".format(instrument)
    params = {"columns": "filename",
              "filters": []}
    response = Mast.service_request_async(service, params)
    results = response[0].json()['data']

    # Parse the results to get the rootnames
    filenames = [result['filename'].split('.')[0] for result in results]

    # Get list of all preview_images
    preview_images = glob.glob(os.path.join(PREVIEW_IMAGE_FILESYSTEM, '*', '*.jpg'))

    # Get subset of preview images that match the filenames
    preview_images = [item for item in preview_images if os.path.basename(item).split('_integ')[0] in filenames]

    return preview_images


def get_preview_images_by_proposal(proposal):
    """Return a list of preview images available in the filesystem for
    the given ``proposal``.

    Parameters
    ----------
    proposal : str
        The five-digit proposal number (e.g. ``88600``).

    Returns
    -------
    preview_images : list
        A list of preview images available in the filesystem for the
        given ``proposal``.
    """

    preview_images = glob.glob(os.path.join(PREVIEW_IMAGE_FILESYSTEM, 'jw{}'.format(proposal), '*'))
    preview_images = [os.path.basename(preview_image) for preview_image in preview_images]

    return preview_images


def get_preview_images_by_rootname(rootname):
    """Return a list of preview images available in the filesystem for
    the given ``rootname``.

    Parameters
    ----------
    rootname : str
        The rootname of interest (e.g. ``jw86600008001_02101_00007_guider2``).

    Returns
    -------
    preview_images : list
        A list of preview images available in the filesystem for the
        given ``rootname``.
    """

    proposal = rootname.split('_')[0].split('jw')[-1][0:5]
    preview_images = sorted(glob.glob(os.path.join(
        PREVIEW_IMAGE_FILESYSTEM,
        'jw{}'.format(proposal),
        '{}*'.format(rootname))))
    preview_images = [os.path.basename(preview_image) for preview_image in preview_images]

    return preview_images


def get_proposal_info(filepaths):
    """Builds and returns a dictionary containing various information
    about the proposal(s) that correspond to the given ``filepaths``.

    The information returned contains such things as the number of
    proposals, the paths to the corresponding thumbnails, and the total
    number of files.

    Parameters
    ----------
    filepaths : list
        A list of full paths to files of interest.

    Returns
    -------
    proposal_info : dict
        A dictionary containing various information about the
        proposal(s) and files corresponding to the given ``filepaths``.
    """

    proposals = list(set([f.split('/')[-1][2:7] for f in filepaths]))
    thumbnail_dir = os.path.join(get_config()['jwql_dir'], 'thumbnails')
    thumbnail_paths = []
    num_files = []
    for proposal in proposals:
        thumbnail_search_filepath = os.path.join(thumbnail_dir, 'jw{}'.format(proposal), 'jw{}*rate*.thumb'.format(proposal))
        thumbnail = glob.glob(thumbnail_search_filepath)
        if len(thumbnail) > 0:
            thumbnail = thumbnail[0]
            thumbnail = '/'.join(thumbnail.split('/')[-2:])
        thumbnail_paths.append(thumbnail)

        fits_search_filepath = os.path.join(FILESYSTEM_DIR, 'jw{}'.format(proposal), 'jw{}*.fits'.format(proposal))
        num_files.append(len(glob.glob(fits_search_filepath)))

    # Put the various information into a dictionary of results
    proposal_info = {}
    proposal_info['num_proposals'] = len(proposals)
    proposal_info['proposals'] = proposals
    proposal_info['thumbnail_paths'] = thumbnail_paths
    proposal_info['num_files'] = num_files

    return proposal_info


def get_thumbnails_by_instrument(inst):
    """Return a list of thumbnails available in the filesystem for the
    given instrument.

    Parameters
    ----------
    inst : str
        The instrument of interest (e.g. ``NIRCam``).

    Returns
    -------
    preview_images : list
        A list of thumbnails available in the filesystem for the
        given instrument.
    """

    # Make sure the instrument is of the proper format (e.g. "Nircam")
    instrument = inst[0].upper() + inst[1:].lower()

    # Query MAST for all rootnames for the instrument
    service = "Mast.Jwst.Filtered.{}".format(instrument)
    params = {"columns": "filename",
              "filters": []}
    response = Mast.service_request_async(service, params)
    results = response[0].json()['data']

    # Parse the results to get the rootnames
    filenames = [result['filename'].split('.')[0] for result in results]

    # Get list of all thumbnails
    thumbnails = glob.glob(os.path.join(THUMBNAIL_FILESYSTEM, '*', '*.thumb'))

    # Get subset of preview images that match the filenames
    thumbnails = [item for item in thumbnails if os.path.basename(item).split('_integ')[0] in filenames]

    return thumbnails


def get_thumbnails_by_proposal(proposal):
    """Return a list of thumbnails available in the filesystem for the
    given ``proposal``.

    Parameters
    ----------
    proposal : str
        The five-digit proposal number (e.g. ``88600``).

    Returns
    -------
    thumbnails : list
        A list of thumbnails available in the filesystem for the given
        ``proposal``.
    """

    thumbnails = glob.glob(os.path.join(THUMBNAIL_FILESYSTEM, 'jw{}'.format(proposal), '*'))
    thumbnails = [os.path.basename(thumbnail) for thumbnail in thumbnails]

    return thumbnails


def get_thumbnails_by_rootname(rootname):
    """Return a list of preview images available in the filesystem for
    the given ``rootname``.

    Parameters
    ----------
    rootname : str
        The rootname of interest (e.g. ``jw86600008001_02101_00007_guider2``).

    Returns
    -------
    thumbnails : list
        A list of preview images available in the filesystem for the
        given ``rootname``.
    """

    proposal = rootname.split('_')[0].split('jw')[-1][0:5]
    thumbnails = sorted(glob.glob(os.path.join(
        THUMBNAIL_FILESYSTEM,
        'jw{}'.format(proposal),
        '{}*'.format(rootname))))

    thumbnails = [os.path.basename(thumbnail) for thumbnail in thumbnails]

    return thumbnails


def split_files(file_list, page_type):
    """JUST FOR USE DURING DEVELOPMENT WITH FILESYSTEM

    Splits the files in the filesystem into "unlooked" and "archived",
    with the "unlooked" images being the most recent 10% of files.
    """
    exp_times = []
    for file in file_list:
        hdr = fits.getheader(file, ext=0)
        exp_start = hdr['EXPSTART']
        exp_times.append(exp_start)

    exp_times_sorted = sorted(exp_times)
    i_cutoff = int(len(exp_times) * .1)
    t_cutoff = exp_times_sorted[i_cutoff]

    mask_unlooked = np.array([t < t_cutoff for t in exp_times])

    if page_type == 'unlooked':
        print('ONLY RETURNING {} "UNLOOKED" FILES OF {} ORIGINAL FILES'.format(len([m for m in mask_unlooked if m]), len(file_list)))
        return [f for i, f in enumerate(file_list) if mask_unlooked[i]]
    elif page_type == 'archive':
        print('ONLY RETURNING {} "ARCHIVED" FILES OF {} ORIGINAL FILES'.format(len([m for m in mask_unlooked if not m]), len(file_list)))
        return [f for i, f in enumerate(file_list) if not mask_unlooked[i]]


def thumbnails(inst, proposal=None):
    """Generate a page showing thumbnail images corresponding to
    activities, from a given ``proposal``

    Parameters
    ----------
    inst : str
        Name of JWST instrument
    proposal : str (optional)
        Number of APT proposal to filter

    Returns
    -------
    dict_to_render : dict
        Dictionary of parameters for the thumbnails
    """

    filepaths = get_filenames_by_instrument(inst)

    # JUST FOR DEVELOPMENT
    # Split files into "archived" and "unlooked"
    if proposal is not None:
        page_type = 'archive'
    else:
        page_type = 'unlooked'
    filepaths = split_files(filepaths, page_type)

    # Determine file ID (everything except suffix)
    # e.g. jw00327001001_02101_00002_nrca1
    full_ids = set(['_'.join(f.split('/')[-1].split('_')[:-1]) for f in filepaths])

    # If the proposal is specified (i.e. if the page being loaded is
    # an archive page), only collect data for given proposal
    if proposal is not None:
        full_ids = [f for f in full_ids if f[2:7] == proposal]

    detectors = []
    proposals = []
    for i, file_id in enumerate(full_ids):
        for file in filepaths:
            if '_'.join(file.split('/')[-1].split('_')[:-1]) == file_id:

                # Parse filename to get program_id
                try:
                    program_id = filename_parser(file)['program_id']
                    detector = filename_parser(file)['detector']
                except ValueError:
                    # Temporary workaround for noncompliant files in filesystem
                    program_id = nfile_id[2:7]
                    detector = file_id[26:]

        # Add parameters to sort by
        if detector not in detectors and not detector.startswith('f'):
            detectors.append(detector)
        if program_id not in proposals:
            proposals.append(program_id)

    # Extract information for sorting with dropdown menus
    # (Don't include the proposal as a sorting parameter if the
    # proposal has already been specified)
    if proposal is not None:
        dropdown_menus = {'detector': detectors}
    else:
        dropdown_menus = {'detector': detectors,
                          'proposal': proposals}

    dict_to_render = {'inst': inst,
                      'tools': MONITORS,
                      'dropdown_menus': dropdown_menus,
                      'prop': proposal}

    return dict_to_render


def thumbnails_ajax(inst, proposal=None):
    """Generate a page that provides data necessary to render the
    ``thumbnails`` template.

    Parameters
    ----------
    inst : str
        Name of JWST instrument
    proposal : str (optional)
        Number of APT proposal to filter

    Returns
    -------
    data_dict : dict
        Dictionary of data needed for the ``thumbnails`` template
    """

    # Get the available files for the instrument
    filepaths = get_filenames_by_instrument(inst)
    if proposal is not None:
        filepaths = split_files(filepaths, 'archive')
    else:
        filepaths = split_files(filepaths, 'unlooked')

    # Get set of unique rootnames
    rootnames = set(['_'.join(f.split('/')[-1].split('_')[:-1]) for f in filepaths])

    # If the proposal is specified (i.e. if the page being loaded is
    # an archive page), only collect data for given proposal
    if proposal is not None:
        rootnames = [rootname for rootname in rootnames if rootname[2:7] == proposal]

    # Initialize dictionary that will contain all needed data
    data_dict = {}
    data_dict['inst'] = inst
    data_dict['file_data'] = {}

    # Gather data for each rootname
    for rootname in rootnames:

        # Parse filename
        try:
            filename_dict = filename_parser(rootname)
        except ValueError:
            # Temporary workaround for noncompliant files in filesystem
            filename_dict = {'activity': file_id[17:19],
                             'detector': file_id[26:],
                             'exposure_id': file_id[20:25],
                             'observation': file_id[7:10],
                             'parallel_seq_id': file_id[16],
                             'program_id': file_id[2:7],
                             'visit': file_id[10:13],
                             'visit_group': file_id[14:16]}

        # Get list of available filenames
        available_files = get_filenames_by_rootname(rootname)

        # Add data to dictionary
        data_dict['file_data'][rootname] = {}
        data_dict['file_data'][rootname]['filename_dict'] = filename_dict
        data_dict['file_data'][rootname]['available_files'] = available_files
        data_dict['file_data'][rootname]['expstart'] = get_expstart(rootname)
        data_dict['file_data'][rootname]['suffixes'] = [filename_parser(filename)['suffix'] for filename in available_files]

    # Extract information for sorting with dropdown menus
    # (Don't include the proposal as a sorting parameter if the
    # proposal has already been specified)
    detectors = [data_dict['file_data'][rootname]['filename_dict']['detector'] for rootname in list(data_dict['file_data'].keys())]
    proposals = [data_dict['file_data'][rootname]['filename_dict']['program_id'] for rootname in list(data_dict['file_data'].keys())]
    if proposal is not None:
        dropdown_menus = {'detector': detectors}
    else:
        dropdown_menus = {'detector': detectors,
                          'proposal': proposals}

    data_dict['tools'] = MONITORS
    data_dict['dropdown_menus'] = dropdown_menus
    data_dict['prop'] = proposal

    return data_dict
