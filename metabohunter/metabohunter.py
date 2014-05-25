# -*- coding: utf-8 -*-

import os
import sys
import re
import logging

from collections import OrderedDict

import csv
import requests

try:
    import IPython
except:
    IPython = False
else:
    from IPython.html import widgets
    from IPython.display import display

METABOTYPE_ALL = 'All'
METABOTYPE_DRUG = 'Drug'
METABOTYPE_FOOD_ADDITIVE = 'Food additive'
METABOTYPE_MAMMALIAN = 'Mammalian'
METABOTYPE_MICROBIAL = 'Microbial'
METABOTYPE_PLANT = 'Plant'
METABOTYPE_SYNTHETIC = 'Synthetic/Industrial chemical'

METABOTYPES = [ METABOTYPE_ALL, METABOTYPE_DRUG, METABOTYPE_FOOD_ADDITIVE, 
                METABOTYPE_MAMMALIAN, METABOTYPE_MICROBIAL, METABOTYPE_PLANT, 
                METABOTYPE_SYNTHETIC ]
    
DATABASE_HMDB = 'HMDB'
DATABASE_MMCD = 'MMCD'

DATABASES = [ DATABASE_HMDB, DATABASE_MMCD ]

SAMPLE_PH10_11 = 'ph7'
SAMPLE_PH7_10 = 'ph7'
SAMPLE_PH6_7 = 'ph6'
SAMPLE_PH5_6 = 'ph5'
SAMPLE_PH4_5 = 'ph4'
SAMPLE_PH3_4 = 'ph3'

SAMPLE_PHS = [ SAMPLE_PH10_11, SAMPLE_PH7_10, SAMPLE_PH6_7, SAMPLE_PH5_6, 
               SAMPLE_PH3_4, SAMPLE_PH4_5 ]

SOLVENT_ALL = 'all'
SOLVENT_WATER = 'water'
SOLVENT_CDCL3 = 'cdcl3'
SOLVENT_CD3OD = '5d30d'
SOLVENT_5PC_DMSO = '5dmso'

SOLVENTS = [ SOLVENT_ALL, SOLVENT_WATER, SOLVENT_CDCL3, SOLVENT_CD3OD, SOLVENT_5PC_DMSO ]

FREQUENCY_ALL = 'all'
FREQUENCY_600 = '600'
FREQUENCY_500 = '500'
FREQUENCY_400 = '400'

FREQUENCIES = [ FREQUENCY_ALL, FREQUENCY_600, FREQUENCY_500, FREQUENCY_400 ]

METHOD_HIGHEST_NUMBER = 'HighestNumber'
METHOD_HIGHEST_NUMBER_NEIGHBOURHOOD = 'HighestNumberNeighbourhood'
METHOD_GREEDY = 'Greedy2'
METHOD_HIGHEST_NUMBER_HEIGHTS = 'HighestNumberHeights'
METHOD_GREEDY_HEIGHTS = 'Greedy2Heights'

METHODS = [ METHOD_HIGHEST_NUMBER, METHOD_HIGHEST_NUMBER_NEIGHBOURHOOD, METHOD_GREEDY, 
            METHOD_HIGHEST_NUMBER_HEIGHTS, METHOD_GREEDY_HEIGHTS ]

DEFAULT_NOISE_THRESHOLD = 0.0
DEFAULT_CONFIDENCE_THRESHOLD = 0.4
DEFAULT_TOLERANCE = 0.1

class MetaboHunterDataIncorrectType(Exception):
    pass

class MetaboHunterDataLengthMismatch(Exception):
    pass

class MetaboHunterInvalidParameter(Exception):
    pass


def request( ppms, 
                  peaks, 
                  metabotype=METABOTYPE_ALL,
                  database=DATABASE_HMDB,
                  ph=SAMPLE_PH7_10,
                  solvent=SOLVENT_WATER,
                  frequency=FREQUENCY_600,
                  method=METHOD_HIGHEST_NUMBER_NEIGHBOURHOOD,
                  noise=DEFAULT_NOISE_THRESHOLD,
                  confidence=DEFAULT_CONFIDENCE_THRESHOLD,
                  tolerance=DEFAULT_TOLERANCE):

    # Validate inputs before we start requesting
    
    # Flattent to 1D
    if type(ppms) != list:
        ppms = ppms.squeeze()
        if len(ppms.shape) != 1:
            raise MetaboHunterDataIncorrectType("Array 'ppm' should be a 1D numpy array or list")
        ppms = list( ppms ) 
    
    
    if type(peaks) != list:
        peaks = peaks.squeeze()
        if len(peaks.shape) != 1:
            raise MetaboHunterDataIncorrectType("Array 'peaks' should be a 1D numpy array or list")
        peaks = list( peaks ) 
    
        
    # Check peaks and ppms are equal length
    if len(ppms) != len(peaks):
        raise MetaboHunterDataLengthMismatch("Inputs 'ppm' and 'peaks' have mismatched dimensions: %s and %s" % ( len(ppms), len(peaks) ) )        

    # We should be good to go from here
    

    ### GLOBAL VARIABLES ###
    remote_data = dict()  # Container for references to metabolite data on remote server

    # Web service peak-list assignment (metabohunter)
    logging.info("Sending peaklist to MetaboHunter...")

    peaks_list = '\n'.join([' '.join([str(a), str(b)]) for a, b in zip(ppms, peaks)])

    url = 'http://www.nrcbioinformatics.ca/metabohunter/post_handler.php'

    values = {  # 'file'          : open('metabid_peaklist.txt', 'rt'),
              'posturl': 'upload_file.php',
              'useall': 'yes',
              'peaks_list': peaks_list,
              'dbsource': database,
              'metabotype': metabotype,
              'sampleph': ph,
              'solvent': solvent,
              'freq': frequency,
              'method': method,
              'noise': noise,
              'thres': confidence,
              'neighbourhood': tolerance,  # tolerance, # Use same tolerance as for shift
              'submit': 'Find matches',
             }

    try:
        r = requests.post(url, data=values)
    except e:
        raise e

    html = r.content

    m = re.search('name="hits" value="(.*?)\n(.*?)\n"', html, re.MULTILINE | re.DOTALL)
    remote_data['metabolite_table'] = m.group(2)

    m = re.search('name="sample_file" value="(.*?)"', html, re.MULTILINE | re.DOTALL)
    remote_data['sample_file'] = m.group(1)

    logging.info("Received analysis from MetaboHunter, interpreting...")

    # Regexp out the matched peaks table from the hidden form field in response (easiest to parse)
    metabolites = OrderedDict()
    
    # Iterate line by line (skip first, header) building a table of the returned metabolites
    for row in remote_data['metabolite_table'].split('\n'):

        fields = row.split('\t')  # split row on tabs
        m = re.match("(.*?) \((\d*?)/(\d*?)\)", fields[3])

        metabolites[fields[1]] = {
            'name': fields[2],
            'score': float(m.group(1)),
            'peaks': "%d/%d" % (int(m.group(2)), int(m.group(3))),
            'taxnomic': fields[4],
            'rank': fields[0],
            }

    logging.info("Retrieving matched peaks to metabolite relationships...")

    values = {'sample_file': remote_data['sample_file'],
              'matched_peaks_file': remote_data['sample_file'] + "_matched_spectra.txt",
              'noise': noise,
              'hits': 'Rank\tID\tMetabolite name\tMatching peaks score\tTaxonomic origin\r\n' + remote_data['metabolite_table'] + '\r\n',
     }

    # Create the Request object
    url = 'http://www.nrcbioinformatics.ca/metabohunter/download_matched_peaks.php'
    try:
        r = requests.post(url, data=values, files={'foo': 'bar'})
    except e:
        raise e

    matched_peaks_text = r.content

    logging.info("Extracting data...")

    # Need to do this in two steps, so they are in the correct order for output
    metabolite_peaks = dict()
    matched_peak_metabolites = dict()

    for row in matched_peaks_text.split('\n'):
        fields = row.split()
        if fields:
            # fields[0] contains the HMDBid plus a colon :(
            fields[0] = fields[0].rstrip(':')
            metabolite_peaks[fields[0]] = fields[1:]

    for metabolite in metabolites:
        if metabolite in metabolite_peaks:
            # Save metabolite for each peak
            for peak in metabolite_peaks[metabolite]:
                #if peak in matched_peak_metabolites:
                #    matched_peak_metabolites[ peak ].append(metabolite)
                #else:
                matched_peak_metabolites[peak] = metabolite

    matched_hmdbs = []
    # Add matched HMDBs to a list to return or None where no match found
    # Returned list will be the same length as input ppms/peaks
    for n, p in enumerate(ppms):
        sp2 = '%.2f' % p # Returned peaks are at 2dp so we need to match on this basis
        if sp2 in matched_peak_metabolites:
            hmdbid = matched_peak_metabolites[sp2]
            matched_hmdbs.append(hmdbid)
        else:
            matched_hmdbs.append(None)

    return matched_hmdbs
    
if IPython:
    # We've found IPython installed. Define a wrapper object to generate widgets for controlling
    # the metabohunter request object
    
    class IPyMetaboHunter(object):
        
        def __init__(self, *args, **kwargs):
            
            self.defaults = dict(
                  metabotype=METABOTYPE_ALL,
                  database=DATABASE_HMDB,
                  ph=SAMPLE_PH7_10,
                  solvent=SOLVENT_WATER,
                  frequency=FREQUENCY_600,
                  method=METHOD_HIGHEST_NUMBER_NEIGHBOURHOOD,
                  noise=DEFAULT_NOISE_THRESHOLD,
                  confidence=DEFAULT_CONFIDENCE_THRESHOLD,
                  tolerance=DEFAULT_TOLERANCE
            )
            
            # Overwrite the settings using kwargs passed on create
            self.settings = dict( self.defaults.items() + kwargs.items() )

            # Set up widgets
            self.form = widgets.ContainerWidget()
            self.form.set_css({'width':'20ex'}, selector='.widget-hlabel')
            
            w_metabotype = widgets.DropdownWidget(description="Metabotype:", value=self.settings['metabotype'], values=METABOTYPES)
            w_metabotype.on_trait_change(lambda t, value: self.set_value('metabotype', value), 'value' )
            
            w_database = widgets.DropdownWidget(description="Database:", value=self.settings['database'], values=DATABASES)
            w_database.on_trait_change(lambda t, value: self.set_value('database', value), 'value' )
            
            w_ph = widgets.DropdownWidget(description="pH:", value=self.settings['ph'], values=SAMPLE_PHS)
            w_ph.on_trait_change(lambda t, value: self.set_value('ph', value), 'value' )

            w_solvent = widgets.DropdownWidget(description="Solvent:", value=self.settings['solvent'], values=SOLVENTS)
            w_solvent.on_trait_change(lambda t, value: self.set_value('solvent', value), 'value' )
            
            w_frequency = widgets.DropdownWidget(description="Frequency MHz:", value=self.settings['frequency'], values=FREQUENCIES)
            w_frequency.on_trait_change(lambda t, value: self.set_value('frequency', value), 'value' )
            
            w_method = widgets.DropdownWidget(description="Method:", value=self.settings['method'], values=METHODS)
            w_method.on_trait_change(lambda t, value: self.set_value('method', value), 'value' )
            
            w_noise = widgets.FloatSliderWidget(description="Noise threshold:", value=self.settings['noise'], min=0., max=10., step=0.01)
            w_noise.on_trait_change(lambda t, value: self.set_value('noise', value), 'value' )

            w_confidence = widgets.FloatSliderWidget(description="Confidence threshold:", value=self.settings['confidence'], min=0., max=1., step=0.1)
            w_confidence.on_trait_change(lambda t, value: self.set_value('confidence', value), 'value' )

            w_tolerance = widgets.FloatSliderWidget(description="Shift tolerance:", value=self.settings['tolerance'], min=0., max=10., step=0.01)
            w_tolerance.on_trait_change(lambda t, value: self.set_value('tolerance', value), 'value' )

            self.form.children = [ w_metabotype, w_database, w_ph, w_solvent, w_frequency, 
                                   w_method, w_noise, w_confidence, w_tolerance]
            
        def set_value(self, name, value):
            self.settings[name] = value

        @property
        def kwargs(self):
            return self.settings
            
        def display(self):
            display( self.form )


