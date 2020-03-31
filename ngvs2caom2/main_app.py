# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#

import importlib
import logging
import os
import sys
import traceback

from caom2 import Observation, CalibrationLevel, DataProductType, ProductType
from caom2 import TemporalWCS, CoordAxis1D, Axis, CoordBounds1D, CoordRange1D
from caom2repo import CAOM2RepoClient
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc


__all__ = ['ngvs_main_app', 'to_caom2', 'update', 'NGVSName', 'COLLECTION',
           'APPLICATION', 'ARCHIVE']


APPLICATION = 'ngvs2caom2'
COLLECTION = 'NGVS'
ARCHIVE = 'NGVS'
INSTRUMENT = 'MegaPrime'

filter_repair_lookup = {'i': 'i_sdss',  # i.MP9703
                        'g': 'g_sdss',  # g.MP9402
                        'r': 'r_sdss',  # r.MP9602
                        'u': 'u_sdss',  # u.MP9302
                        'z': 'z_sdss'}  # z.MP9901
filter_cache = ac.FilterMetadataCache(
    filter_repair_lookup, {}, 'CFHT', {}, 'NONE')


class NGVSName(mc.StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support uncompressed files in storage
    """

    NGVS_NAME_PATTERN = '*'

    def __init__(self, obs_id=None, fname_on_disk=None, file_name=None):
        self.fname_in_ad = file_name
        obs_id = mc.StorageName.remove_extensions(file_name)
        super(NGVSName, self).__init__(
            obs_id, COLLECTION, NGVSName.NGVS_NAME_PATTERN, fname_on_disk)
        self._product_id = NGVSName.get_product_id(file_name)

    @property
    def product_id(self):
        return self._product_id

    def is_valid(self):
        return True

    @staticmethod
    def get_filter_name(f_name):
        bits = f_name.split('.')
        return bits[2]

    @staticmethod
    def get_product_id(f_name):
        if NGVSName.is_catalog(f_name):
            result = 'catalog'
        else:
            bits = f_name.split('.')
            result = f'{bits[1]}.{bits[2]}.{bits[3]}'
        return result

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        return name.replace('.fits', '').replace('.fz', '').replace(
            '.header', '').replace('.sig', '').replace('.weight', '').replace(
            '.cat', '').replace('.mask.rd.reg', '')

    @staticmethod
    def is_catalog(name):
        return '.cat' in name


def get_artifact_product_type(uri):
    result = ProductType.SCIENCE
    if 'weight' in uri or '.sig' in uri:
        result = ProductType.WEIGHT
    elif 'mask.rd.reg' in uri or '.flag' in uri:
        result = ProductType.AUXILIARY
    return result


def get_calibration_level(uri):
    result = CalibrationLevel.PRODUCT
    if NGVSName.is_catalog(uri):
        result = CalibrationLevel.ANALYSIS_PRODUCT
    return result


def get_data_product_type(uri):
    result = DataProductType.IMAGE
    if NGVSName.is_catalog(uri):
        result = DataProductType.CATALOG
    return result


def get_proposal_id(header):
    # JJK - 23-03-20
    # Replace with value of ProposalID from the first input ‘member’  from list
    # TODO
    return None


def get_provenance_version(uri):
    bits = mc.CaomName(uri).file_name.split('.')
    return bits[3]


def get_target_name(uri):
    return mc.CaomName(uri).file_name.split('.')[0]


def _informative_uri(uri):
    result = False
    if ('.weight' not in uri and '.sig' not in uri and '.cat' not in uri and
            '.mask' not in uri):
        # all the excluded names have fewer useful keywords
        result = True
    return result


def accumulate_bp(bp, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug(f'Begin accumulate_bp for {uri}.')
    bp.configure_position_axes((1, 2))

    # they're all CompositeObservations
    bp.set('CompositeObservation.members', {})

    # JJK - comments by email - 28-03-20
    if _informative_uri(uri):
        bp.clear('Observation.environment.seeing')
        bp.add_fits_attribute('Observation.environment.seeing', 'FINALIQ')
        bp.set('Observation.proposal.id', 'get_proposal_id(header)')

        bp.clear('Plane.metrics.magLimit')
        bp.add_fits_attribute('Plane.metrics.magLimit', 'MAGLIM')

        bp.clear('Plane.provenance.keywords')
        bp.add_fits_attribute('Plane.provenance.keywords', 'COMBINET')

        bp.clear('Chunk.position.resolution')
        bp.add_fits_attribute('Chunk.position.resolution', 'FINALIQ')

    bp.set('Observation.instrument.name', 'MegaPrime')

    bp.set('Observation.proposal.pi', 'Laura Ferrarese')
    bp.set('Observation.proposal.project', 'NGVS')
    bp.set('Observation.proposal.title',
           'Next Generation Virgo Cluster Survey')
    bp.set('Observation.proposal.keywords', 'Galaxy Cluster Dwarfs')

    bp.set('Observation.target.name', 'get_target_name(uri)')
    bp.set('Observation.target.type', 'field')

    bp.set('Observation.telescope.name', 'CFHT 3.6m')
    x, y, z = ac.get_geocentric_location('cfht')
    bp.set('Observation.telescope.geoLocationX', x)
    bp.set('Observation.telescope.geoLocationY', y)
    bp.set('Observation.telescope.geoLocationZ', z)

    bp.set('Plane.calibrationLevel', 'get_calibration_level(uri)')
    bp.set('Plane.dataProductType', 'get_data_product_type(uri)')
    bp.set('Plane.dataRelease', '2022-01-01T00:00:00')

    bp.set('Plane.provenance.name', 'MEGAPIPE')
    bp.set('Plane.provenance.producer', 'CADC')
    bp.set('Plane.provenance.project', 'NGVS')
    bp.set('Plane.provenance.reference',
           'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/ngvs/docs/'
           'ngvsdoc.html')
    bp.set('Plane.provenance.version', 'get_provenance_version(uri)')

    bp.set('Artifact.productType', 'get_artifact_product_type(uri)')

    logging.debug('Done accumulate_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes, must
    have this signature for import_module loading and execution.

    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.debug('Begin update.')
    mc.check_param(observation, Observation)
    fqn = kwargs.get('fqn')
    headers = kwargs.get('headers')
    uri = kwargs.get('uri')

    if uri is not None:
        f_name = mc.CaomName(uri).file_name
        ngvs_name = NGVSName(file_name=f_name)
    elif fqn is not None:
        ngvs_name = NGVSName(file_name=os.path.basename(fqn))
    else:
        raise mc.CadcException(f'Cannot define a NGVSName instance for '
                               f'{observation.observation_id}')

    for plane in observation.planes.values():
        if plane.product_id != ngvs_name.product_id:
            continue
        if _informative_uri(fqn):
            cc.update_plane_provenance_single(plane, headers, 'HISTORY',
                                              'CFHT',
                                              _repair_history_provenance_value,
                                              observation.observation_id)
            for artifact in plane.artifacts.values():
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        _update_energy(chunk, artifact.uri,
                                       observation.observation_id)
                        # _update_time(chunk, headers[0], plane.provenance,
                        #              observation.observation_id)

    cc.update_observation_members(observation)
    logging.debug('Done update.')
    return observation


def _update_energy(chunk, uri, obs_id):
    logging.debug(f'Begin _update_energy for {obs_id}')
    filter_name = NGVSName.get_filter_name(uri)
    filter_md = filter_cache.get_svo_filter(INSTRUMENT, filter_name)
    cc.build_chunk_energy_range(chunk, filter_name, filter_md)
    logging.debug(f'End _update_energy.')


def _update_time(chunk, header, provenance, obs_id):
    logging.debug(f'Begin _update_time for {obs_id}')

    if chunk is not None and provenance is not None:
        # bounds = ctor
        config = mc.Config()
        config.get_executors()
        subject = mc.define_subject(config)
        client = CAOM2RepoClient(
            subject, config.logging_level, config.resource_id)
        metrics = mc.Metrics(config)
        bounds = CoordBounds1D()
        min_date = 0
        max_date = sys.float_info.max
        exposure = 0
        for entry in provenance.inputs:
            ip_obs_id, ip_product_id = mc.CaomName.decompose_provenance_input(
                entry)
            ip_obs = mc.repo_get(client, COLLECTION, ip_obs_id, metrics)
            ip_plane = ip_obs.planes.get(ip_product_id)
            if (ip_plane is not None and ip_plane.time is not None and
                    ip_plane.time.bounds is not None):
                bounds.samples.append(CoordRange1D(ip_plane.time.bounds.lower,
                                                   ip_plane.time.bounds.upper))
                min_date = min(ip_plane.time.bounds.lower, min_date)
                max_date = max(ip_plane.time.bounds.upper, max_date)
                exposure += ip_plane.time.exposure
        axis = Axis(ctype='TIME', cunit='mjd')
        time_axis = CoordAxis1D(axis=axis,
                                error=None,
                                range=None,
                                bounds=bounds,
                                function=None)
        temporal_wcs = TemporalWCS(axis=time_axis, timesys=None, trefpos=None,
                                   mjdref=None, exposure=exposure,
                                   resolution=None)
        chunk.time_axis = 3
        chunk.time = temporal_wcs
    logging.debug(f'End _update_time.')


def _repair_history_provenance_value(value, obs_id):
    logging.debug(f'Begin _repair_history_provenance_value for {obs_id}')
    results = []
    # HISTORY headers with provenance:
    # HISTORY = input image 973887p.fits; phot ref: SDSS; IQ=0.55; Sky= 1741.0
    # HISTORY = input image 973888p.fits; phot ref: SDSS; IQ=0.56; Sky= 1740.0
    # HISTORY = input image 973889p.fits; phot ref: SDSS; IQ=0.54; Sky= 1728.0
    # HISTORY = input image 973890p.fits; phot ref: SDSS; IQ=0.50; Sky= 1700.0
    # HISTORY = input image 973891p.fits; phot ref: SDSS; IQ=0.53; Sky= 1675.0
    if 'input image' in str(value):
        for entry in value:
            if 'input image' in entry:
                temp = str(entry).split('input image ')
                # TODO TODO TODO - this is the old archive product id and
                #  plane id - whoops, when it comes to un-intended consequences
                prov_prod_id = temp[1].split(';')[0].replace('.fits', '')
                prov_obs_id = prov_prod_id[:-1]
                # 0 - observation
                # 1 - plane
                results.append([prov_obs_id, prov_prod_id])
    logging.debug(f'End _repair_history_provenance_value')
    return results


def _build_blueprints(uris):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param uris The artifact URIs for the files to be processed."""
    module = importlib.import_module(__name__)
    blueprints = {}
    for uri in uris:
        blueprint = ObsBlueprint(module=module)
        if not mc.StorageName.is_preview(uri):
            accumulate_bp(blueprint, uri)
        blueprints[uri] = blueprint
    return blueprints


def _get_uris(args):
    result = []
    if args.lineage:
        for ii in args.lineage:
            result.append(ii.split('/', 1)[1])
    elif args.local:
        for ii in args.local:
            file_uri = mc.build_uri(ARCHIVE, os.path.basename(ii))
            result.append(file_uri)
    else:
        raise mc.CadcException(
            f'Could not define uri from these args {args}')
    return result


def to_caom2():
    args = get_gen_proc_arg_parser().parse_args()
    uris = _get_uris(args)
    blueprints = _build_blueprints(uris)
    result = gen_proc(args, blueprints)
    logging.debug(f'Done {APPLICATION} processing.')
    return result
           

def ngvs_main_app():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        result = to_caom2()
        sys.exit(result)
    except Exception as e:
        logging.error(f'Failed {APPLICATION} execution for {args}.')
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
