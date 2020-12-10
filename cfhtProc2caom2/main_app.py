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

from caom2repo import CAOM2RepoClient
from caom2 import Observation, CalibrationLevel, DataProductType, ProductType
from caom2 import TemporalWCS, CoordAxis1D, Axis, CoordBounds1D, CoordRange1D
from caom2 import RefCoord, CoordFunction1D, RefCoord
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from caom2pipe import translate_composable as tc
from cfhtProc2caom2 import storage_names as sn


__all__ = ['cfht_proc_main_app', 'to_caom2', 'update', 'APPLICATION']


APPLICATION = 'cfhtProc2caom2'
INSTRUMENT = 'MegaPrime'

# filter names have to conform to those of the SVO filter service, because
# they're going to be used to create a URL to query from there
FILTER_REPAIR_LOOKUP = {'i': 'i_sdss',  # i.MP9703
                        'I2': 'i',      # i.MP9702
                        'I': 'i1',      # i.MP9701
                        'g': 'g_sdss',  # g.MP9402
                        'G': 'g',       # g.MP9401
                        'GRI': 'gri',   # gri.MP9605
                        'r': 'r_sdss',  # r.MP9602
                        'R': 'r',       # r.MP9601
                        'u': 'u_sdss',  # u.MP9302
                        'U': 'u',       # u.MP9301
                        'z': 'z_sdss',  # z.MP9901
                        'Z': 'z'}       # z.MP9801
filter_cache = ac.FilterMetadataCache(
    FILTER_REPAIR_LOOKUP, {}, 'CFHT', {}, 'NONE')

# because 1 is never enough ... and also because, in an effort to keep
# the number of filter names smaller, they have to have the case of
# CFHT, which is filter name letters are lower-case, but CFHT parts
# numbers are upper-case.
CAOM_FILTER_REPAIR_LOOKUP = {'I.MP9703': 'i.MP9703',
                             'I.MP9702': 'i.MP9702',
                             'I.MP9701': 'i.MP9701',
                             'G.MP9402': 'g.MP9402',
                             'G.MP9401': 'g.MP9401',
                             'GRI.MP9605': 'gri.MP9605',
                             'R.MP9602': 'r.MP9602',
                             'R.MP9601': 'r.MP9601',
                             'U.MP9302': 'u.MP9302',
                             'U.MP9301': 'u.MP9301',
                             'Z.MP9901': 'z.MP9901',
                             'Z.MP9801': 'z.MP9801'}


def accumulate_bp(bp, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model
    Observation level."""
    logging.debug(f'Begin accumulate_bp for {uri}.')
    bp.configure_position_axes((1, 2))

    scheme, archive, file_name = mc.decompose_uri(uri)
    storage_name = sn.get_storage_name(file_name, file_name)
    if sn.is_ngvs(uri):
        _accumulate_ngvs_bp(bp, storage_name)
    else:
        _accumulate_mp_bp(bp, storage_name)

    # they're all DerivedObservations
    bp.set('DerivedObservation.members', {})
    bp.set('Observation.type', 'OBJECT')

    bp.set('Observation.proposal.id', 'get_proposal_id(header)')

    bp.clear('Plane.metaRelease')
    bp.add_fits_attribute('Plane.metaRelease', 'REL_DATE')

    bp.clear('Chunk.position.resolution')
    bp.add_fits_attribute('Chunk.position.resolution', 'FINALIQ')

    bp.set('Observation.instrument.name', INSTRUMENT)
    bp.set('Observation.telescope.name', 'CFHT 3.6m')
    x, y, z = ac.get_geocentric_location('cfht')
    bp.set('Observation.telescope.geoLocationX', x)
    bp.set('Observation.telescope.geoLocationY', y)
    bp.set('Observation.telescope.geoLocationZ', z)

    bp.set('Plane.calibrationLevel', 'get_calibration_level(uri)')
    bp.set('Plane.dataProductType', 'get_data_product_type(uri)')
    bp.set('Plane.provenance.producer', 'CADC')

    bp.set('Artifact.productType', 'get_artifact_product_type(uri)')

    if storage_name.collection == sn.MP_COLLECTION:
        _accumulate_mp_bp(bp, storage_name)
    else:
        _accumulate_ngvs_bp(bp, storage_name)

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
        storage_name = sn.get_storage_name(uri, uri)
    elif fqn is not None:
        temp = os.path.basename(fqn)
        storage_name = sn.get_storage_name(temp, temp)
    else:
        raise mc.CadcException(f'Cannot define a MEGAPIPEName instance for '
                               f'{observation.observation_id}')
    if headers is None:
        logging.warning(f'No metadata for {storage_name.file_name}')
        return observation

    logging.debug(f'Update for {observation.observation_id} with '
                  f'{storage_name.file_name}.')

    max_meta_release = observation.meta_release
    min_seeing = None
    if (observation.environment is not None and
            observation.environment.seeing is not None):
        min_seeing = observation.environment.seeing
    if not storage_name.is_catalog:
        for plane in observation.planes.values():
            max_meta_release = _update_release_date(
                plane, max_meta_release, headers)
            if plane.product_id != storage_name.product_id:
                continue
            min_seeing = _minimize(min_seeing,
                                   _get_keyword(headers, 'FINALIQ'))
            for artifact in plane.artifacts.values():
                if artifact.uri != storage_name.file_uri:
                    continue
                if (artifact.product_type is ProductType.WEIGHT and
                        storage_name.collection == sn.MP_COLLECTION):
                    artifact.parts = None
                    continue
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        if _informative_uri(storage_name.file_name):
                            _update_energy(chunk, headers, storage_name,
                                           observation.observation_id)
                        if not storage_name.is_catalog:
                            if storage_name.collection == sn.MP_COLLECTION:
                                if chunk.position is not None:
                                    chunk.position.resolution = None

                        # SGw - 24-11-20 - set observable/axis to None: there
                        # is no such information in the files
                        if chunk.observable_axis is not None:
                            chunk.observable_axis = None
                        if chunk.observable is not None:
                            chunk.observable = None

            if (_informative_uri(storage_name.file_name) and
                    plane.provenance is not None):
                cc.append_plane_provenance_single(
                    plane, headers, 'HISTORY', 'CFHT',
                    _repair_history_provenance_value,
                    observation.observation_id)
                if plane.provenance.run_id == 'None':
                    plane.provenance.run_id = None
                if (plane.provenance.keywords is not None and
                        'None' in plane.provenance.keywords):
                    plane.provenance.keywords.remove('None')

            # _update_ngvs_time is dependent on provenance information that is
            # generated right before this
            if (storage_name.collection == sn.NGVS_COLLECTION and
                    not storage_name.is_catalog):
                for artifact in plane.artifacts.values():
                    if artifact.uri != storage_name.file_uri:
                        continue
                    for part in artifact.parts.values():
                        for chunk in part.chunks:
                            _update_ngvs_time(chunk, plane.provenance,
                                              observation.observation_id)
    observation.meta_release = max_meta_release
    if observation.environment is not None:
        observation.environment.seeing = min_seeing
    if (observation.target is not None and
            storage_name.collection == sn.MP_COLLECTION):
        observation.target.standard = False
    cc.update_observation_members(observation)
    logging.debug('Done update.')
    return observation


def _accumulate_mp_bp(bp, storage_name):
    logging.debug(f'Begin _accumulate_mp_bp for {storage_name.file_name}')

    bp.clear('Observation.metaRelease')
    bp.add_fits_attribute('Observation.metaRelease', 'DATE')
    bp.add_fits_attribute('Observation.metaRelease', 'REL_DATE')

    bp.set('Observation.algorithm.name', 'MEGAPIPE')

    bp.set('Observation.environment.photometric', True)
    bp.clear('Observation.environment.seeing')
    bp.add_fits_attribute('Observation.environment.seeing', 'FINALIQ')

    bp.set('Observation.proposal.id', 'get_proposal_id(header)')

    bp.clear('Plane.metaRelease')
    bp.add_fits_attribute('Plane.metaRelease', 'DATE')
    bp.add_fits_attribute('Plane.metaRelease', 'REL_DATE')

    bp.clear('Chunk.position.resolution')
    bp.add_fits_attribute('Chunk.position.resolution', 'FINALIQ')

    bp.clear('Plane.dataRelease')
    bp.add_fits_attribute('Plane.dataRelease', 'DATE')
    bp.add_fits_attribute('Plane.dataRelease', 'REL_DATE')

    bp.set('Plane.provenance.name', 'MEGAPIPE')
    bp.set('Plane.provenance.producer', 'CADC')
    bp.set('Plane.provenance.version', '2.0')
    bp.set('Plane.provenance.reference',
           'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/')

    if storage_name.is_preview:
        bp.set('Artifact.productType', ProductType.PREVIEW)

    bp.clear('Chunk.position.coordsys')
    bp.add_fits_attribute('Chunk.position.coordsys', 'RADECSYS')

    bp.clear('Observation.target.name')
    if not storage_name.is_weight:
        bp.clear('Plane.provenance.lastExecuted')
        bp.add_fits_attribute('Plane.provenance.lastExecuted', 'DATE')
        bp.add_fits_attribute('Observation.target.name', 'OBJECT')
        bp.clear('Plane.metrics.magLimit')
        bp.add_fits_attribute('Plane.metrics.magLimit', 'ML_5SIGA')

    logging.debug(f'End _accumulate_mp_bp.')


def _accumulate_ngvs_bp(bp, storage_name):
    logging.debug(f'Begin _accumulate_ngvs_bp for {storage_name.file_name}.')
    bp.clear('Observation.algorithm.name')

    # NGVS
    # JJK - comments by email - 28-03-20
    # if _informative_uri(uri):
    # make sure information set from header keywords is only set for
    # the fits files where it's accessible
    bp.set('Observation.metaRelease', '2022-01-01 00:00:00')

    bp.set('Observation.proposal.pi', 'Laura Ferrarese')
    bp.set('Observation.proposal.project', 'NGVS')
    bp.set('Observation.proposal.title',
           'Next Generation Virgo Cluster Survey')
    bp.set('Observation.proposal.keywords', 'Galaxy Cluster Dwarfs')

    bp.set('Observation.target.name', 'get_target_name(uri)')
    if 'weight' not in storage_name.file_uri:
        bp.set('Observation.target.type', 'field')

    bp.set('Plane.dataRelease', '2022-01-01T00:00:00')
    bp.set('Plane.metaRelease', '2022-01-01T00:00:00')
    # SGw 28-07-20
    # The last 4 characters in the ID column (l128 vs g002 vs g004) refer to
    # variation in the pipeline. If they could be captured as
    # plane.provenance.name = MegaPipe_g004 that would be great.
    bp.set('Plane.provenance.name', f'MegaPipe_{storage_name.version}')
    bp.clear('Plane.provenance.keywords')
    bp.add_fits_attribute('Plane.provenance.keywords', 'COMBINET')
    bp.set('Plane.provenance.project', 'NGVS')
    bp.set('Plane.provenance.reference',
           'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/ngvs/'
           'docs/ngvsdoc.html')
    bp.set('Plane.provenance.version', 'get_provenance_version(uri)')

    bp.set('Artifact.productType', 'get_artifact_product_type(uri)')

    bp.set('Chunk.energy.bandpassName', 'get_ngvs_bandpass_name(uri)')

    bp.clear('Plane.metrics.magLimit')
    bp.add_fits_attribute('Plane.metrics.magLimit', 'MAGLIM')

    logging.debug(f'End _accumulate_ngvs_bp.')


def _get_keyword(headers, keyword):
    result = headers[0].get(keyword)
    if result is None and len(headers) > 1:
        result = headers[1].get(keyword)
    return result


def get_artifact_product_type(uri):
    result = ProductType.SCIENCE
    if 'weight' in uri or '.sig' in uri:
        result = ProductType.WEIGHT
    elif 'mask.rd.reg' in uri or '.flag' in uri:
        result = ProductType.AUXILIARY
    return result


def get_calibration_level(uri):
    result = CalibrationLevel.PRODUCT
    storage_name = sn.get_storage_name(uri, uri)
    if storage_name.is_catalog:
        result = CalibrationLevel.ANALYSIS_PRODUCT
    return result


def get_data_product_type(uri):
    result = DataProductType.IMAGE
    storage_name = sn.get_storage_name(uri, uri)
    if storage_name.is_catalog:
        # PD 09-12-20
        # I need to modify the ObsCore view to filter the observations with
        # DataProductType.catalog out (not compliant to spec) but there is a
        # different value measurements that means roughly the same thing.
        # There should be a DataProductType constant declared in the py library
        # for this value.
        result = DataProductType.MEASUREMENTS
    return result


def get_ngvs_bandpass_name(uri):
    reverse_filter_lookup = {
        'i': 'i.MP9703',
        'g': 'g.MP9402',
        'r': 'r.MP9602',
        'u': 'u.MP9302',
        'z': 'z.MP9901'}
    storage_name = sn.get_storage_name(uri, uri)
    result = None
    if storage_name.filter_name is not None:
        result = reverse_filter_lookup.get(storage_name.filter_name)
    return result


def get_proposal_id(header):
    # JJK - 23-03-20
    # Replace with value of ProposalID from the first input ‘member’  from list
    # TODO
    return None


def get_provenance_version(uri):
    storage_name = sn.get_storage_name(uri, uri)
    return storage_name.version


def get_target_name(uri):
    return mc.CaomName(uri).file_name.split('.')[0]


def _informative_uri(uri):
    result = False
    if ('.weight' not in uri and '.sig' not in uri and '.cat' not in uri and
            '.mask' not in uri and '.flag' not in uri):
        # all the excluded names have fewer useful keywords
        result = True
    return result


def _finish_catalog_plane(observation, plane):
    logging.debug(f'Begin _finish_catalog_plane for '
                  f'{observation.observation_id}.')
    plane.meta_release = observation.meta_release
    logging.debug('Done _finish_catalog_plane.')


def _minimize(x, candidate):
    result = x
    if candidate is not None:
        if x is None:
            result = candidate
        else:
            result = min(x, candidate)
    return result


def _update_energy(chunk, headers, storage_name, obs_id):
    logging.debug(f'Begin _update_energy for {obs_id}')
    filter_name = storage_name.filter_name
    filter_md = filter_cache.get_svo_filter(INSTRUMENT, filter_name)
    cc.build_chunk_energy_range(chunk, filter_name, filter_md)
    if storage_name.collection == sn.MP_COLLECTION:
        temp = _get_keyword(headers, 'FILTER')
        # an attempt to keep the number of unique filter names lower
        chunk.energy.bandpass_name = CAOM_FILTER_REPAIR_LOOKUP.get(temp)
    if storage_name.collection == sn.MP_COLLECTION:
        chunk.energy.resolving_power = None
    logging.debug(f'End _update_energy.')


def _update_observation_metadata(obs, headers, ngvs_name, uri):
    """
    Why this method exists:

    The NGVS weight files have almost no metadata in the primary HDU, but
    all the needed metadata in subsequent HDUs.

    It's not possible to apply extension
    numbers for non-chunk blueprint entries, so that means that to use the
    information captured in the blueprint, the header that's provided
    must be manipulated instead. There is only access to the header
    information in this extension of the fitscaom2 module (i.e. this file)
    during the execution of the 'update' step of fits2caom2 execution.
    Hence the self-referential implementation. Maybe it will result in an
    infinite loop and I'll be sad.
    """
    n_axis = headers[0].get('NAXIS')
    if n_axis == 0:
        logging.warning(f'Resetting the header/blueprint relationship for '
                        f'{ngvs_name.file_name} in {obs.observation_id}')
        module = importlib.import_module(__name__)
        blueprint = ObsBlueprint(module=module)
        accumulate_bp(blueprint, uri)
        tc.add_headers_to_obs_by_blueprint(
            obs, [headers[1]], blueprint, uri, ngvs_name.product_id)


def _update_ngvs_time(chunk, provenance, obs_id):
    logging.debug(f'Begin _update_ngvs_time for {obs_id}')
    if (chunk is not None and provenance is not None and
            len(provenance.inputs) > 0):
        # bounds = ctor
        config = mc.Config()
        config.get_executors()
        subject = mc.define_subject(config)
        client = CAOM2RepoClient(
            subject, config.logging_level, 'ivo://cadc.nrc.ca/ams')
        metrics = mc.Metrics(config)
        bounds = CoordBounds1D()
        min_date = 0
        max_date = sys.float_info.max
        exposure = 0
        for entry in provenance.inputs:
            ip_obs_id, ip_product_id = mc.CaomName.decompose_provenance_input(
                entry.uri)
            logging.info(f'Retrieving provenance metadata for {ip_obs_id}.')
            ip_obs = mc.repo_get(client, 'CFHT', ip_obs_id, metrics)
            if ip_obs is not None:
                ip_plane = ip_obs.planes.get(ip_product_id)
                if (ip_plane is not None and ip_plane.time is not None and
                        ip_plane.time.bounds is not None):
                    bounds.samples.append(CoordRange1D(
                        RefCoord(pix=0.5, val=ip_plane.time.bounds.lower),
                        RefCoord(pix=1.5, val=ip_plane.time.bounds.upper)))
                    min_date = min(ip_plane.time.bounds.lower, min_date)
                    max_date = max(ip_plane.time.bounds.upper, max_date)
                    exposure += ip_plane.time.exposure
        axis = Axis(ctype='TIME', cunit='d')
        time_axis = CoordAxis1D(axis=axis,
                                error=None,
                                range=None,
                                bounds=bounds,
                                function=None)
        temporal_wcs = TemporalWCS(axis=time_axis, timesys=None, trefpos=None,
                                   mjdref=None, exposure=mc.to_float(exposure),
                                   resolution=None)
        chunk.time = temporal_wcs
    logging.debug(f'End _update_ngvs_time.')


def _update_release_date(plane, max_meta_release, headers):
    logging.debug(f'Begin _update_release_date for {plane.product_id}')
    if plane.meta_release is None:
        plane.meta_release = mc.make_time(_get_keyword(headers, 'DATE'))
        if plane.meta_release is None:
            plane.meta_release = mc.make_time(
                _get_keyword(headers, 'REL_DATE'))

    if plane.meta_release is not None:
        max_meta_release = max(max_meta_release, plane.meta_release)

    if plane.data_release is None and plane.meta_release is not None:
        plane.data_release = plane.meta_release
    logging.debug('End _update_release_date')
    return max_meta_release


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
        accumulate_bp(blueprint, uri)
        blueprints[uri] = blueprint
    return blueprints


def _filter_args(args):
    uris_for_later = []
    result = []
    if args.lineage:
        for ii in args.lineage:
            uri = ii.split('/', 1)[1]
            result.append(uri)
            storage_name = sn.get_storage_name(uri, uri)
            if not storage_name.use_metadata:
                uris_for_later.append(uri)
    else:
        raise mc.CadcException(
            f'Could not define uri from these args {args}')
    return result, uris_for_later


def to_caom2():
    parser = get_gen_proc_arg_parser()
    args = parser.parse_args()
    # set the arguments for those files that, despite being fits files,
    # are processed with a generic parser - do this because fits2caom2
    # manages parser creation based on file names, mostly
    #
    uris, generic_uris = _filter_args(args)
    blueprints = _build_blueprints(uris)
    if len(generic_uris) > 0:
        sys.argv.append('--use_generic_parser')
        for ii in generic_uris:
            sys.argv.append(ii)
        args = parser.parse_args()
    result = gen_proc(args, blueprints)
    logging.debug(f'Done {APPLICATION} processing.')
    return result
           

def cfht_proc_main_app():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        result = to_caom2()
        sys.exit(result)
    except Exception as e:
        logging.error(f'Failed {APPLICATION} execution for {args}.')
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
