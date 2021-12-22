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

import logging
import sys

from caom2 import Observation, CalibrationLevel, DataProductType, ProductType
from caom2 import TemporalWCS, CoordAxis1D, Axis, CoordBounds1D, CoordRange1D
from caom2 import RefCoord
from caom2utils import update_artifact_meta
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import client_composable as clc
from caom2pipe import manage_composable as mc
from cfhtProc2caom2 import storage_names as sn


__all__ = ['APPLICATION', 'CFHTProductMapping']


APPLICATION = 'cfhtProc2caom2'
INSTRUMENT = 'MegaPrime'

# filter names have to conform to those of the SVO filter service, because
# they're going to be used to create a URL to query from there
FILTER_REPAIR_LOOKUP = {
    'i': 'i_sdss',  # i.MP9703
    'I2': 'i',  # i.MP9702
    'I': 'i1',  # i.MP9701
    'g': 'g_sdss',  # g.MP9402
    'G': 'g',  # g.MP9401
    'GRI': 'gri',  # gri.MP9605
    'r': 'r_sdss',  # r.MP9602
    'R': 'r',  # r.MP9601
    'u': 'u_sdss',  # u.MP9302
    'U': 'u',  # u.MP9301
    'z': 'z_sdss',  # z.MP9901
    'Z': 'z',
}  # z.MP9801
filter_cache = ac.FilterMetadataCache(
    FILTER_REPAIR_LOOKUP, {}, 'CFHT', {}, 'NONE'
)

# because 1 is never enough ... and also because, in an effort to keep
# the number of filter names smaller, they have to have the case of
# CFHT, which is filter name letters are lower-case, but CFHT parts
# numbers are upper-case.
CAOM_FILTER_REPAIR_LOOKUP = {
    'I.MP9703': 'i.MP9703',
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
    'Z.MP9801': 'z.MP9801',
}


class CFHTProductMapping(cc.TelescopeMapping):
    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def accumulate_blueprint(self, bp, application=None):
        super().accumulate_blueprint(bp, APPLICATION)
        bp.configure_position_axes((1, 2))
        bp.set('DerivedObservation.members', {})
        bp.set('Observation.type', 'OBJECT')

        bp.set('Observation.proposal.id', 'get_proposal_id()')

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

        bp.set('Plane.calibrationLevel', 'get_calibration_level()')
        bp.set('Plane.dataProductType', 'get_data_product_type()')
        bp.set('Plane.provenance.producer', 'CADC')

        bp.set('Artifact.productType', 'get_artifact_product_type()')
        bp.set('Artifact.releaseType', 'data')

    def update(self, observation, file_info, caom_repo_client):
        """
        Called to fill multiple CAOM model elements and/or attributes.
        """
        mc.check_param(observation, Observation)
        self._logger.debug(
            f'Update for {observation.observation_id} with '
            f'{self._storage_name.file_name}.'
        )

        max_meta_release = observation.meta_release
        min_seeing = None
        if (
            observation.environment is not None
            and observation.environment.seeing is not None
        ):
            min_seeing = observation.environment.seeing
        for plane in observation.planes.values():
            max_meta_release = self._update_release_date(
                plane, max_meta_release
            )
            if plane.product_id != self._storage_name.product_id:
                continue
            min_seeing = self._minimize(
                min_seeing, self._get_keyword('FINALIQ')
            )
            for artifact in plane.artifacts.values():
                if artifact.uri != self._storage_name.file_uri:
                    self._logger.debug(
                        f'Artifact uri is {artifact.uri} but working on '
                        f'{self._storage_name.file_uri}. Continuing.'
                    )
                    continue
                update_artifact_meta(artifact, file_info)

                if not self._storage_name.is_catalog:
                    if self._update_parts(artifact):
                        continue
                    for part in artifact.parts.values():
                        for chunk in part.chunks:
                            if self._informative_uri():
                                self._update_energy(
                                    chunk, observation.observation_id
                                )
                            self._update_position(chunk)
                            # SGw - 24-11-20 - set observable/axis to None:
                            # there is no such information in the files
                            cc.reset_observable(chunk)

                if self._informative_uri() and plane.provenance is not None:
                    # SGw - 22-01-21
                    # When re-processing, I sometimes find that an image
                    # wasn't as well calibrated as I thought it was and it
                    # gets removed from the next generation data products. So
                    # I would like the ability to remove inputs
                    cc.update_plane_provenance_single(
                        plane,
                        self._headers,
                        'HISTORY',
                        'CFHT',
                        _repair_history_provenance_value,
                        observation.observation_id,
                    )
                    if plane.provenance.run_id == 'None':
                        plane.provenance.run_id = None
                    if (
                        plane.provenance.keywords is not None
                        and 'None' in plane.provenance.keywords
                    ):
                        plane.provenance.keywords.remove('None')

            # _update_time is dependent on provenance information
            # that is generated right before this
            self._update_time(
                plane, observation.observation_id, caom_repo_client
            )
        observation.meta_release = max_meta_release
        if observation.environment is not None:
            observation.environment.seeing = min_seeing
        cc.update_observation_members(observation)
        self._logger.debug('Done update.')
        return observation

    def _get_keyword(self, keyword):
        result = None
        if len(self._headers) > 0:
            result = self._headers[0].get(keyword)
            if result is None and len(self._headers) > 1:
                result = self._headers[1].get(keyword)
        return result

    def get_artifact_product_type(self, ext):
        result = ProductType.SCIENCE
        if (
            'weight' in self._storage_name.file_uri
            or '.sig' in self._storage_name.file_uri
        ):
            result = ProductType.WEIGHT
        elif (
            'mask.rd.reg' in self._storage_name.file_uri
            or '.flag' in self._storage_name.file_uri
        ):
            result = ProductType.AUXILIARY
        return result

    def get_calibration_level(self, ext):
        result = CalibrationLevel.PRODUCT
        if self._storage_name.is_catalog:
            result = CalibrationLevel.ANALYSIS_PRODUCT
        return result

    def get_data_product_type(self, ext):
        result = DataProductType.IMAGE
        if self._storage_name.is_catalog:
            # PD 09-12-20
            # I need to modify the ObsCore view to filter the observations
            # with DataProductType.catalog out (not compliant to spec) but
            # there is a different value measurements that means roughly the
            # same thing. There should be a DataProductType constant declared
            # in the py library for this value.
            result = DataProductType.MEASUREMENTS
        return result

    def get_proposal_id(self, ext):
        # JJK - 23-03-20
        # Replace with value of ProposalID from the first input ‘member’
        # from list
        # TODO
        return None

    def get_provenance_version(self, ext):
        return self._storage_name.version

    def get_target_name(self, ext):
        return self._storage_name.file_name.split('.')[0]

    def _informative_uri(self):
        result = False
        if (
            '.weight' not in self._storage_name.file_uri
            and '.sig' not in self._storage_name.file_uri
            and '.cat' not in self._storage_name.file_uri
            and '.mask' not in self._storage_name.file_uri
            and '.flag' not in self._storage_name.file_uri
            and '.gif' not in self._storage_name.file_uri
        ):
            # all the excluded names have fewer useful keywords
            result = True
        return result

    def _minimize(self, x, candidate):
        result = x
        if candidate is not None:
            if x is None:
                result = candidate
            else:
                result = min(x, candidate)
        return result

    def _update_energy(self, chunk, obs_id):
        self._logger.debug(f'Begin _update_energy for {obs_id}')
        filter_md = filter_cache.get_svo_filter(
            INSTRUMENT, self._storage_name.filter_name
        )
        cc.build_chunk_energy_range(
            chunk, self._storage_name.filter_name, filter_md
        )
        self._logger.debug(f'End _update_energy.')

    def _update_parts(self, artifact):
        return False

    def _update_position(self, chunk):
        pass

    def _update_release_date(self, plane, max_meta_release):
        self._logger.debug(
            f'Begin _update_release_date for {plane.product_id}'
        )
        if plane.meta_release is None:
            plane.meta_release = mc.make_time(self._get_keyword('DATE'))
            if plane.meta_release is None:
                plane.meta_release = mc.make_time(
                    self._get_keyword('REL_DATE')
                )

        if plane.meta_release is not None:
            max_meta_release = max(max_meta_release, plane.meta_release)

        if plane.data_release is None and plane.meta_release is not None:
            plane.data_release = plane.meta_release
        self._logger.debug('End _update_release_date')
        return max_meta_release

    def _update_time(self, plane, obs_id, caom_repo_client):
        pass


class NGVSProductMapping(CFHTProductMapping):
    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def accumulate_blueprint(self, bp, application=None):
        super().accumulate_blueprint(bp, APPLICATION)
        bp.clear('Observation.algorithm.name')

        # NGVS
        # JJK - comments by email - 28-03-20
        # if _informative_uri(uri):
        # make sure information set from header keywords is only set for
        # the fits files where it's accessible
        bp.set('Observation.metaRelease', '2022-01-01 00:00:00')

        bp.set('Observation.proposal.pi', 'Laura Ferrarese')
        bp.set('Observation.proposal.project', 'NGVS')
        bp.set(
            'Observation.proposal.title',
            'Next Generation Virgo Cluster Survey',
        )
        bp.set('Observation.proposal.keywords', 'Galaxy Cluster Dwarfs')

        bp.set('Observation.target.name', 'get_target_name()')
        if 'weight' not in self._storage_name.file_uri:
            bp.set('Observation.target.type', 'field')

        bp.set('Plane.dataRelease', '2022-01-01T00:00:00')
        bp.set('Plane.metaRelease', '2022-01-01T00:00:00')
        # SGw 28-07-20
        # The last 4 characters in the ID column (l128 vs g002 vs g004)
        # refer to variation in the pipeline. If they could be captured as
        # plane.provenance.name = MegaPipe_g004 that would be great.
        bp.set(
            'Plane.provenance.name', f'MegaPipe_{self._storage_name.version}'
        )
        bp.clear('Plane.provenance.keywords')
        bp.add_fits_attribute('Plane.provenance.keywords', 'COMBINET')
        bp.set('Plane.provenance.project', 'NGVS')
        bp.set(
            'Plane.provenance.reference',
            'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/ngvs/'
            'docs/ngvsdoc.html',
        )
        bp.set('Plane.provenance.version', 'get_provenance_version()')
        bp.clear('Plane.metrics.magLimit')
        bp.add_fits_attribute('Plane.metrics.magLimit', 'MAGLIM')

        bp.set('Artifact.productType', 'get_artifact_product_type()')

        bp.set('Chunk.energy.bandpassName', 'get_bandpass_name()')

        self._logger.debug(f'End accumulate_bp.')

    def get_bandpass_name(self, ext):
        reverse_filter_lookup = {
            'i': 'i.MP9703',
            'g': 'g.MP9402',
            'r': 'r.MP9602',
            'u': 'u.MP9302',
            'z': 'z.MP9901',
        }
        result = None
        if self._storage_name.filter_name is not None:
            result = reverse_filter_lookup.get(self._storage_name.filter_name)
        return result

    def _update_time(self, plane, obs_id, caom_repo_client):
        self._logger.debug(f'Begin _update_time for {obs_id}')
        if (
            not self._storage_name.is_catalog
            and plane.provenance is not None
            and len(plane.provenance.inputs) > 0
        ):
            for artifact in plane.artifacts.values():
                if artifact.uri != self._storage_name.file_uri:
                    continue
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        config = mc.Config()
                        metrics = mc.Metrics(config)
                        bounds = CoordBounds1D()
                        min_date = 0
                        max_date = sys.float_info.max
                        exposure = 0
                        for entry in plane.provenance.inputs:
                            (
                                ip_obs_id,
                                ip_product_id,
                            ) = mc.CaomName.decompose_provenance_input(
                                entry.uri,
                            )
                            self._logger.info(
                                f'Retrieving provenance metadata for '
                                f'{ip_obs_id}.'
                            )
                            ip_obs = clc.repo_get(
                                caom_repo_client, 'CFHT', ip_obs_id, metrics
                            )
                            if ip_obs is not None:
                                ip_plane = ip_obs.planes.get(ip_product_id)
                                if (
                                    ip_plane is not None
                                    and ip_plane.time is not None
                                    and ip_plane.time.bounds is not None
                                ):
                                    bounds.samples.append(
                                        CoordRange1D(
                                            RefCoord(
                                                pix=0.5,
                                                val=ip_plane.time.bounds.lower,
                                            ),
                                            RefCoord(
                                                pix=1.5,
                                                val=ip_plane.time.bounds.upper,
                                            ),
                                        ),
                                    )
                                    min_date = min(
                                        ip_plane.time.bounds.lower, min_date
                                    )
                                    max_date = max(
                                        ip_plane.time.bounds.upper, max_date
                                    )
                                    exposure += ip_plane.time.exposure
                        axis = Axis(ctype='TIME', cunit='d')
                        time_axis = CoordAxis1D(
                            axis=axis,
                            error=None,
                            range=None,
                            bounds=bounds,
                            function=None,
                        )
                        temporal_wcs = TemporalWCS(
                            axis=time_axis,
                            timesys=None,
                            trefpos=None,
                            mjdref=None,
                            exposure=mc.to_float(exposure),
                            resolution=None,
                        )
                        chunk.time = temporal_wcs
        self._logger.debug(f'End _update_time.')


class MegapipeProductMapping(CFHTProductMapping):
    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def accumulate_blueprint(self, bp, application=None):
        super().accumulate_blueprint(bp, APPLICATION)
        self._logger.debug(
            f'Begin _accumulate_mp_bp for {self._storage_name.file_name}'
        )
        bp.clear('Observation.metaRelease')
        bp.add_fits_attribute('Observation.metaRelease', 'DATE')
        bp.add_fits_attribute('Observation.metaRelease', 'REL_DATE')

        bp.set('Observation.algorithm.name', 'MEGAPIPE')

        bp.set('Observation.environment.photometric', True)
        bp.clear('Observation.environment.seeing')
        bp.add_fits_attribute('Observation.environment.seeing', 'FINALIQ')
        bp.set('Observation.target.standard', False)

        bp.clear('Plane.metaRelease')
        bp.add_fits_attribute('Plane.metaRelease', 'DATE')
        bp.add_fits_attribute('Plane.metaRelease', 'REL_DATE')
        bp.clear('Plane.dataRelease')
        bp.add_fits_attribute('Plane.dataRelease', 'DATE')
        bp.add_fits_attribute('Plane.dataRelease', 'REL_DATE')

        bp.set('Plane.provenance.name', 'MEGAPIPE')
        bp.set('Plane.provenance.producer', 'CADC')
        bp.set('Plane.provenance.version', '2.0')
        bp.set(
            'Plane.provenance.reference',
            'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/',
        )

        if self._storage_name.is_preview:
            bp.set('Artifact.productType', ProductType.PREVIEW)

        bp.clear('Chunk.position.resolution')
        bp.add_fits_attribute('Chunk.position.resolution', 'FINALIQ')

        bp.clear('Chunk.position.coordsys')
        bp.add_fits_attribute('Chunk.position.coordsys', 'RADECSYS')

        bp.clear('Observation.target.name')
        if not self._storage_name.is_weight:
            bp.clear('Plane.provenance.lastExecuted')
            bp.add_fits_attribute('Plane.provenance.lastExecuted', 'DATE')
            bp.add_fits_attribute('Observation.target.name', 'OBJECT')
            bp.clear('Plane.metrics.magLimit')
            bp.add_fits_attribute('Plane.metrics.magLimit', 'ML_5SIGA')

        self._logger.debug(f'End accumulate_blueprint.')

    def _update_energy(self, chunk, obs_id):
        super()._update_energy(chunk, obs_id)
        temp = self._get_keyword('FILTER')
        # an attempt to keep the number of unique filter names lower
        chunk.energy.bandpass_name = CAOM_FILTER_REPAIR_LOOKUP.get(temp)
        chunk.energy.resolving_power = None
        self._logger.debug(f'End _update_energy.')

    def _update_parts(self, artifact):
        if artifact.product_type is ProductType.WEIGHT:
            artifact.parts = None
            return True
        return False

    def _update_position(self, chunk):
        if chunk.position is not None:
            chunk.position.resolution = None


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
