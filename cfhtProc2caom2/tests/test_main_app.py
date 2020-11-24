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
import os
import sys

from mock import patch

from astropy.io.votable import parse_single_table
from cfhtProc2caom2 import main_app, storage_names
from caom2pipe import manage_composable as mc

import test_storage_name

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
PLUGIN = os.path.join(os.path.dirname(THIS_DIR), 'main_app.py')


def pytest_generate_tests(metafunc):
    obs_id_list = []
    for key, value in test_storage_name.LOOKUP.items():
        obs_id_list.append(key)
    metafunc.parametrize('test_name', obs_id_list)


@patch('caom2pipe.astro_composable.get_vo_table')
@patch('caom2pipe.manage_composable.repo_get')
@patch('caom2utils.fits2caom2.CadcDataClient')
@patch('caom2utils.fits2caom2.Client')
def test_main_app(
        vo_client, data_client_mock, repo_get_mock, vo_mock, test_name):
    obs_id = os.path.basename(test_name)
    storage_name = storage_names.get_storage_name(
        test_storage_name.LOOKUP[obs_id][0],
        test_storage_name.LOOKUP[obs_id][0])
    working_dir = get_work_dir(test_name)
    output_file = f'{TEST_DATA_DIR}/{working_dir}/{obs_id}.actual.xml'
    input_file = f'{TEST_DATA_DIR}/{working_dir}/{obs_id}.in.xml'
    obs_path = f'{TEST_DATA_DIR}/{working_dir}/{obs_id}.expected.xml'
    data_client_mock.return_value.get_file_info.side_effect = get_file_info
    vo_client.return_value.get_node.side_effect = _get_node_mock
    repo_get_mock.side_effect = _repo_read_mock
    vo_mock.side_effect = _vo_mock

    sys.argv = \
        (f'{main_app.APPLICATION} --no_validate --local '
         f'{_get_local(test_name)} -i {input_file} -o {output_file} --plugin '
         f'{PLUGIN} --module {PLUGIN} --lineage '
         f'{test_storage_name.get_lineage(test_name)}').split()
    print(sys.argv)
    main_app.to_caom2()

    compare_result = mc.compare_observations(output_file, obs_path)
    if compare_result is not None:
        raise AssertionError(compare_result)
    # assert False  # cause I want to see logging messages


def get_work_dir(value):
    working_dir = 'megaprime'
    if 'NGVS' in value:
        working_dir = 'ngvs'
    return working_dir


def get_file_info(archive, file_id):
    if '.cat' in file_id:
        return {'type': 'text/plain'}
    elif '.gif' in file_id:
        return {'type': 'image/gif'}
    else:
        return {'type': 'application/fits'}


def _get_node_mock(uri, **kwargs):
    node = type('', (), {})()
    node.props = {'length': 42,
                  'MD5': '1234'}
    return node


def _get_local(obs_id):
    result = ''
    for ii in test_storage_name.LOOKUP[obs_id]:
        work_dir = get_work_dir(ii)
        result = f'{result} {TEST_DATA_DIR}/{work_dir}/{ii}'
    return result


def _repo_read_mock(ignore1, ignore2, obs_id, ignore4):
    work_dir = get_work_dir(obs_id)
    fqn = f'{TEST_DATA_DIR}/{work_dir}/{obs_id}.xml'
    return mc.read_obs_from_file(fqn)


def _vo_mock(url):
    try:
        x = url.split('/')
        filter_name = x[-1].replace('&VERB=0', '')
        votable = parse_single_table(
            f'{TEST_DATA_DIR}/{filter_name}.xml')
        return votable, None
    except Exception as e:
        logging.error(f'get_vo_table failure for url {url}')
