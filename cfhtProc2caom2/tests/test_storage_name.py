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

import sys

from mock import patch

from cfhtProc2caom2 import main_app, storage_names

LOOKUP = {
      'MegaPipe.358.122': ['MegaPipe.358.122.G.MP9401.fits.header',
                           'MegaPipe.358.122.GRI.MP9605.fits.header',
                           'MegaPipe.358.122.I.MP9702.fits.header',
                           'MegaPipe.358.122.R.MP9601.fits.header',
                           'MegaPipe.358.122.G.MP9401.fits.gif',
                           'MegaPipe.358.122.G.MP9401.weight.fits.header',
                           'MegaPipe.358.122.GRI.MP9605.fits.gif',
                           'MegaPipe.358.122.GRI.MP9605.weight.fits.header',
                           'MegaPipe.358.122.I.MP9702.fits.gif',
                           'MegaPipe.358.122.I.MP9702.weight.fits.header',
                           'MegaPipe.358.122.R.MP9601.fits.gif',
                           'MegaPipe.358.122.R.MP9601.weight.fits.header'],
      'MegaPipe.209.126': ['MegaPipe.209.126.R.MP9602.fits.header'],
      'NGVS+0+0': ['NGVS+0+0.l.i.Mg002.fits.header',
                   'NGVS+0+0.l.i.Mg002.cat',
                   'NGVS+0+0.l.i.Mg002.fits.mask.rd.reg',
                   'NGVS+0+0.l.i.Mg002.flag.fits.fz',
                   'NGVS+0+0.l.i.Mg002.sig.fits.header',
                   'NGVS+0+0.l.i.Mg002.weight.fits.fz.header',
                   'vos:ngvs/masks/NGVS+0+0.l.i.Mg002.flag.fits.fz']}


def test_single():
    test_entry = 'MegaPipe.358.122.G.MP9401.fits'
    test_subject = storage_names.get_storage_name(test_entry, test_entry)
    assert test_subject.obs_id == 'MegaPipe.358.122', 'wrong obs id'
    assert test_subject.product_id == 'MegaPipe.358.122.G.MP9401', \
        'wrong product id'
    assert test_subject.filter_name == 'G', 'wrong filter name'


def test_is_valid():
    for key, value in LOOKUP.items():
        for entry in value:
            sn = storage_names.get_storage_name(entry, entry)
            assert sn.is_valid()
            assert sn.obs_id == key, f'wrong obs id {sn.obs_id}'


@patch('cfhtProc2caom2.main_app.gen_proc')
def test_build_uris(gen_proc_mock):
    test_obs_id = 'NGVS+0+0'
    test_lineage = get_lineage(test_obs_id)
    test_name = 'consistent_local_lineage'
    sys.argv = (f'{main_app.APPLICATION} '
                f'--observation COLLECTION {test_name} --lineage '
                f'{test_lineage}').split()
    main_app.to_caom2()
    assert gen_proc_mock.called, 'should be called'
    args, kwargs = gen_proc_mock.call_args
    generic_parser = vars(args[0]).get('use_generic_parser')
    assert generic_parser is not None, 'expect a generic_parser'
    generic_parser_text = ' '.join(ii for ii in generic_parser)
    assert LOOKUP[test_obs_id][0] not in generic_parser_text, 'expect product'
    assert LOOKUP[test_obs_id][5] not in generic_parser_text, 'expect weight'
    assert LOOKUP[test_obs_id][4] not in generic_parser_text, 'expect sig'
    assert LOOKUP[test_obs_id][1] not in generic_parser_text, 'expect cat'
    assert LOOKUP[test_obs_id][2] in generic_parser_text, 'no mask'
    assert LOOKUP[test_obs_id][3] in generic_parser_text, 'no flag'
    assert LOOKUP[test_obs_id][6] in generic_parser_text, 'no vos flag'


def get_lineage(obs_id):
    result = ''
    for ii in LOOKUP[obs_id]:
        storage_name = storage_names.get_storage_name(ii, ii)
        result = f'{result} {storage_name.lineage}'
    result = result.replace('.header', '')
    return result
