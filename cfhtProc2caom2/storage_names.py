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
#  : 4 $
#
# ***********************************************************************
#

import logging
from caom2pipe import manage_composable as mc


__all__ = ['get_storage_name', 'is_ngvs', 'MEGAPIPEName', 'MP_ARCHIVE',
           'MP_COLLECTION', 'NGVSName', 'NGVS_ARCHIVE', 'NGVS_COLLECTION']

NGVS_COLLECTION = 'NGVS'
NGVS_ARCHIVE = 'NGVS'
MP_COLLECTION = 'CFHTMEGAPIPE'
MP_ARCHIVE = 'CFHTSG'


def get_storage_name(file_name, entry):
    """
    :param file_name may be a file name or a URI. The parameter has to be
        called 'file_name', because there's a constructor call that
        expects that.
    :param entry str what came from the todo.txt file
    :return: The StorageName extension class for the entry
    """
    uri = None
    temp = file_name
    scheme = 'ad'
    if '/' in file_name:
        # parameter is a URI, adjust accordingly
        scheme, archive, f_name = mc.decompose_uri(file_name)
        temp = f_name
        uri = file_name
    if is_ngvs(temp):
        result = NGVSName(file_name=temp, scheme=scheme, entry=entry)
    else:
        result = MEGAPIPEName(file_name=temp, entry=entry)
    result.file_uri = uri
    return result


def is_ngvs(file_name):
    return 'NGVS' in file_name


class CFHTAdvancedProduct(mc.StorageName):

    NAME_PATTERN = '*'

    def __init__(self, obs_id, file_name, collection, archive, scheme, entry):
        super(CFHTAdvancedProduct, self).__init__(
            obs_id, collection=collection, archive=archive,
            collection_pattern=CFHTAdvancedProduct.NAME_PATTERN,
            fname_on_disk=file_name, compression='', scheme=scheme,
            entry=entry)
        self._file_name = file_name
        self._file_uri = None

    @property
    def file_name(self):
        return self._file_name

    @property
    def file_uri(self):
        result = self._file_uri
        if result is None:
            result = super(CFHTAdvancedProduct, self).file_uri
        return result

    @file_uri.setter
    def file_uri(self, value):
        self._file_uri = value

    @property
    def product_id(self):
        return self._product_id

    def is_valid(self):
        return True

    @property
    def is_catalog(self):
        return 'cat' in self._file_name

    @property
    def is_weight(self):
        return 'weight' in self._file_name

    @property
    def use_metadata(self):
        return not ('mask.rd.reg' in self._file_name or
                    '.flag' in self._file_name or '.gif' in self._file_name or
                    '.weight' in self._file_name)


class NGVSName(CFHTAdvancedProduct):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support uncompressed files in storage

    SGw 28-07-20
    NGVS Results Set:
    https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/gsky/
    NGVSgetsets.pl?&ra1=153.004811000824&ra2=222.395188999176&
    dec1=-9.330781193814353&dec2=32.647815298665755&F=U&F=G&F=R&F=I&F=Z&
    T=long&T=short&T=merged&P=Ml128&P=Mg002&P=Mg004&E=none&E=normal&E=LSB&
    grid=true&stacked=true&point=true&galex=false&single=true

    All the rows with the same pointing (e.g. NGVS+0+0) should be in the
    same observation. One plane per row, one artifact per link.

    File names look like:
    NGVS+0+0.l.i.Mg002.sig.fits
    NGVS+0+0.l.i.Mg002.fits
    NGVS+0+0.l.i.Mg002.weight.fits.fz
    NGVS+0+0.l.i.Mg002.fits.mask.rd.reg
    NGVS+0+0.l.i.Mg002.cat
    NGVS+0+0.l.i.Mg002.flag.fits.fz
    """

    def __init__(self, file_name, scheme, entry):
        self.fname_in_ad = file_name
        obs_id = NGVSName.get_obs_id(file_name)
        super(NGVSName, self).__init__(
            obs_id, file_name, NGVS_COLLECTION, NGVS_ARCHIVE, scheme,
            entry=entry)
        self._product_id = NGVSName.get_product_id(file_name)
        self._logger = logging.getLogger(__name__)
        self._logger.debug(self)

    @property
    def filter_name(self):
        bits = self._file_name.split('.')
        if len(bits) == 2:
            bits = self._file_name.replace('.psf', '').split('_')
        result = bits[2]
        if self._file_name.startswith('psfex'):
            result = bits[3]
        return result

    @property
    def use_metadata(self):
        return not ('mask.rd.reg' in self._file_name or
                    '.flag' in self._file_name)

    @property
    def version(self):
        bits = self._file_name.split('.')
        if len(bits) == 2:
            bits = self._file_name.replace('.psf', '').split('_')
        result = bits[3]
        if self._file_name.startswith('psfex'):
            result = bits[4]
        return result

    @staticmethod
    def get_product_id(f_name):
        bits = f_name.split('.')
        if len(bits) == 2:
            bits = f_name.replace('.psf', '').split('_')
        result = f'{bits[1]}.{bits[2]}.{bits[3]}'
        if f_name.startswith('psfex'):
            result = f'{bits[2]}.{bits[3]}.{bits[4]}'
        return result

    @staticmethod
    def get_obs_id(f_name):
        bits = f_name.split('.')
        if len(bits) == 2:
            bits = f_name.split('_')
        result = bits[0]
        if f_name.startswith('psfex'):
            result = bits[1]
        return result

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        return name.replace('.fits', '').replace('.fz', '').replace(
            '.header', '').replace('.sig', '').replace('.weight', '').replace(
            '.cat', '').replace('.mask.rd.reg', '').replace(
            '.flag', '').replace('.psf', '')

    @staticmethod
    def use_later_extensions(name):
        return '.weight' in name


class MEGAPIPEName(CFHTAdvancedProduct):
    """
    Compression is varied, so handle it on a case-by-case basis.

    This example from the story CADC-8590:
    https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ui/view?ID=
        ivo%3A%2F%2Fcadc.nrc.ca%2FCFHTMEGAPIPE%3FMegaPipe.358.122

    Is NOT like this example (from the previous excursion into building a
    data engineering pipeline for MegaPrime):
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ui/view?ID=
        ivo%3A%2F%2Fcadc.nrc.ca%2FCFHTMEGAPIPE%3FW3%2B2-3

    SGw 03-11-20
    The first observation is from MegaPipe 2.0. The naming convention changed
    between the observations. There is a possibility that I may create data
    with another naming convention, but the G024.234.... and W+1+1
    conventions are definitely a thing of the past. There probably won’t be a
    time when an observation like W3+2-3 will need to be re-created.
    """

    def __init__(self, file_name, entry):
        obs_id = MEGAPIPEName.get_obs_id(file_name)
        super(MEGAPIPEName, self).__init__(
            obs_id, file_name, MP_COLLECTION, MP_ARCHIVE, scheme='ad',
            entry=entry)
        self._product_id = MEGAPIPEName.get_product_id(file_name)
        self._logger = logging.getLogger(__name__)
        self._logger.debug(self)

    def __str__(self):
        return f'{super(MEGAPIPEName, self).__str__()} {self.filter_name}'

    @property
    def filter_name(self):
        bits = self._file_name.split('.')
        return bits[3]

    @property
    def is_preview(self):
        return '.gif' in self._file_name

    @staticmethod
    def get_obs_id(f_name):
        bits = f_name.split('.')
        return f'{bits[0]}.{bits[1]}.{bits[2]}'

    @staticmethod
    def get_product_id(f_name):
        result = MEGAPIPEName.remove_extensions(f_name)
        if 'weight' in result:
            result = result.replace('.weight', '')
        return result

    @staticmethod
    def remove_extensions(f_name):
        return f_name.replace('.fits', '').replace('.fz', '').replace(
            '.header', '').replace('.gif', '')
