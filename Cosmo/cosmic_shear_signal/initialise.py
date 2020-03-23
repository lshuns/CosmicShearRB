"""
.. module:: initialise
    :synopsis: intialisation
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

.. modified to cosmic shear prediction by Shun-Sheng Li
"""
from __future__ import print_function
import sys
import os

from data import Data

def initialise(paths):
    """
    Initialisation routine

    This function take the paths,
    initialise a :class:`data` instance, a cosmological code instance.

    Parameters
    ----------
        paths: dict
            paths to all the necessary files
    """

    # Report the cosmological code version.
    # Official version number
    common_file_path = os.path.join(
        paths['cosmo'], 'include', 'common.h')
    with open(common_file_path, 'r') as common_file:
        for line in common_file:
            if line.find('_VERSION_') != -1:
                version = line.split()[-1].replace('"', '')
                break
        print('with CLASS %s' % version)


    # initialise data class
    data = Data(paths)

    # Loading up the cosmological backbone. 
    # For the moment, only CLASS has been wrapped.
    classy_path = ''
    for elem in os.listdir(os.path.join(
            paths['cosmo'], "python", "build")):
        if elem.find("lib.") != -1:
            classy_path = os.path.join(
                paths['cosmo'], "python", "build", elem)
            break

    # Inserting the previously found path into the list of folders to
    # search for python modules.
    sys.path.insert(1, classy_path)
    import classy
    cosmo = classy.Class()
    print('classy location', classy.__file__)

    return cosmo, data
