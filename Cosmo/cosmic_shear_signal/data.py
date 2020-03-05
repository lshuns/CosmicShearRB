"""
.. module:: data
   :synopsis: Define the Data and Parameter classes

.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

.. modified to cosmic shear prediction by Shun-Sheng Li
"""

from __future__ import print_function

import os


class Data(object):
    """
    Store all relevant data to communicate between the different modules.

    """

    def __init__(self, paths):
        """
        The Data class is a collection of all the necessary information. 

        It defines one useful **methods** which is called
        just once at beginning:

        * :func:`read_file`

        It has a number of different **attributes**, where contain all the input:

        # :attr:`paths`
        * :attr:`cosmo_arguments`
        * :attr:`nuisance_parameters`
        * :attr:`const`
        * :attr:`conf`

        To create an instance of this class, one must feed the following
        parameters and keyword arguments:

        Parameters
        ----------
        paths : dict
            Contains a dictionary of important local paths. 
        """

        # Initialise the path dictionnary.
        self.paths = paths


        self.cosmo_arguments = {}
        """
        Simple dictionary that will serve as a communication interface with the
        cosmological code. It contains all the parameters for the code that
        will not be set to their default values.

        :rtype:   dict
        """

        self.nuisance_parameters = {}
        """
        Nuisance paramters needed by cosmic shear prediction.
        
        :rtype: dict
        """

        self.const = {}
        """
        Hardly changed parameters 

        :rtype: dict
        """

        self.conf = {}
        """
        Parameters related to configuration
        
        :rtype: dict
        """


    def read_file(self, param, structure, field='', separate=False):
        """
        Execute all lines concerning the Data class from a parameter file

        All lines starting with `data.` will be replaced by `self.`, so the
        current instance of the class will contain all the information.

        .. note::

            A rstrip() was added at the end, because of an incomprehensible bug
            on some systems that imagined some inexistent characters at the end
            of the line... Now should work

        .. note::

            A security should be added to protect from obvious attacks.

        Parameters
        ----------
        param : str
            Name of the parameter file
        structure : str
            Name of the class entries we want to execute (mainly, data, or any
            other likelihood)

        Keyword Arguments
        -----------------
        field : str
            If nothing is specified, this routine will execute all the lines
            corresponding to the `structure` parameters. If you specify a
            specific field, like `path`, only this field will be read and
            executed.
        separate : bool
            If this flag is set to True, a container class will be created for
            the structure field, so instead of appending to the namespace of
            the data instance, it will append to a sub-namespace named in the
            same way that the desired structure. This is used to extract custom
            values from the likelihoods, allowing to specify values for the
            likelihood directly in the parameter file.

        """
        if separate:
            exec("self.%s = Container()" % structure)

        param_path = os.path.join(self.paths['param'], param)
        with open(param_path, 'r') as param_file:
            for line in param_file:
                if line.find('#') == -1 and line:
                    lhs = line.split('=')[0]
                    if lhs.find(structure+'.') != -1:
                        if field:
                            # If field is not an empty string, you want to skip
                            # the execution of the line (exec statement) if you
                            # do not find the exact searched field
                            if lhs.find('.'.join([structure, field])) == -1:
                                continue
                        if not separate:
                            exec(line.replace(structure+'.', 'self.').rstrip())
                        else:
                            exec(line.replace(
                                structure+'.', 'self.'+structure+'.').rstrip())

class Container(object):
    """Dummy class to act as a namespace for data"""
    pass
