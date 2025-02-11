# encoding: utf-8

import os
import logging

log = logging.getLogger()

def fileenv(varname, fallback=''):
    """
    Mimics the behaviour of the eponym bash function used, for instance, in the postgresql docker image
    Given an environment variable name (`varname`), it first adds the '_FILE' suffix and looks if it point to a file,
    containing a secret. If yes, this will be the value returned.
    If not, it looks for the environment variable's content.
    You can set a default value (`fallback`)
    :param varname:
    :param fallback:
    :return:
    """
    value=''
    filename = os.getenv('{}_FILE'.format(varname), '')
    if filename:
        # try to read from secret file
        try:
            with open(filename) as secret_file:
                value = secret_file.read().rstrip('\n')
                return value
        except IOError as e:
            log.warn("invalid path to secret file {}, {}".format(filename, e))
    else:
        value = os.getenv(varname, fallback)
    return value
