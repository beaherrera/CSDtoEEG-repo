# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 10:36:18 2021

@author: Beatriz Herrera
"""

from __future__ import division
import os
import posixpath


def get_templatename(f):
    """
    Assess from hoc file the templatename being specified within.

    Arguments
    ---------
    f : file, mode 'r'

    Returns
    -------
    templatename : str

    """
    for line in f.readlines():
        if 'begintemplate' in line.split():
            templatename = line.split()[-1]
            print('template {} found!'.format(templatename))
            continue

    return templatename


def posixpth(pth):
    """Replace Windows path separators with posix style separators."""
    return pth.replace(os.sep, posixpath.sep)
