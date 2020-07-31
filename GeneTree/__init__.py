#!/usr/bin/env python3
# coding: utf-8

__version__='0.0.1'

symb = u"\u26A0"

import sys
sys.dont_write_bytecode = True
import logging
logstream = sys.stdout
logfile = False
logger = None

from GeneTree.timer import Timer
from GeneTree import create_tree
from GeneTree import matrice
from GeneTree import clean_name
from GeneTree import cli

class MyFormatter(logging.Formatter):
    info_fmt = '> %(msg)s'
    dbg_fmt = 'DEBUG: %(module)s: %(lineno)d: %(msg)s'
    crit_fmt = '# CRITICAL: %(msg)s'
    warn_fmt = f'{symb} Warning: %(msg)s'

    def __init__(self, fmt='%(levelno)s: %(msg)s'):
        super().__init__(fmt, datefmt=None, style='%')

    def format(self, record):
        orig = self._style._fmt
        if record.levelno == logging.DEBUG:
            self._style._fmt = MyFormatter.dbg_fmt
        elif record.levelno == logging.INFO:
            self._style._fmt = MyFormatter.info_fmt
        elif record.levelno == logging.CRITICAL:
            self._style._fmt = MyFormatter.crit_fmt
        elif record.levelno == logging.WARNING:
            self._style._fmt = MyFormatter.warn_fmt
        
        res = logging.Formatter.format(self, record)
        self._style._fmt = orig
        return res

def setup_logger(name, log_s, log_f):
    
    l = logging.getLogger(name)
    l.setLevel(logging.DEBUG)
    formatter = MyFormatter()
    
    if log_f:
        fileHandler = logging.FileHandler(log_f, mode="w")
        fileHandler.setFormatter(formatter)
        l.addHandler(fileHandler)

    if log_s:
        streamHandler = logging.StreamHandler(logstream)
        streamHandler.setFormatter(formatter)
        l.addHandler(streamHandler)




def info(msg: str) -> None:
    logger.info(msg)

def warn(msg: str) -> None:
    logger.warning(msg)

def crit(msg: str, code: int) -> None:
    logger.critical(msg)
    sys.exit(code)

def err(msg: str) -> None:
    logger.error(msg)



