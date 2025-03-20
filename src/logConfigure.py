import logging
import sys

logging.basicConfig(level=logging.WARNING, stream=sys.stdout,
                    format='%(asctime)s - %(name)s - %(message)s')

def getLogger(name):
    return logging.getLogger(name)

def setLogLevel(level):
    logLevel = getattr(logging, level, logging.WARNING)
    logging.getLogger().setLevel(logLevel)