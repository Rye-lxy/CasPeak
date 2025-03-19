import logging
import sys

logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format='%(asctime)s - %(name)s - %(message)s')

def getLogger(name):
    return logging.getLogger(name)