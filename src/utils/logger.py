from prody import LOGGER
import time
import os 
import logging 
from prody.utilities.logger import now, LOGGING_LEVELS

# Add a _times dict if not already present
if not hasattr(LOGGER, '_times'):
    LOGGER._times = {}

if not hasattr(LOGGER, '_reports'):
    LOGGER._reports = {}

if not hasattr(LOGGER, '_report_times'):
    LOGGER._report_times = {}

# Define your replacement method
def custom_report(self, msg='Completed in %.2fs.', label=None):
    if label not in self._times:
        self.warning(f"No timing info for label '{label}'")
        return
    elapsed = time.time() - self._times[label]
    self.debug(msg % elapsed)

    if label not in self._reports:
        self._reports[label] = elapsed
        self._report_times[label] = 1
    else:
        self._reports[label] += elapsed
        self._report_times[label] += 1
        
def start(self, filename, **kwargs):
    filename = str(filename)
    if os.path.splitext(filename)[1] == '':
        filename += '.log'
    rollover = False
    if os.path.isfile(filename) and kwargs.get('mode', None) != 'a':
        rollover = True

    logfile = logging.handlers.RotatingFileHandler(
        filename,
        mode=kwargs.get('mode', 'a'),
        maxBytes=0,
        backupCount=kwargs.get('backupcount', 1)
    )
    logfile.setLevel(LOGGING_LEVELS[kwargs.get('loglevel', 'debug')])
    logfile.setFormatter(logging.Formatter('%(message)s'))

    # Attach to the root logger so all loggers propagate here
    root_logger = logging.getLogger()
    root_logger.addHandler(logfile)
    root_logger.setLevel(logging.INFO)

    self._logger.info(f"Logging into file: {filename}")
    if rollover:
        logfile.doRollover()
    self._logger.info(f"Logging started at {str(now())}")

def close(self, filename):
    filename = str(filename)
    if os.path.splitext(filename)[1] == '':
        filename += '.log'
    root_logger = logging.getLogger()
    for index, handler in enumerate(root_logger.handlers):
        if isinstance(handler, logging.handlers.RotatingFileHandler):
            if handler.stream.name in (filename, os.path.abspath(filename)):
                self.info("Logging stopped at {0}".format(str(now())))
                handler.close()
                root_logger.removeHandler(handler)
                self.info("Closing logfile: {0}".format(filename))
                return
    self.warning("Logfile '{0}' was not found.".format(filename))

# Monkey-patch: bind the method to LOGGER
import types
LOGGER.report = types.MethodType(custom_report, LOGGER)
LOGGER.start = types.MethodType(start, LOGGER)
LOGGER.close = types.MethodType(close, LOGGER)

