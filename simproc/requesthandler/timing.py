"""For time stamping and time deltas"""

#Standard library
from datetime import datetime

#Standard time format
TIMEFORMAT="%a %d-%b-%Y %H:%M:%S.%f"

def format_delta(tdelta):
  "Return a string in standard format for a time delta"
  hour,rem=divmod(tdelta.seconds,60)
  min,sec=divmod(rem,60)
  return "%dD %02d:%02d:%02d (%d.%d sec)"%(tdelta.days,hour,min,sec,tdelta.total_seconds(),tdelta.microseconds)

class Timer(object):
  """A simple timer"""
  def __init__(self):
    self.start=datetime.now()
  def stop(self):
    self.end=datetime.now()
    self.delta=self.end-self.start
    self.delta_str=format_delta(self.delta)
    return self.delta_str

def timestamp(dt=None):
  if dt is None:
    dt=datetime.now()
  return dt.strftime(TIMEFORMAT)