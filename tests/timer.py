import datetime
import time

from math import floor

class Chronograph:
    def __init__(self, title):
        self.title = title
        self.start_time = datetime.datetime.now()

    def __str__(self):
        elapsed_time = datetime.datetime.now() - self.start_time
        tot_secs = floor(elapsed_time.total_seconds())
        days = floor(tot_secs / 86400)
        hours = floor((tot_secs % 86400) / 3600)
        minutes = floor((tot_secs % 3600) / 60)
        seconds = (tot_secs % 60) + elapsed_time.total_seconds() - tot_secs
        return "\"{}\" lasted {}s ({})".format(
            self.title,
            elapsed_time.total_seconds(),
            "{:d}d {:02d}h {:02d}m {:5.03f}s".format(
                days,
                hours,
                minutes,
                seconds
        ))
