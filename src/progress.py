"""Small progress protocol shared by the CLI, data loaders and Tk launcher.

Percent is work completed within the named stage, never a guessed percentage
of a first installation. ETA is estimated only after measured work advances.
This module depends only on the Python standard library.
"""
import json
import math
import time

PREFIX = '@@PROGRESS '


def parse_event(line):
    if not line.startswith(PREFIX):
        return None
    try:
        event = json.loads(line[len(PREFIX):])
        if event.get('protocol') != 'kinetic-progress-1' or not isinstance(event.get('stage'), str):
            return None
        for key in ('fraction', 'eta_s', 'elapsed_s'):
            value = event.get(key)
            if value is not None and (not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0):
                return None
        if event.get('fraction') is not None and event['fraction'] > 1:
            return None
        return event
    except (ValueError, TypeError, AttributeError):
        return None


def duration(seconds):
    if seconds is None:
        return 'оценивается'
    seconds = max(0, int(round(seconds)))
    if seconds < 60:
        return f'{seconds} с'
    if seconds < 3600:
        return f'{seconds // 60} мин {seconds % 60:02d} с'
    return f'{seconds // 3600} ч {(seconds % 3600) // 60:02d} мин'


def describe(event):
    fraction = event.get('fraction')
    percent = '' if fraction is None else f' · {100 * fraction:.1f}%'
    detail = f" · {event['detail']}" if event.get('detail') else ''
    eta = ''
    if event.get('status') == 'running':
        eta = f" · до конца этапа ≈ {duration(event.get('eta_s'))}" if event.get('eta_s') is not None else ' · время оценивается'
    return f"{event['stage']}{percent} · {event.get('module', '')}{detail}{eta}"


class Progress:
    def __init__(self, mode='off', sink=None, interval_s=.4, now=time.monotonic):
        if mode not in ('off', 'text', 'json'):
            raise ValueError('progress: off, text или json')
        self.mode = mode
        self.sink = sink or (lambda line: print(line, flush=True))
        self.interval = interval_s
        self.now = now
        self.stage = ''
        self.module_name = ''
        self.detail = ''
        self.fraction = None
        self.started = now()
        self.last_emit = -math.inf
        self.rate = None
        self.last_work = 0.
        self.last_work_time = self.started
        self.advances = 0

    def start(self, stage, module='', detail=''):
        self.stage, self.module_name, self.detail = stage, module, detail
        self.started = self.now()
        self.fraction = None
        self.rate = None
        self.last_work = 0.
        self.last_work_time = self.started
        self.advances = 0
        self.emit(force=True)

    def module(self, module, detail=None, force=False):
        self.module_name = module
        if detail is not None:
            self.detail = detail
        self.emit(force=force)

    def update(self, completed, total, detail=None, force=False):
        if not math.isfinite(total) or total <= 0 or not math.isfinite(completed) or completed < 0:
            raise ValueError('Неверный объём работы progress')
        fraction = min(1., completed / total)
        if self.fraction is not None and fraction < self.fraction - 1e-14:
            raise ValueError('Прогресс этапа не должен уменьшаться')
        current = self.now()
        advance = fraction - self.last_work
        # A sample endpoint can differ from the last accepted t by one ULP.
        # Dividing its elapsed wall time by that ~1e-17 fraction produced
        # millions of hours of ETA. Accumulate tiny advances until meaningful.
        if advance > 1e-12:
            measured = (current - self.last_work_time) / advance
            self.rate = measured if self.rate is None else .3 * measured + .7 * self.rate
            self.last_work, self.last_work_time = fraction, current
            self.advances += 1
        self.fraction = fraction
        if detail is not None:
            self.detail = detail
        self.emit(force=force)

    def emit(self, force=False, status='running'):
        if self.mode == 'off':
            return
        current = self.now()
        if not force and current - self.last_emit < self.interval:
            return
        self.last_emit = current
        remaining = None
        if self.fraction is not None and self.advances >= 3 and self.rate is not None:
            # A long operation/rejected step increases this lower estimate too.
            rate = max(self.rate, (current - self.started) / max(self.fraction, 1e-300))
            remaining = max(0., (1 - self.fraction) * rate)
        event = dict(protocol='kinetic-progress-1', stage=self.stage,
                     module=self.module_name, fraction=self.fraction,
                     eta_s=remaining, elapsed_s=max(0., current - self.started),
                     detail=self.detail, status=status)
        self.sink(PREFIX + json.dumps(event, ensure_ascii=False, allow_nan=False)
                  if self.mode == 'json' else describe(event))

    def done(self, detail=''):
        self.fraction = 1.
        if detail:
            self.detail = detail
        self.emit(force=True, status='done')

    def fail(self, detail):
        self.detail = str(detail)
        self.emit(force=True, status='failed')


def progress_from_config(config, mode=None):
    options = config.get('progress', {})
    if not isinstance(options, dict) or set(options) - {'mode', 'interval_s'}:
        raise ValueError('progress: mode и interval_s')
    interval = options.get('interval_s', .4)
    if isinstance(interval, bool) or not isinstance(interval, (int, float)) or not math.isfinite(interval) or interval < .05:
        raise ValueError('progress.interval_s >= 0.05')
    return Progress(mode or options.get('mode', 'text'), interval_s=interval)
