"""
A simple and short Redis performance test.
"""
__author__ = 'Dan Gunter <dkgunter@lbl.gov>'
__date__ = '8/8/16'

import argparse
import logging
import os
import redis
import subprocess
import sys
import time

_log = logging.getLogger(__name__)
_h = logging.StreamHandler()
_h.setFormatter(logging.Formatter('%(asctime)s %(levelname)10s - %(message)s'))
_log.addHandler(_h)

def run_server(binpath=None):
    _log.info("Run Redis server")
    if binpath:
        server_cmd = os.path.join(binpath, 'redis-server')
    else:
        server_cmd = 'redis-server'
    retcode = subprocess.Popen([server_cmd])
    return retcode

def run_performance_test(cmd='set', num_items=1000, list_len=5):
    _log.info("Run Performance test")
    r = redis.StrictRedis(host='localhost', port=6379, db=0)
    data = ['bar'] * list_len
    if cmd == 'set':
        t0 = time.time()
        t1 = redis_set(r, num_items, data)
    elif cmd == 'get':
        redis_set(r, num_items, data)
        t0 = time.time()
        t1 = redis_get(r, num_items)
    elif cmd == 'mix':
        t0 = time.time()
        t1 = redis_getset(r, num_items, data)
    else:
        _log.error('Bad command: {}'.format(cmd))
        return
    report_timing(True, cmd, t1 - t0, num_items, ['list-length'], ['{}'.format(list_len)])


def redis_set(r, num_items, data):
    i = 0
    while i < num_items:
        key = 'foo' + str(i)
        r.set(key, data)
        i += 1
    return time.time()

def redis_get(r, num_items):
    i = 0
    while i < num_items:
        key = 'foo' + str(i)
        data = r.get(key)
        i += 1
    return time.time()

def redis_getset(r, num_items, data):
    i = 0
    while i < num_items:
        key = 'foo' + str(i)
        r.set(key, data)
        data2 = r.get(key)
        i += 1
    return time.time()

def report_timing(readable, mode, dt, n, info_hdr, info):
    rate = 1. * n / dt
    gap = 1. * dt / n
    if readable:
        kvp = ', '.join(['{}={}'.format(k, v) for k, v in zip(info_hdr, info)])
        print("{}: Processed {:d} items in {:.3f} seconds: {:.1f} items/sec <-> {:.6f} seconds/item. {}"
              .format(mode, n, dt, rate, gap, kvp))
    else:
        print('blah')

def verbose_add(parser):
    """Add a verbosity argument to an ArgumentParser.
    """
    parser.add_argument('-v', '--verbose', dest='vb',
                        action='count', default=0)

def verbose_set_log(vb, log):
    """Set logging level from verbosity level.
    """
    if vb >= 2:
        log.setLevel(logging.DEBUG)
    elif vb >= 1:
        log.setLevel(logging.INFO)
    else:
        log.setLevel(logging.WARN)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode', dest='mode', help='mode: get, set, mix')
    parser.add_argument('-n', '--count', dest='count', type=int, help='iterations', default=1000)
    parser.add_argument('-z', '--length', dest='len', type=int, help='list length', default=5)
    parser.add_argument('-s', '--server', dest='server', action='store_true')
    verbose_add(parser)
    args = parser.parse_args()
    verbose_set_log(args.vb, _log)
    if args.server:
        retcode = run_server()
        _log.info("Redis server stopped")
        return retcode
    else:
        run_performance_test(cmd=args.mode, num_items=args.count, list_len=args.len)

if __name__ == '__main__':
    sys.exit(main())
