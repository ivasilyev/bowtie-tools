#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import argparse
import multiprocessing
import subprocess
import redis
import uuid
import hashlib
import time
import json


def parse_args():
    starting_parser = argparse.ArgumentParser(description="This script builds redis-based queue to launch the filtering pipeline")
    starting_parser.add_argument("-f", "--filter", required=True,
                                 help="Filtering DNA sequence REFDATA")
    starting_parser.add_argument("-c", "--coverage", required=True,
                                 help="DNA sequence REFDATA to calculate coverage")
    starting_parser.add_argument("-s", "--sampledata", required=True,
                                 help="Input list containing two tab-delimited columns for colorspace or non-colorspace sequences and three for paired-end sequences: sample name and absolute path(s). May contain a header")
    starting_parser.add_argument("-m", "--mask", default=None,
                                 help="(Optional) Mask to be added to resulting files. Automtically apended by both REFDATA file names")
    starting_parser.add_argument("-t", "--threads", default=None, type=int,
                                 help="(Optional) Number of CPU cores to use, maximal by default")
    starting_parser.add_argument("-o", "--output", required=True,
                                 help="Output directory")
    return starting_parser.parse_args()


def parse_namespace():
    namespace = parse_args()
    for file_name in [namespace.filter, namespace.coverage, namespace.sampledata]:
        if not os.path.isfile(file_name):
            raise ValueError("Not found: '" + file_name + "'\nIf you're using Docker, please make sure you have mounted required volume with the '-v' flag.")
    default_threads = int(subprocess.getoutput("nproc"))
    if not namespace.threads or default_threads < namespace.threads:
        namespace.threads = default_threads
    return namespace.filter, namespace.coverage, namespace.sampledata, namespace.mask, str(namespace.threads), namespace.output


def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError:
        pass

# Based on http://peter-hoffmann.com/2012/python-simple-queue-redis-queue.html
# and the suggestion in the redis documentation for RPOPLPUSH, at
# http://redis.io/commands/rpoplpush, which suggests how to implement a work-queue.

class RedisWQ(object):
    """Simple Finite Work Queue with Redis Backend

    This work queue is finite: as long as no more work is added
    after workers start, the workers can detect when the queue
    is completely empty.

    The items in the work queue are assumed to have unique values.

    This object is not intended to be used by multiple threads
    concurrently.
    """
    def __init__(self, name, **redis_kwargs):
       """The default connection parameters are: host='localhost', port=6379, db=0

       The work queue is identified by "name".  The library may create other
       keys with "name" as a prefix.
       """
       self._db = redis.StrictRedis(**redis_kwargs)
       # The session ID will uniquely identify this "worker".
       self._session = str(uuid.uuid4())
       # Work queue is implemented as two queues: main, and processing.
       # Work is initially in main, and moved to processing when a client picks it up.
       self._main_q_key = name
       self._processing_q_key = name + ":processing"
       self._lease_key_prefix = name + ":leased_by_session:"

    def sessionID(self):
        """Return the ID for this session."""
        return self._session

    def _main_qsize(self):
        """Return the size of the main queue."""
        return self._db.llen(self._main_q_key)

    def _processing_qsize(self):
        """Return the size of the main queue."""
        return self._db.llen(self._processing_q_key)

    def empty(self):
        """Return True if the queue is empty, including work being done, False otherwise.

        False does not necessarily mean that there is work available to work on right now,
        """
        return self._main_qsize() == 0 and self._processing_qsize() == 0

# TODO: implement this
#    def check_expired_leases(self):
#        """Return to the work queueReturn True if the queue is empty, False otherwise."""
#        # Processing list should not be _too_ long since it is approximately as long
#        # as the number of active and recently active workers.
#        processing = self._db.lrange(self._processing_q_key, 0, -1)
#        for item in processing:
#          # If the lease key is not present for an item (it expired or was
#          # never created because the client crashed before creating it)
#          # then move the item back to the main queue so others can work on it.
#          if not self._lease_exists(item):
#            TODO: transactionally move the key from processing queue to
#            to main queue, while detecting if a new lease is created
#            or if either queue is modified.

    def _itemkey(self, item):
        """Returns a string that uniquely identifies an item (bytes)."""
        return hashlib.sha224(item).hexdigest()

    def _lease_exists(self, item):
        """True if a lease on 'item' exists."""
        return self._db.exists(self._lease_key_prefix + self._itemkey(item))

    def lease(self, lease_secs=60, block=True, timeout=None):
        """Begin working on an item the work queue.
        Lease the item for lease_secs.  After that time, other
        workers may consider this client to have crashed or stalled
        and pick up the item instead.
        If optional args block is true and timeout is None (the default), block
        if necessary until an item is available."""
        if block:
            item = self._db.brpoplpush(self._main_q_key, self._processing_q_key, timeout=timeout)
        else:
            item = self._db.rpoplpush(self._main_q_key, self._processing_q_key)
        if item:
            # Record that we (this session id) are working on a key.  Expire that
            # note after the lease timeout.
            # Note: if we crash at this line of the program, then GC will see no lease
            # for this item a later return it to the main queue.
            itemkey = self._itemkey(item)
            self._db.setex(self._lease_key_prefix + itemkey, lease_secs, self._session)
        return item

    def complete(self, value):
        """Complete working on the item with 'value'.
        If the lease expired, the item may not have completed, and some
        other worker may have picked it up.  There is no indication
        of what happened."""
        self._db.lrem(self._processing_q_key, 0, value)
        # If we crash here, then the GC code will try to move the value, but it will
        # not be here, which is fine.  So this does not need to be a transaction.
        itemkey = self._itemkey(value)
        self._db.delete(self._lease_key_prefix + itemkey, self._session)

# TODO: add functions to clean up all keys associated with "name" when
# processing is complete.
# TODO: add a function to add an item to the queue.  Atomically
# check if the queue is empty and if so fail to add the item
# since other workers might think work is done and be in the process
# of exiting.
# TODO(etune): move to my own github for hosting, e.g. github.com/erictune/rediswq-py and
# make it so it can be pip installed by anyone (see
# http://stackoverflow.com/questions/8247605/configuring-so-that-pip-install-can-work-from-github)
# TODO(etune): finish code to GC expired leases, and call periodically
#  e.g. each time lease times out.




def parse_args():
    starting_parser = argparse.ArgumentParser(description="This tool automatically fixes, cuts and indexes DNA sequences in FASTA format. \nRequired software: bowtie, bowtie2, samtools.")
    starting_parser.add_argument("-i", "--input", required=True,
                                 help="Reference DNA FASTA file;")
    starting_parser.add_argument("-s", "--size", type=float, default=3.6,
                                 help="(Optional) Sequence chunk size in billion characters, 3.6 ;")
    starting_parser.add_argument("-t", "--threads", type=int, default=int(multiprocessing.cpu_count()),
                                 help="(Optional) Number of CPU cores to speed up the chunks processing, maximal by default;")
    starting_parser.add_argument("-l", "--large_index", default=False,
                                 help="(Optional) Force usage of a 'large' index, disables sequence splitting;")
    starting_parser.add_argument("-o", "--output", required=True,
                                 help="Output directory. Must not exist and shall be created.")
    return starting_parser.parse_args()


def parse_namespace():
    namespace = parse_args()
    namespace.input = str(os.path.abspath(namespace.input))
    is_path_exists(namespace.output)
    namespace.output = ends_with_slash(os.path.abspath(namespace.output))
    for path in [namespace.input, namespace.output]:
        if ' ' in path:
            print("The path must not contain spaces: " + path + "\n Exiting...")
            sys.exit(2)
    return str(os.path.abspath(namespace.input)), int(float(namespace.size) * 10 ** 9), int(namespace.threads), namespace.large_index, namespace.output


def var_to_file(var_to_write, file_to_write):
    file = open(file_to_write, 'w')
    file.write(var_to_write)
    var_to_write = None
    file.close()


def get_time():
    from datetime import datetime
    now = datetime.now()
    output_list = []
    for time_unit in [now.year, now.month, now.day, now.hour, now.minute, now.second]:
        time_unit = str(time_unit)
        if len(time_unit) < 2:
            time_unit = '0' + time_unit
        output_list.append(time_unit)
    return '-'.join(output_list)


def ends_with_slash(string):
    if string.endswith("/"):
        return string
    else:
        return str(string + "/")


def dump_sampledata(jsons_list):
    output_buffer = ""
    for i in jsons_list:
        output_buffer += i["sampledata"]["sample_name"] + "\t" + i["sampledata"]["sample_path"] + "\n"
    output_file = ends_with_slash(jsons_list[0]["output"]) + "sampledata_" + jsons_list[0]["mask"] + "/" + hostNameString + "_" + get_time() + ".sampledata"
    is_path_exists(ends_with_slash(jsons_list[0]["output"]) + "sampledata_" + jsons_list[0]["mask"])
    var_to_file(output_buffer, output_file)
    return output_file


def external_route(input_direction_list):
    process = subprocess.Popen(input_direction_list, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    if error:
        print(error)
    return output.decode("utf-8")


if __name__ == '__main__':
    # The script part is based on: https://raw.githubusercontent.com/kubernetes/website/master/docs/tasks/job/fine-parallel-processing-work-queue/worker.py
    # Uncomment next two lines if you do not have Kube-DNS working.
    # import os
    # host = os.getenv("REDIS_SERVICE_HOST")
    hostNameString = subprocess.getoutput("hostname")
    scriptDir = ends_with_slash('/'.join(os.path.abspath(sys.argv[0]).split('/')[:-1]))

    q = RedisWQ(name="bwt_filtering_pipeline_docker", host="redis")
    print("Worker with sessionID: " + q.sessionID())
    print("Initial queue state: empty=" + str(q.empty()))
    sampledata_queue_list = []
    while not q.empty():
        item = q.lease(lease_secs=10, block=True, timeout=2)
        if item is not None:
            itemstr = item.decode("utf=8")
            try:
                json_single_queue = json.loads(itemstr)
                # Note that all entries must have the same refdata
                # json example: {"filter": "", "coverage": "", "sampledata": {"sample_name": "", "sample_path": ""}, "mask": "", "threads": "", "output": ""}
                sampledata_queue_list.append(json_single_queue)
                if len(sampledata_queue_list) == int(json_single_queue["threads"]):
                    print("Loaded full queue on:", hostNameString)
                    sampleDataFileName = dump_sampledata(sampledata_queue_list)
                    external_route(["python3", scriptDir + "pipeline_wrapper.py", "-f", json_single_queue["filter"], "-c", json_single_queue["coverage"], "-m", json_single_queue["mask"], "-o", json_single_queue["output"]])
                    sampledata_queue_list = []
            except ValueError:
                print("Cannot parse JSON:", itemstr)
            q.complete(item)
        else:
            print("Waiting for work")
    print("Queue empty, exiting")
