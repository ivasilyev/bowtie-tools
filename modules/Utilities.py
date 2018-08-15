#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import logging
import subprocess


class Utilities:
    @staticmethod
    def ends_with_slash(string):
        return (string + "/", string)[string.endswith("/")]

    @staticmethod
    def remove_empty_values(input_list):
        output_list = []
        if input_list:
            for i in input_list:
                if i:
                    try:
                        if len(i) > 0:
                            output_list.append(i)
                    except TypeError:
                        continue
        return output_list

    @staticmethod
    def log_and_raise(message):
        logging.critical(message)
        print(message)
        raise ValueError(message)

    @staticmethod
    def filename_only(string: str):
        return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))

    @staticmethod
    def load_string(file: str):
        with open(file=file, mode="r", encoding="utf-8") as f:
            s = f.read()
        return s

    @staticmethod
    def load_list(file: str):
        import re
        return Utilities.remove_empty_values([re.sub("[\r\n]+", "", i) for i in Utilities.load_string(file).split("\n")])

    @staticmethod
    def load_2d_array(file: str):
        return [i.split("\t") for i in Utilities.load_list(file)]

    @staticmethod
    def dump_string(string: str, file: str):
        with open(file=file, mode="w", encoding="utf-8") as f:
            f.write(string)

    @staticmethod
    def dump_list(lst: list, file: str):
        Utilities.dump_string(string="\n".join(lst) + "\n", file=file)

    @staticmethod
    def dump_2d_array(array: list, file: str):
        Utilities.dump_list(lst=["\t".join(i) for i in array], file=file)

    @staticmethod
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

    @staticmethod
    def external_route(input_direction: list, output_direction: str):
        # cmd = input_direction.split(" ")
        process = subprocess.Popen(input_direction, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (output, error) = process.communicate()
        process.wait()
        if error:
            print(error)
        if not output_direction:
            print(output.decode("utf-8"))
        else:
            Utilities.dump_string(string=output.decode("utf-8"), file=output_direction)

    @staticmethod
    def batch_remove(*args):
        import os
        lst = [i for i in args if os.path.isfile(i) or os.path.isdir(i)]
        if "/" not in lst:  # I know it is useless
            s = " ".join(lst)
            subprocess.getoutput("rm -rf {}".format(s))
            print("Removed items: {}".format(s))

    @staticmethod
    def single_core_queue(func, queue: list):
        return [func(i) for i in queue]

    @staticmethod
    def multi_core_queue(func, queue: list, processes: int = int(subprocess.getoutput("nproc").strip())):
        import multiprocessing
        pool = multiprocessing.Pool(processes=processes)
        output = pool.map_async(func, queue)
        pool.close()
        pool.join()
        return output

    @staticmethod
    def extracting_wrapper(func, *args):
        return func(args)
