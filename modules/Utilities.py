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
    def split_lines(string: str):
        import re
        return Utilities.remove_empty_values([i.strip() for i in re.sub("[\r\n]+", "\n", string).split("\n")])

    @staticmethod
    def load_list(file: str):
        return Utilities.split_lines(Utilities.load_string(file))

    @staticmethod
    def string_to_2d_array(string: str):
        return Utilities.remove_empty_values([[j.strip() for j in i.split("\t")] for i in Utilities.split_lines(string)])

    @staticmethod
    def _2d_array_to_dicts_list(arr: list, names: list):
        if len(arr[0]) != len(names):
            raise ValueError("Cannot parse dictionary: keys number is not equal!")
        out = []
        for row_list in arr:
            counter = 0
            while counter < len(row_list):
                out.append({names[counter]: row_list[counter]})
                counter += 1
        return [[j for j in i] for i in arr]

    @staticmethod
    def load_2d_array(file: str):
        return Utilities.string_to_2d_array(Utilities.load_string(file))

    @staticmethod
    def load_dicts_list(file: str, names: list):
        arr = Utilities.load_2d_array(file)
        if len(arr[0]) != len(names):
            raise ValueError("Cannot parse dictionary: keys number is not equal!")

    @staticmethod
    def dump_string(string: str, file: str):
        with open(file=file, mode="w", encoding="utf-8") as f:
            f.write(string)

    @staticmethod
    def dump_list(lst: list, file: str):
        Utilities.dump_string(string="\n".join([str(i) for i in lst]) + "\n", file=file)

    @staticmethod
    def dump_2d_array(array: list, file: str):
        Utilities.dump_list(lst=["\t".join([str(j) for j in i]) for i in array], file=file)

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
            logging.critical("""
Error report:
{}
""".format(error))
        exit_code = process.returncode
        stdout = output.decode("utf-8")
        if not output_direction:
            logging.debug("""
Completed command arguments: '{a}'; exit code: {b}; STDOUT details: 
{c}
""".format(a=input_direction, b=exit_code, c=stdout))
        else:
            logging.debug("Completed command arguments: '{a}'; exit code: {b}; log file: '{c}'").format(a=input_direction, b=exit_code, c=output_direction)
            Utilities.dump_string(string=stdout, file=output_direction)

    @staticmethod
    def batch_remove(*args):
        import os
        items = [i for i in args if os.path.isfile(i) or (os.path.isdir(i) and i != "/")]  # I know it is useless
        if len(items) > 0:
            s = " ".join(items)
            subprocess.getoutput("rm -rf {}".format(s))
            logging.info("Removed items: {}".format(s))

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
