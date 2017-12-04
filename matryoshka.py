#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import getopt
import sys
import re
import os
import multiprocessing
import subprocess
import logging
import paramiko


def usage():
    print("\nUsage: " + sys.argv[0] + " -i <file> -o <dir> -w -t <int> -n <str> -s <str>" + "\n\n" +
          "-i/--input <file> \tA tab-delimited table without a header containing all required variables" + "\n" +
          "-o/--output <dir> \tOutput directory" + "\n" +
          "-w/--wait <bool> \t(Optional) If specified, program will wait for completion of all tasks, always enabled for local use" + "\n" +
          "-t/--threads <int> \t(Optional) Number of threads to create" + "\n" +
          "-n/--nodes <str> \t(Optional) Comma-divided list of nodes or text file with authentication data one per line in format: \"node:username:password:port\"" + "\n" +
          "-s/--string <str> \tA constant string containing zero-based indexes of required table columns supplied by the \"\$\" symbols" + "\n\n" +
          "The example table: " + "\n" +
          "A1\tB1\nA2\tB2" + "\n\n" +
          "The example constant string: " + "\n" +
          "\"program constant_arg1 \$0 constant_arg2 \$1\"" + "\n\n" +
          "The example bash analogue: " + "\n" +
          "program constant_arg1 A1 constant_arg2 B1\nprogram constant_arg1 A2 constant_arg2 B2" + "\n\n")
    sys.exit(2)


def main():
    opts = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:wt:n:s:", ["help", "input=", "output=", "wait=", "threads=", "nodes=", "string="])
    except getopt.GetoptError as arg_err:
        print(str(arg_err))
        usage()
    i = None
    o = None
    w = False
    t = None
    n = None
    s = None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
        elif opt in ("-i", "--input"):
            i = str(arg)
        elif opt in ("-o", "--output"):
            o = str(arg)
        elif opt in ("-w", "--wait"):
            w = True
        elif opt in ("-t", "--threads"):
            try:
                t = int(arg)
            except ValueError:
                print("Incorrect threads number!")
                usage()
        elif opt in ("-n", "--nodes"):
            n = str(arg)
        elif opt in ("-s", "--string"):
            s = str(arg)
    if not any(var is None for var in [i, o, s]):
        return i, o, w, t, n, s
    print("The parameters are not yet specified!")
    usage()


def file_to_str(file):
    file_parsed = open(file, 'rU').read()
    return file_parsed


def list_to_file(header, list_to_write, file_to_write):
    header += "".join(str(i) for i in list_to_write if i if i is not None)
    file = open(file_to_write, 'w')
    file.write(header)
    file.close()


def var_to_file(var_to_write, file_to_write):
    file = open(file_to_write, 'w')
    file.write(var_to_write)
    file.close()


def is_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        print(exception)


def ends_with_slash(string):
    if string.endswith("/"):
        return string
    else:
        return str(string + "/")


def stdin_assembler(const_string, var_strings_list):
    const_string_indexes = sorted(list(map(lambda x: int(x), list(filter(None, re.findall('\$([0-9]*)', const_string))))))
    output_strings_list = []
    for var_string in list(filter(None, var_strings_list)):
        var_string_splitted = list(filter(None, var_string.replace('\r', '\n').replace('\n', '').split('\t')))
        output_string = str(const_string) + " "
        for const_string_index in const_string_indexes:
            output_string = output_string.replace("$" + str(const_string_index) + " ", str(var_string_splitted[const_string_index]) + " ")
        if output_string.endswith(' '):
            output_string = output_string[:-1]
        output_strings_list.append(output_string)
    return output_strings_list


def list_chop(input_list, number_of_chunks):
    chunks_length = int(len(input_list) / number_of_chunks)
    if sys.version_info >= (3, 0):
        chunks_list = [input_list[i:i + chunks_length] for i in range(0, (number_of_chunks - 1) * chunks_length, chunks_length)]
    else:
        chunks_list = [input_list[i:i + chunks_length] for i in xrange(0, (number_of_chunks - 1) * chunks_length, chunks_length)]
    chunks_list.append(input_list[(number_of_chunks - 1) * int(chunks_length):])
    return chunks_list


def make_single_core_dummy(dummy_commands_list, dummy_mask):
    dummy_body = str("#!/usr/bin/env python3)" + "\n" +
                     "# -*- coding: utf-8 -*- " + "\n" +
                     "" + "\n" +
                     "import os" + "\n" +
                     "import sys" + "\n" +
                     "import subprocess" + "\n" +
                     "import logging" + "\n" +
                     "" + "\n" +
                     "" + "\n" +
                     "def external_route(input_direction, output_direction):" + "\n" +
                     "\tlogging.info(\"Processing: \" + input_direction + \" with logging to \" + output_direction)" + "\n" +
                     "\tcmd = input_direction.split(\" \")" + "\n" +
                     "\tprocess = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)" + "\n" +
                     "\t(output, error) = process.communicate()" + "\n" +
                     "\tfile_append(output.decode(\"utf-8\"), output_direction)" + "\n" +
                     "\tprocess.wait()" + "\n" +
                     "" + "\n" +
                     "" + "\n" +
                     "def file_append(string, file_to_append):" + "\n" +
                     "\tfile = open(file_to_append, \"a+\")" + "\n" +
                     "\tfile.write(string)" + "\n" +
                     "\tfile.close()" + "\n" +
                     "" + "\n" +
                     "" + "\n" +
                     "commands = sorted(list(filter(None, " + str(dummy_commands_list) + ")))" + "\n" +
                     "mask = \"" + str(dummy_mask) + "\"" + "\n" +
                     "" + "\n" +
                     "logging.basicConfig(format=u\'%(levelname)-8s [%(asctime)s] %(message)s\', level=logging.DEBUG, filename=str(mask + \"_master.log\"))" + "\n" +
                     "logging.info(\"Launching command: \" + \" \".join(str(i) for i in sys.argv))" + "\n" +
                     "logging.info(\"Main process ID: \" + str(os.getpid()))" + "\n" +
                     "" + "\n" +
                     "iteration = 0" + "\n" +
                     "for command in commands:" + "\n" +
                     "\ttry:" + "\n" +
                     "\t\texternal_route(command, mask + \"_stdout_\" + str(iteration) + \".log\")" + "\n" +
                     "\t\tlogging.info(\"Successfully processed: \" + command + \" iteration \" + str(iteration) + \" with logging to \" + mask + \"_stdout_\" + str(iteration) + \".log\")" + "\n" +
                     "\t\titeration += 1" + "\n" +
                     "\texcept Exception as exception:" + "\n" +
                     "\t\tlogging.critical(\"COMMAND \" + command + \" HAS CRASHED AT ITERATION \" + str(iteration) + \" WITH LOG FILE \" + mask + \"_stdout_\" + str(iteration) + \".log AND ERROR CODE \" + str(exception))" + "\n" +
                     "" + "\n" +
                     "logging.info(\"COMPLETED\")" + "\n"
                     )
    var_to_file(dummy_body, dummy_mask + ".py")
    return dummy_mask + ".py"


def make_multi_core_dummy(dummy_commands_list, dummy_mask, cores_number):
    dummy_body = str("#!/usr/bin/env python3)" + "\n" +
                     "# -*- coding: utf-8 -*- " + "\n" +
                     "" + "\n" +
                     "import os" + "\n" +
                     "import sys" + "\n" +
                     "import subprocess" + "\n" +
                     "import logging" + "\n" +
                     "import multiprocessing" + "\n" +
                     "" + "\n" +
                     "" + "\n" +
                     "def external_route(input_direction, output_direction):" + "\n" +
                     "\tlogging.info(\"Processing: \" + input_direction + \" with logging to \" + output_direction)" + "\n" +
                     "\tcmd = input_direction.split(\" \")" + "\n" +
                     "\tprocess = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)" + "\n" +
                     "\t(output, error) = process.communicate()" + "\n" +
                     "\tfile_append(output.decode(\"utf-8\"), output_direction)" + "\n" +
                     "\tprocess.wait()" + "\n" +
                     "" + "\n" +
                     "" + "\n" +
                     "def file_append(string, file_to_append):" + "\n" +
                     "\tfile = open(file_to_append, \"a+\")" + "\n" +
                     "\tfile.write(string)" + "\n" +
                     "\tfile.close()" + "\n" +
                     "" + "\n" +
                     "" + "\n" +
                     "def multi_core_queue(function):" + "\n" +
                     "\tpool = multiprocessing.Pool(cores)" + "\n" +
                     "\tpool.map(function, commands)" + "\n" +
                     "\tpool.close()" + "\n" +
                     "\tpool.join()" + "\n" +
                     "" + "\n" +
                     "" + "\n" +
                     "def external_wrapper(command):" + "\n" +
                     "\tlogging.info(\"Created process with ID: \" + str(os.getpid()))" + "\n" +
                     "\titeration = [i for i, x in enumerate(commands) if x == command][0]" + "\n" +
                     "\ttry:" + "\n" +
                     "\t\texternal_route(command, mask + \"_stdout_\" + str(iteration) + \".log\")" + "\n" +
                     "\t\tlogging.info(\"Successfully processed: \" + command + \" iteration \" + str(iteration) + \" with logging to \" + mask + \"_stdout_\" + str(iteration) + \".log\")" + "\n" +
                     "\texcept Exception as exception:" + "\n" +
                     "\t\tlogging.critical(\"COMMAND \" + command + \" HAS CRASHED AT ITERATION \" + str(iteration) + \" WITH LOG FILE \" + mask + \"_stdout_\" + str(iteration) + \".log AND ERROR CODE \" + str(exception))" + "\n" +
                     "" + "\n" +
                     "" + "\n" +
                     "commands = sorted(list(filter(None, " + str(dummy_commands_list) + ")))" + "\n" +
                     "mask = \"" + str(dummy_mask) + "\"" + "\n" +
                     "cores = " + str(cores_number) + "\n" +
                     "" + "\n" +
                     "logging.basicConfig(format=u\'%(levelname)-8s [%(asctime)s] %(message)s\', level=logging.DEBUG, filename=str(mask + \"_master.log\"))" + "\n" +
                     "logging.info(\"Launching command: \" + \" \".join(str(i) for i in sys.argv))" + "\n" +
                     "logging.info(\"Main process ID: \" + str(os.getpid()))" + "\n" +
                     "" + "\n" +
                     "multi_core_queue(external_wrapper)" + "\n" +
                     "" + "\n" +
                     "logging.info(\"COMPLETED\")" + "\n"
                     )
    var_to_file(dummy_body, dummy_mask + ".py")
    return dummy_mask + ".py"


def dummy_multiprocessing_check(commands, mask):
    if inputThreadsNumber is not None and inputThreadsNumber > 1:
        # if inputThreadsNumber > int(multiprocessing.cpu_count()):
        #     logging.critical("The threads number is too large! Available CPU cores: " + str(int(multiprocessing.cpu_count())))
        #     logging.info("Exiting...")
        #     sys.exit(2)
        dummy_name = make_multi_core_dummy(commands, mask, inputThreadsNumber)
        logging.info("The multi-threading wrappers were created!")
    else:
        dummy_name = make_single_core_dummy(commands, mask)
        logging.info("The single-threading wrappers were created!")
    return dummy_name


def local_use():
    commands_total_list = stdin_assembler(inputString, file_to_str(inputFile).split("\n"))
    dummy = dummy_multiprocessing_check(commands_total_list, outputDir + "dummy_local")
    dummy_launch_command = "nohup nice -n 19 python3 " + dummy
    return dummy_launch_command


def remote_use(nodes_2d_array):
    commands_total_list = stdin_assembler(inputString, file_to_str(inputFile).split("\n"))
    commands_splitted_list = list_chop(commands_total_list, len(nodes_2d_array))
    dummies_launch_commands = []
    node_index = 0
    for node_list in nodes_2d_array:
        dummy = dummy_multiprocessing_check(commands_splitted_list[node_index], outputDir + "dummy_" + node_list[0])
        dummies_launch_commands.append([node_list, "nohup nice -n 19 python3 " + dummy])
        node_index += 1
    return dummies_launch_commands


def get_available_nodes():
    alive_nodes = []
    for node_list in input2nodes_lists(nodesInputList):
        try:
            check_node(node_list)
            alive_nodes.append(node_list)
        except:
            continue
    if len(alive_nodes) == 0:
        logging.critical("No alive nodes! Exiting...")
        sys.exit(2)
    return alive_nodes


def input2nodes_lists(string_or_file):
    if os.path.isfile(string_or_file):
        nodes = sorted(list(filter(None, file_to_str(string_or_file).split('\n'))), reverse=True)  # our newest nodes are faster
    else:
        nodes = sorted(list(filter(None, string_or_file.split(","))), reverse=True)
    nodes_2d_array = []
    for node in nodes:
        node_list = node2list(node)
        if node_list is not None:
            if len(node_list) > 0:
                nodes_2d_array.append(node_list)
    return nodes_2d_array


def node2list(node):
    node_list = []
    iteration = 0
    while iteration <= 3:
        try:
            node_list.append(node.split(":")[iteration])
        except IndexError:
            node_list.append('')
        iteration += 1
    if len(node_list[3]) == 0:
        node_list[3] = 22
    else:
        try:
            node_list[3] = int(node_list[3])
        except ValueError:
            node_list[3] = 22
    if node_list[0]:
        return node_list


def check_node(destination):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        client.connect(hostname=destination[0], username=destination[1], password=destination[2], port=destination[3])
    except paramiko.ssh_exception.AuthenticationException:
        try:
            client.connect(hostname=destination[0], username=destination[1], port=destination[3])
        except paramiko.ssh_exception.AuthenticationException:
            try:
                client.connect(hostname=destination[0], port=destination[3])
            except paramiko.ssh_exception.AuthenticationException:
                raise
    client.close()


def execute_via_ssh(destination, remote_command):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        client.connect(hostname=destination[0], username=destination[1], password=destination[2], port=destination[3])
    except paramiko.ssh_exception.AuthenticationException:
        try:
            client.connect(hostname=destination[0], username=destination[1], port=destination[3])
        except paramiko.ssh_exception.AuthenticationException:
            try:
                client.connect(hostname=destination[0], port=destination[3])
            except paramiko.ssh_exception.AuthenticationException:
                raise
    stdin, stdout, stderr = client.exec_command("cd " + outputDir[:-1] + "; " + remote_command)
    if waitBoolean:
        data = stdout.read() + stderr.read()
        return data.decode("utf-8")


def external_route(input_direction, output_direction):
    logging.info("Processing: " + input_direction + " with logging to " + output_direction)
    cmd = input_direction.split(" ")
    process = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output, error) = process.communicate()
    process.wait()
    file_append(output.decode("utf-8"), output_direction)


def file_append(string, file_to_append):
    file = open(file_to_append, 'a+')
    file.write(string)
    file.close()


def filename_only(string):
    if len(str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))) == 0:
        return ".".join(string.rsplit("/", 1)[-1].split("."))
    return str(".".join(string.rsplit("/", 1)[-1].split(".")[:-1]))


def multi_core_queue(function, queue):
    pool = multiprocessing.Pool()
    pool.map(function, queue)
    pool.close()
    pool.join()


def launch_facility(commands_list):
    logging.info("Created process with ID: " + str(os.getpid()) + " for processing: " + commands_list[1] + " via SSH on " + "".join(str(i) for i in commands_list[0]))
    try:
        execute_via_ssh(commands_list[0], commands_list[1])
        logging.info("Succesfully processed: " + commands_list[1] + " via SSH on " + "".join(str(i) for i in commands_list[0]))
    except Exception as exception:
        logging.critical("PROCESS WITH ID " + str(os.getpid()) + " HAS CRASHED WHILE PROCESSING: " + commands_list[1] + " VIA SSH on " + "".join(str(i) for i in commands_list[0]) + " WITH ERROR CODE " + exception)


########################################
inputFile, outputDir, waitBoolean, inputThreadsNumber, nodesInputList, inputString = main()

if sys.version_info < (3, 3):
    print("This program requires Python 3.3+ interpreter!\nExiting...")
    sys.exit(2)

is_path_exists(outputDir)
outputDir = ends_with_slash(os.path.abspath(outputDir))
logging.basicConfig(format=u'%(levelname)-8s [%(asctime)s] %(message)s', level=logging.DEBUG, filename=str(outputDir + filename_only(sys.argv[0]) + "_master.log"))
logging.info("Main process ID: " + str(os.getpid()))

if nodesInputList is None:
    logging.info("Using local strategy!")
    finalCommand = local_use()
    os.chdir(outputDir)
    external_route(finalCommand, outputDir + filename_only(sys.argv[0]) + ".log")
else:
    logging.info("Using remote strategy!")
    multi_core_queue(launch_facility, remote_use(get_available_nodes()))

logging.info("COMPLETED")
