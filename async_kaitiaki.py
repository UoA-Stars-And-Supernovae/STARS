"""
██   ██  █████  ██ ████████ ██  █████  ██   ██ ██ 
██  ██  ██   ██ ██    ██    ██ ██   ██ ██  ██  ██ 
█████   ███████ ██    ██    ██ ███████ █████   ██ 
██  ██  ██   ██ ██    ██    ██ ██   ██ ██  ██  ██ 
██   ██ ██   ██ ██    ██    ██ ██   ██ ██   ██ ██ 

                Kaitiaki [kaɪ.ti.a.ki]
                   by Sean Richards

Kaitiaki ("watcher" in Te Reo Māori) is a gatekeeper program
to monitor the STARS code and control the models it produces.
"""

import numpy as np
import os
import sys
import subprocess
import select
import time
from tqdm import tqdm
from decimal import Decimal
from queue import Queue
import imf

# All features except IMF currently completely unused
n    = 2    # Number of stars to sample
q    = 0.5  # Mass ratio: m2 = q * m1
IMF  = imf.chabrier # IMF to select
logP = np.arange(4, 4.2, 0.2)

threshold = 0.01

DEBUG_MODE = False
Queue_Size = 5 # The queue of the last few statuses to monitor
               # if the last QUEUE_SIZE statuses are *all* SNAFU,
               # then we assume fail state and kill it.

save_plotfiles = False

Z_array_BPASS = ['zem5', 'zem4', 'z001', 'z002',
                 'z003', 'z004', 'z005', 'z006',
                 'z008', 'z010', 'z014', 'z020',
                 'z030', 'z040']

Z_array = [1e-5, 1e-4, 0.001, 0.002,
           0.003, 0.004, 0.005, 0.006,
           0.008, 0.010, 0.014, 0.020,
           0.030, 0.040]

Z = 0.020
Z_fmt = 'z020'

if len(sys.argv > 0):
    i = sys.argv[1]
else:
    raise ValueError(f"Must pass a mass index to compute -- e.g., {sys.argv[0]} 4 to compute the 4th mass index.")

Kaitiaki_Header = __doc__

def modify_datafile(line, cursor, new_value):
    with open("data", 'r+') as f:
        contents = f.readlines()

        L = len(new_value)

        newline = contents[line][:cursor] + str(new_value) + contents[line][cursor+L:]

        contents[line] = newline

        f.seek(0)
        f.write("".join(contents))
        f.truncate()

        f.close()

def Python_to_FORTRAN(number):
    return '%.2E' % Decimal(str(number))

def FORTRAN_to_Python(number):
    parts = number.strip().split("E")
    exp = int(parts[1])
    base = float(parts[0])

    return base * 10 ** exp

def respect_queue(queue, item):
    if queue.full():
        queue.get()
        queue.put(item)
    else:
        queue.put(item)

def bash_command(command, timeout=5*60):
    cmd = command.split(" ")

    cplobj = None

    try:
        cplobj = subprocess.run(cmd, timeout=timeout, capture_output=True)
    except subprocess.TimeoutExpired:
        reason = "timeout"
        stdout = b""
        stderr = b""

    if cplobj != None:
        # Finished!

        stdout = cplobj.stdout
        stderr = cplobj.stderr
        reason = "finished"

    return stdout.decode('utf-8').strip(), stderr.decode('utf-8').strip(), reason

def debug(msgtype, message, critical=False):
    if msgtype == 'warning':
        header = "\033[93m\033[1m[WARNING]\033[0m "
    elif msgtype == 'error':
        header = "\033[91m\033[1m[ERROR]\033[0m "
    elif msgtype == 'info':
        header = "\033[96m\033[1m[INFO]\033[0m "
    elif msgtype == 'status':
        header = "\033[92m\033[1m[STATUS]\033[0m "

    if DEBUG_MODE or critical:
        print(header + str(message))

def write_file(fname, contents, prepend=False):
    if prepend:
        with open(fname, 'r') as original:
            data = original.read()

        with open(fname, 'w') as modified:
            modified.write(contents + data)
    else:
        with open(fname, 'w+') as f:
            f.seek(0)
            f.write(contents)
            f.truncate()

def get_hydrogen_abundance(hydrogen_queue=None):
    out, _, _ = bash_command('tail -1 hydrogen')
    if out.strip() == '':
        # first run, just pretend its infinite
        ret = np.infty
    else:
        if 'E' in out:
            ret = FORTRAN_to_Python(out.strip())
        else:
            ret = float(out.strip())

    if hydrogen_queue is not None:
        respect_queue(hydrogen_queue, ret)

    return ret

def communicate(message):
    bash_command('rm instatus')
    write_file('instatus', message)

def prepend_eof_to_file(f_from, f_to, n_lines):
    bash_command(f'echo -e "$(tail -{n_lines} {f_from})\n$(cat {f_to})" > {f_to}')

def save_plot(Z_fmt):
    if save_plotfiles:
        with open(f'plot', 'r+') as fplot:
            with open(f'model_plots/{Z_fmt}/plot', 'a') as fsave:
                fsave.writelines(fplot.readlines())
                fsave.truncate()

def get_mass():
    out, _, _ = bash_command('tail -399 modout')
    model = out.split("\n")[0]
    parameters = model.split()

    mass = FORTRAN_to_Python(parameters[0])

    out, _, _ = bash_command('cat out')
    parts = out.split("\n\n")
    mass = -np.infty

    for part in parts[::-1]:
        if 'dt/age/MH/MHe' in part:
            mass = part.strip().split("\n")[1]
            mass = mass.strip().split()[0]
            break

    return float(mass)

def mark_iteration():
    if DEBUG_MODE:
        print()
        print("  \033[94m\033[1m---[NEW ITERATION]---\033[0m  ")

IMF_Masses = 10**np.arange(0, 3.1, 0.1)
m1 = IMF_Masses[i]

if DEBUG_MODE:
    print(Kaitiaki_Header)

mark_iteration()

m1 = IMF_Masses[i]

if os.path.exists(f'models/modout_{Z_fmt}_{round(np.log10(m1), 2)}'):
    debug("info", f"Z={Z}, m={round(m1, 2)} has already been computed. Skipping!")
    continue

# Restore the original datafile
bash_command('rm data')
bash_command('cp data.bak data')

modify_datafile(0,  31, '1') # ITH
modify_datafile(0,  35, '0') # IX
modify_datafile(0,  39, '0') # IY
modify_datafile(0,  43, '0') # IZ
modify_datafile(1,  31, '0') # ML prescription
modify_datafile(17, 19, '0.00E-01') # RE-ML factor
modify_datafile(17, 1,  '1.00E+03') # RCD
modify_datafile(17, 10, '0.12E+00') # Convective overshooting

# Create the output directory if it doesn't exist already

bash_command('rm status')
bash_command('rm hydrogen')
bash_command('touch status')
bash_command('touch hydrogen')

# wipe the nucmodin file
bash_command('rm nucmodin')
bash_command('touch nucmodin')

# Move the COtables around
if os.path.exists(f"../COTables/COtables_{Z_fmt}"):
    bash_command(f"cp ../COTables/COtables_{Z_fmt} dat/COtables")

    # Convert it to a fortran-formatted exponent
    num = '%.2E' % Decimal(str(Z))

    # And now insert it into the datafile.
    modify_datafile(16, 1, num)
else:
    # We are missing an opacity table... skip this one
    debug("warning", f"Missing opacity table for {Z_fmt}, and cannot evaluate the grid there. Skipping.", True)
    continue

file = 'modin.zams'

# if os.path.exists('modin.last'):
#     file = 'modin.last'

bash_command('rm modin')
bash_command(f'cp {file} modin')

# Run STARS
k = 0

# check if we reached the target mass:
modin, _, _ = bash_command('head -1 modin')
mass = float(modin.split(" ")[0])

debug('status', f'Current mass: {mass}')

# Turn on RE-ML
modify_datafile(1, 31, '9')
debug("info", f"Registered target mass of {Python_to_FORTRAN(m1)}")
modify_datafile(17, 19, Python_to_FORTRAN(m1))

# Turn off thermal energy generation
modify_datafile(0, 31, '0')

last_mass = Queue(maxsize=Queue_Size)
respect_queue(last_mass, mass)

_, _, reason = bash_command('./run_bs')

if reason == 'finished':
    # STARS terminated successfully

    # Step 1: check if we reached the target mass:
    mass = get_mass()
    respect_queue(last_mass, mass)

    if list(last_mass.queue) == [mass for _ in range(Queue_Size)]:
        # not converging
        debug('warning', 'Mass is not converging -- aborting.')
        break

    debug('status', f'Current mass: {mass}')
else:
    # Something fucked up
    debug("warning", f"run_bs died -- error code \"{reason}\".")
    break

if mass > m1 + threshold:
    # FUCK, we overshot. Using RE-ML we should never reach this but
    # just in case.
    debug("warning", (f"Overshot target ZAMS mass! "
                      f"(target: {m1}, "
                      f"max: {m1+threshold}, "
                      f"achieved: {mass})"))

# ok now we should have a ZAMS star.
out, err, _ = bash_command(f"tail -399 modout")

write_file(f'models/modout_{Z_fmt}_{round(np.log10(m1), 2)}', out)
write_file(f'modin.last', out)

debug('status', f'Finished Z={Z}, target mass={m1}, actual mass={mass}', True)