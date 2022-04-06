#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2022-04-06 13:25:46
LastEditTime : 2022-04-06 14:26:13
LastEditors  : windz
FilePath     : pbs.py
'''

#!python

#!/usr/bin/env python3
import os
import subprocess
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
# jobscript = '/public/home/mowp/test/singularity/.snakemake/tmp.yjfn0i2u/snakejob.NAME.0.sh'
job_properties = read_job_properties(jobscript)

j = '-j oe'
env = '-V'
jobname = '-N ' + job_properties['rule']

# resources
if "threads" in job_properties:
    ppn = "ppn=" + str(job_properties["threads"])

if "resources" in job_properties:
    resources = job_properties["resources"]

    if "nodes" in resources:
        nodes="nodes=" + str(resources["nodes"])
    else:
        nodes="nodes=1"

    if "mem" in resources:
        mem="mem=" + str(resources["mem"])
    else:
        mem="mem=10g"

    if "walltime" in resources:
        walltime="walltime=" + str(resources["walltime"])
    else:
        walltime="walltime=900:00:00"

resourceparams = f'-l "{nodes}:{ppn},{mem},{walltime}"'

if "cluster" in job_properties:
    cluster = job_properties["cluster"]
    if "error" in cluster:
        os.makedirs(os.path.dirname(cluster["error"]), exist_ok=True)
        se = "-e " + cluster["error"]
    else:
        if not os.path.exists('log'):
            os.makedirs('log', exist_ok=True)
        se = "-e log"

    if "output" in cluster:
        os.makedirs(os.path.dirname(cluster["output"]), exist_ok=True)
        so = "-o " + cluster["output"]
    else:
        if not os.path.exists('log'):
            os.makedirs('log', exist_ok=True)
        so = "-o log"

cmd = f'qsub {env} {jobname} {resourceparams} {j} {so} {se} {jobscript}'

try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

res = res.stdout.decode()
print(res.strip())