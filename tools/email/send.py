#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2022-03-04 14:52:38
LastEditTime : 2022-04-06 12:59:09
LastEditors  : windz
FilePath     : send.py
'''


import sys
sys.path.append('./')
import sendEmail
import subprocess
import click


def run_cmd(cmd, verbose=True):
    '''It runs a command and returns the output
    
    Parameters
    ----------
    cmd
        the command to run
    verbose, optional
        If True, print the command to the screen before running it.
    
    '''

    run_cmd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = run_cmd.communicate()
    return_code = run_cmd.returncode

    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')
    if verbose:
        if stdout != '':
            print(stdout)
        if stderr != '':
            print(stderr, file=sys.stderr)

    if len(stdout) > 10000:
        stdout = '<i>Output exceeds the size limit.</i><br>' + stdout[-10000: ]

    if len(stderr) > 10000:
        stderr = '<i>Output exceeds the size limit.</i><br>' + stderr[-10000: ]
    
    return stdout, stderr, return_code


###################
# cli1
###################
@click.group()
def cli1():
    pass

@cli1.command(context_settings=dict(ignore_unknown_options=True))
@click.argument('args', nargs=-1, type=click.UNPROCESSED)
@click.option('--verbose / --no-verbose', default=True)
def cmd(args, verbose):
    """Send CLI result"""
    args = ' '.join(args)
    stdout, stderr, return_code = run_cmd(args, verbose)

    # send email
    mail = sendEmail.MAIL_CMD(cmd=args, return_code=return_code, stdout=stdout, stderr=stderr)
    mail.send()


###################
# cli2
###################
@click.group()
def cli2():
    pass

@cli2.command()
@click.argument('filename', nargs=1, type=click.Path(exists=True))
def file(filename):
    """Send file attachment"""
    mail = sendEmail.MAIL_ATTACHMENT(file=filename)
    mail.send()


###################
# cli3
###################
@click.group()
def cli3():
    pass

@cli3.command()
def qsub():
    """Send qsub result"""
    jobid = subprocess.Popen('echo $PBS_JOBID', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0].decode('utf-8').rstrip()
    jobname = subprocess.Popen('echo $PBS_JOBNAME', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0].decode('utf-8').rstrip()

    if jobid != '':
        cmd = f'tail -n 100 /var/spool/torque/spool/{jobid}.OU'
    else:
        print('[warning] No Jobid found!')
        exit(0)

    stdout, stderr, return_code = run_cmd(cmd, verbose=False)

    # send email
    mail = sendEmail.MAIL_CMD(cmd=f'{jobname} ({jobid})', return_code=return_code, stdout=stdout, stderr=stderr)
    mail.send()


###################
# cli4
###################
@click.group()
def cli4():
    pass

@cli4.command()
@click.option('-l', '--log_file', required=True, type=str)
@click.option('--return_code', required=False, default=0)
def log(log_file, return_code):
    """Send snakemake log"""
    cmd = f'cat {log_file}'

    stdout, stderr, _ = run_cmd(cmd, verbose=False)
        
    mail = sendEmail.MAIL_CMD(cmd='snakemake', return_code=return_code, stdout=stdout, stderr=stderr)
    mail.send()
    

send = click.CommandCollection(sources=[cli1, cli2, cli3, cli4])

if __name__ == '__main__':
    send()