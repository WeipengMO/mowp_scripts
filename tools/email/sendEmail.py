#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2020-01-19 15:25:08
LastEditTime : 2022-04-01 19:51:38
'''


from email.mime.text import MIMEText
from email.header import Header
from smtplib import SMTP_SSL
import os
import sys
import datetime
import yaml


def load_config():
    script_path = os.path.split(os.path.realpath(__file__))[0]
    # load config.ymal
    with open(script_path+'/config.yaml', 'r') as f:
        config = yaml.safe_load(f)

    return config



class MAIL():
    def __init__(self, message=None, file=None, cmd=None, return_code=0, stdout=None, stderr=None):
        self.message = message
        self.file = file
        self.cmd = cmd
        self.return_code = return_code
        self.stdout = stdout
        self.stderr = stderr

        # host info
        self.hostname = os.popen('echo $HOSTNAME').read().rstrip()
        self.workdir = os.popen('echo $PWD').read().rstrip()
        self.jobid = ''
        self.jobname = ''

        # html format
        self.html_header = '<font face="Courier New, Courier, monospace" size=3><pre>'
        self.html_footer = '</pre></font>'


    def get_pbs_job_info(self,):
        self.jobid = os.popen('echo $PBS_JOBID').read().rstrip()
        self.jobname = os.popen('echo $PBS_JOBNAME').read().rstrip()
        
    
    def get_content(self,):
        title = f'【IPF】 Message'

        body = self.message if self.message is not None else 'No message'
        e = datetime.datetime.today()
        signature = f'<br>Cheers,<br>{e:%H:%M:%S %Y-%m-%d (%a)}'

        content = self.html_header+body+'<br>'+signature+self.html_footer

        return title, content
    
    def send(self,):
        config = load_config()
        #qq邮箱smtp服务器
        host_server = config['host_server']
        #sender_qq为发件人的qq号码
        sender_qq = config['sender_qq']
        #pwd为qq邮箱的授权码
        pwd = config['pwd']
        #发件人的邮箱
        sender_qq_mail = config['sender_qq_mail']
        #收件人邮箱
        receiver = config['receiver']

        mail_title, mail_content = self.get_content()

        if self.file:
            from email.mime.application import MIMEApplication
            from email.mime.multipart import MIMEMultipart

            msg = MIMEMultipart()
            # 构建正文
            part_text = MIMEText(mail_content, "html", 'utf-8')
            msg.attach(part_text)  # 把正文加到邮件体里面去

            # 构建邮件附件
            part_attach1 = MIMEApplication(open(self.file, 'rb').read())  # 打开附件
            filename = os.path.basename(self.file)
            part_attach1.add_header('Content-Disposition', 'attachment', filename=filename)  # 为附件命名
            msg.attach(part_attach1)
        else:
            msg = MIMEText(mail_content, "html", 'utf-8') # 添加附件
        
        msg["Subject"] = Header(mail_title, 'utf-8')
        msg["From"] = sender_qq_mail
        msg["To"] = receiver

        #ssl登录
        smtp = SMTP_SSL(host_server)
        smtp.set_debuglevel(0)
        smtp.ehlo(host_server)
        smtp.login(sender_qq, pwd)

        smtp.sendmail(sender_qq_mail, receiver, msg.as_string())
        smtp.quit()


class MAIL_CMD(MAIL):
    def get_content(self,):
        # Title
        if not self.return_code:
            title = f'【IPF】Job at {self.hostname} Successfully completed'
        else:
            title = f'【IPF】Job at {self.hostname} Failed'

        # Content
        body = f'''<b>Job path:</b>
{self.workdir}

************************************************************

<b>Your job looked like:</b>
{self.cmd}

************************************************************

'''

        if self.stdout is not None and self.stdout != '':
            html_stdout = f'''<b>Stdout looked like:</b>
{self.stdout}

************************************************************

'''
        else:
            html_stdout = ''
    
        if self.stderr is not None and self.stderr != '':
            html_stderr = f'''<b>Stderr looked like:</b>
{self.stderr}

************************************************************

'''
        else:
            html_stderr = ''

        e = datetime.datetime.today()
        if not self.return_code:
            signature = f'<br>Successfully completed :)<br>{e:%H:%M:%S %Y-%m-%d (%a)}'
        else:
            signature = f'<br>Job Error :(<br>{e:%H:%M:%S %Y-%m-%d (%a)}'

        content = self.html_header+body+html_stdout+html_stderr+signature+self.html_footer

        return title, content


class MAIL_ATTACHMENT(MAIL):
    def get_content(self,):
        if not os.path.exists(self.file):
            print('No such file or directory')
            exit(1)
        
        if os.path.getsize(self.file) > 1048576:
            print('File is too large')
            exit(1)
        
        title = f'【IPF】File Attachment'

        filepath = os.path.realpath(self.file)
        e = datetime.datetime.today()
        body = f'''<b>File path:</b>
{filepath}

The file is attached to this email.
'''
        signature = f'<br>Cheers,<br>{e:%H:%M:%S %Y-%m-%d (%a)}'
        content = self.html_header+body+signature+self.html_footer

        return title, content


if __name__ == '__main__':
    if len(sys.argv) > 1:
        message = ' '.join(sys.argv[1:])
    else:
        message=None
    mail = MAIL(message=message)
    mail.send()