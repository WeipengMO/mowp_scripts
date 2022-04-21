# Send email with QQ mailbox

This is a tools for sending all kinds of logs to you QQ mailbox.

```
Usage: send.py [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  cmd   Send CLI result
  file  Send file attachment
  log   Send snakemake log
  qsub  Send qsub result
```

## Features

- `cmd`: send the results of command line

- `file`: send the file as attachment

- `log`: send snakemake log, add the below fields to your `snakemake` file:

    ```
    onsuccess:
        shell("send.py log -l {log}")

    onerror:
        shell("send.py log -l {log} --return_code 1")
    ```

- `qsub`: add below line to the end of the qsub script, and it will send the qsub log when the job is completed.

    `send.py qsub`


## Configs

Here you need to do two things: 

1. Open the mailbox SMTP function 

2. Get the authorization code

3. Write the authorization code to the configuration file [config.yaml](config.yaml)

For more detail, please refer to:
https://blog.csdn.net/MATLAB_matlab/article/details/106240424