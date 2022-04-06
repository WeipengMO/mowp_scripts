# pbs-torque

This profile configures Snakemake to run on the [Torque Scheduler](http://www.adaptivecomputing.com/products/open-source/torque/).

## Setup

### Deploy profile

To deploy this profile, copy this folder to `~/.config/snakemake`

    mkdir -p ~/.config/snakemake
    cd ~/.config/snakemake

Then, you can run Snakemake with

    snakemake --profile pbs ...
