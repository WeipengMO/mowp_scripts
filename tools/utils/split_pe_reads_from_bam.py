import pysam
import click


@click.command()
@click.argument('infile', type=click.Path(exists=True))
@click.argument('fwd', type=click.Path())
@click.argument('rev', type=click.Path())
def main(infile: str, fwd: str, rev: str):
    inbam = pysam.AlignmentFile(infile, 'rb')
    if infile == fwd or infile == rev:
        raise ValueError('Input and output files must be different')
        
    out_fwd = pysam.AlignmentFile(fwd, 'wb', template=inbam)
    out_rev = pysam.AlignmentFile(rev, 'wb', template=inbam)

    for read in inbam:
        if read.is_read1:
            if read.is_forward:
                out_fwd.write(read)
            else:
                out_rev.write(read)
        else:
            if read.is_forward:
                out_rev.write(read)
            else:
                out_fwd.write(read)

    inbam.close()
    out_fwd.close()
    out_rev.close()

    pysam.index(fwd, '-@ 30')
    pysam.index(rev, '-@ 30')


if __name__ == '__main__':
    main()