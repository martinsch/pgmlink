#!/usr/bin/python

import h5py
import numpy as np
import glob


def create_dataset(shape,
                   filenames,
                   f,
                   dtype=np.uint32,
                   label_path='segmentation/labels'):
    filenames = sorted(filenames)
    shape = (len(filenames),) + shape
    dset = np.empty(shape, dtype=dtype)
    for idx, fn in enumerate(filenames):
        print "Adding label image @ timestep %-2d (%s)" % (idx, fn)
        with h5py.File(fn, 'r') as f_in:
            dset[idx,...,0] = f_in[label_path][...]
    f.create_dataset('exported_data', data=dset)

def append_tracking_group(filenames,
                          f,
                          tracking_group='tracking',
                          ts_name='%04d'):
    filenames = sorted(filenames)
    tg = f.create_group('tracking')
    for idx, fn in enumerate(filenames):
        print "Adding tracking group @ timestep %-2d (%s)" % (idx, fn)
        with h5py.File(fn, 'r') as f_in:
            tg_at = tg.create_group(ts_name % idx)
            for t_ds in f_in[tracking_group].keys():
                tg_at.create_dataset(t_ds, data = f_in[tracking_group][t_ds])

    


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--label-directory', '-d', required=True, help='Path to directory containing h5 files with label image.')
    parser.add_argument('--tracking-directory', '-D', required=True, help='Path to directory containing h5 files with tracking group.')
    parser.add_argument('--label-path', '-l', default='segmentation/labels', help='Path to label image inside a h5 file. [default=%(default)s]')
    parser.add_argument('--out', '-o', required=True, help='Path to result h5 file.')
    parser.add_argument('--dtype', '-t', default='uint32', help='numpy dtype of the result data set. [default=%(default)s]')
    parser.add_argument('--tracking-group', '-g', default='tracking', help='Path to tracking group. [default=%(default)s]')
    parser.add_argument('--ts-name', '-n', default='%04d', help='Tracking subgroup pattern. [default=%(default)s]')
    parser.add_argument('--shape', '-s', required=True, help='Shape of data at one time frame.')
    args = parser.parse_args()
    filenames1 = sorted(glob.glob('%s/*.h5' % args.label_directory))
    filenames2 = sorted(glob.glob('%s/*.h5' % args.tracking_directory))
    label_path = args.label_path
    shape = tuple(int(i) for i in args.shape.replace('(','').replace(')','').split(','))
    with h5py.File(args.out, 'w') as f:
          create_dataset(shape,
                         filenames1,
                         f,
                         label_path=label_path,
                         dtype=np.dtype(args.dtype))
          append_tracking_group(filenames2,
                                f,
                                args.tracking_group,
                                args.ts_name)
